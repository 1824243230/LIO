#include "lio/plane_fit.h"
#include "lio/deskew.h"

#include "lio/super_lio.h"

#include <sys/resource.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>


using namespace BASIC;

namespace LI2Sup{

inline bool compute_error(
  const std::array<double, 4>& abcd, const V3& point, 
  const float length, scalar& error)
{
  error = abcd[0] * point[0] + abcd[1] * point[1] + abcd[2] * point[2] + abcd[3];
  return length > 81 * error * error;
}


void SuperLIO::init(){
  ros::NodeHandle nh;
  nh.param("/lio/imu_init/require_stationary", require_stationary_init_, false);
  nh.param("/lio/imu_init/enable_zupt", enable_zupt_, false);
  ROS_INFO_STREAM("IMU stationary-window initialization: " << require_stationary_init_);
  ROS_INFO_STREAM("IMU zero-velocity update: " << enable_zupt_);
  ivox_.reset(new OctVoxMapType(OctVoxMapType::Options{g_ivox_resolution, g_ivox_capacity}));
  kf_.reset(new ESKF());
  InitGeometry();
  data_wrapper_->setESKF(kf_);
  
  scan_undistort_full_.reset(new PointCloudType());
  ds_undistort_.reset(new PointCloudType());
  world_pc_.reset(new PointCloudType());
  ds_world_.reset(new PointCloudType());

  if(g_save_map){
    point_map_.reset(new PointCloudType());
  }
  
  points_world_v3_.reserve(21000);
  abcd_vec_.resize(20000);
  effect_knn_idxs_.resize(20000);
  voxel_grid_fliter_.setLeafSize(g_voxel_fliter_size);

  state_fn_ = &SuperLIO::stateWaitKFInit;

  LOG(INFO) << GREEN << " ---> [SuperLIO]: initialized." << RESET;
}


void SuperLIO::stateWaitKFInit()
{
  if (kf_init()) {
    state_fn_ = &SuperLIO::stateWaitMapInit;
    LOG(INFO) << GREEN << " ---> [SuperLIO]: KF init done" << RESET;
  }
}

void SuperLIO::stateWaitMapInit()
{
  if (map_init()) {
    kf_->init_ = true;
    state_fn_ = &SuperLIO::stateProcess;
    LOG(INFO) << GREEN << " ---> [SuperLIO]: Map init done" << RESET;
  }
}

void SuperLIO::process(){
  if(!data_wrapper_->sync_measure(measures_)){
    return;
  }
  if (enable_zupt_) {
    // The synchronizer lends an IMU sample beyond the scan end. Only use
    // measurements already available at this state time; add() deduplicates
    // the borrowed sample when it appears in the next scan.
    for (const auto& imu : measures_.imu)
      if (imu.secs <= measures_.lidar.end_time) zupt_imu_window_.add(imu);
  }
  (this->*state_fn_)();
}


bool SuperLIO::kf_init(){
  if (measures_.imu.empty()) return false;
  for (const auto& imu : measures_.imu) {
    imu_initialization_.add(imu);
    if (require_stationary_init_) stationary_imu_window_.add(imu);
  }

  /// 100 Hz for 1 second.
  if(imu_initialization_.count < 50){
    return false;
  }

  V3 mean_gyro = imu_initialization_.mean_gyro;
  V3 mean_acce = imu_initialization_.mean_acc;
  if (require_stationary_init_) {
    const auto estimate = stationary_imu_window_.estimate();
    if (!estimate.ready) {
      // In a run without a quiet interval, keep waiting instead of fabricating
      // a gyro bias from motion. The optional mode will produce no odometry.
      ROS_WARN_THROTTLE(5.0, "Waiting for a 1 s stationary IMU window before LIO initialization");
      return false;
    }
    mean_gyro = estimate.mean_gyro;
    mean_acce = estimate.mean_acc;
    ROS_INFO_STREAM("Stationary IMU initialization: " << estimate.count
                    << " samples in " << estimate.duration << " s; gyro bias="
                    << mean_gyro.transpose() << "; gyro std="
                    << estimate.std_gyro.transpose() << "; acc std="
                    << estimate.std_acc.transpose());
  }

  // 零加速度均值无法确定重力方向，继续等待数据，避免归一化产生 NaN。
  if (!mean_acce.allFinite() || !mean_gyro.allFinite() || mean_acce.norm() < 1e-6 ||
      !std::isfinite(g_gravity_norm) || g_gravity_norm <= 0) return false;
  V3 gravity = - mean_acce * g_gravity_norm / mean_acce.norm();
  V3 ref_gravity(0, 0, - g_gravity_norm);
  M3 init_rot = Quat::FromTwoVectors(gravity, ref_gravity).toRotationMatrix();
  V3 n = init_rot.col(0);
  double yaw = atan2(n(1), n(0));

  M3 R_yaw_inv = Eigen::AngleAxis<scalar>(-yaw, V3::UnitZ()).toRotationMatrix(); 

  // init_rot represents the IMU orientation after gravity alignment (level orientation).
  // Perform LiDAR leveling correction, then transform the orientation into the robot frame.
  M3 rot = g_lidar_robo_yaw * R_yaw_inv * init_rot;  

  ESKF::Options options;
  options.smooth_motion_ = smooth_motion_;
  options.gyro_var_ = g_imu_ng;
  options.acce_var_ = g_imu_na;
  options.bias_gyro_var_ = g_imu_nbg;
  options.bias_acce_var_ = g_imu_nba;
  options.num_iterations_ = g_kf_max_iterations;
  options.quit_eps_ = g_kf_quit_eps;

  float imu_scale = g_gravity_norm / mean_acce.norm();
  kf_->SetInitialConditions(options, mean_gyro, V3::Zero(), imu_scale, ref_gravity);
  // 初始平移沿用机器人参考原点约定；时间由 IMU 锚点统一设置。
  kf_->SetPoseAtImu(SE3(SO3(rot), g_odom_robo.t_), measures_.imu.back());
  sys_init_pose_ = kf_->GetSE3();
  return true;
}


bool SuperLIO::map_init(){
  // 初始化后的数帧也会运动。必须把状态和点云推进到各自扫描末端，
  // 不能固定初始位姿拼图后仅改 last_obs_time，留下旧状态/新时间的组合。
  if (!Propagation_Undistort()) return false;
  const auto pose = kf_->GetSE3();
  points_world_v3_.resize(scan_undistort_full_->size());
  for (size_t i=0;i<scan_undistort_full_->size();++i) {
    const auto& p = scan_undistort_full_->points[i];
    points_world_v3_[i] = pose*V3(p.x,p.y,p.z);
  }
  ivox_->insert(points_world_v3_);
  UpdateBumpMap();
  ++frame_num_;
  if(frame_num_ > 3){
    g_flg_map_init = false;
    return true;
  }
  return false;
}


// 每帧：IMU传播/去畸变 → 采样 → 迭代观测 → 成功后入图。失败预测不得污染历史地图。
void SuperLIO::stateProcess(){
  observation_valid_ = false;
  bool propagated = false;
  if (g_time_eva)
    time_record_.Evaluate([&](){ propagated = Propagation_Undistort(); }, "Undistort");
  else
    propagated = Propagation_Undistort();
  // 未达到末端的状态不能配准本帧，也不能发布带有错误时刻的点云。
  // 已成功积分的前缀仍保留，后续扫描可继续使用有效 IMU 端点。
  if (!propagated) return;
  ++frame_num_;
  if(g_time_eva){
    time_record_.Evaluate([this]() { DownSample(); }, "DownSample");
    time_record_.Evaluate([this]() { Observe(); }, "Observe");
    if (enable_zupt_) MaybeApplyZupt();
    if (observation_valid_) time_record_.Evaluate([this]() { UpdateMap(); }, "UpdateMap");
  }else{
    DownSample();
    Observe();
    if (enable_zupt_) MaybeApplyZupt();
    if (observation_valid_) UpdateMap();
  }
  PublishGeometry();
  Output();
  if (observation_valid_) caceData();
}

bool SuperLIO::MaybeApplyZupt() {
  if (!enable_zupt_ || !kf_) return false;
  const auto estimate = zupt_imu_window_.estimate();
  const auto velocity = kf_->GetNavState().v;
  const bool stationary = estimate.ready && estimate.mean_acc.allFinite() &&
    std::isfinite(g_gravity_norm) &&
    std::abs(double(estimate.mean_acc.norm())-g_gravity_norm) <= 0.5 &&
    velocity.allFinite() && velocity.norm() < 0.25;
  const bool applied = stationary && kf_->UpdateZeroVelocity();
  if (applied) {
    if (!zupt_active_)
      ROS_INFO_STREAM("ZUPT active at t=" << kf_->GetTime()
                      << ": stationary IMU window, velocity=" << velocity.norm() << " m/s");
    ++zupt_streak_;
    ++zupt_total_;
    if (zupt_streak_ % 100 == 0)
      ROS_INFO_STREAM("ZUPT active: " << zupt_streak_ << " frames in this interval, "
                      << zupt_total_ << " total");
  } else if (zupt_active_) {
    ROS_INFO_STREAM("ZUPT ended at t=" << kf_->GetTime() << ": " << zupt_streak_
                    << " frames in this interval, " << zupt_total_ << " total");
    zupt_streak_ = 0;
  }
  zupt_active_ = applied;
  // observation_valid_ remains the LiDAR result. ZUPT must never authorize
  // insertion of a scan whose LiDAR update failed.
  return applied;
}


bool SuperLIO::Propagation_Undistort(){
  motion_covariances_.clear();
  scan_undistort_full_->clear();
  const auto& scan = measures_.lidar;
  if (!scan.pc || scan.pc->empty() || !std::isfinite(scan.start_time) ||
      !std::isfinite(scan.end_time) || scan.end_time < scan.start_time ||
      scan.end_time <= kf_->GetTime()) return false;
  std::vector<motion::Knot> knots;
  const auto append_knot = [&]() {
    const auto state = kf_->GetDynamicState();
    motion::Knot knot;
    knot.time = state.time; knot.R = state.R.cast<double>(); knot.p = state.p.cast<double>();
    if (knots.empty()) knot.covariance = kf_->GetCov().cast<double>();
    else {
      knot.transition = kf_->GetTransition().cast<double>();
      knot.noise = kf_->GetProcessNoise().cast<double>();
      knot.covariance = knot.transition*knots.back().covariance*knot.transition.transpose()+knot.noise;
    }
    knots.push_back(knot);
  };
  propagate_states_.clear();
  propagate_states_.emplace_back(kf_->GetDynamicState());
  if (uncertainty_motion_) append_knot();
  kf_->SetObsTime(measures_.lidar.end_time);
  for (auto &imu : measures_.imu) {
    if (kf_->Predict(imu) && kf_->GetTime() > propagate_states_.back().time) {
      propagate_states_.emplace_back(kf_->GetDynamicState());
      if (uncertainty_motion_) append_knot();
    }
  }

  // Predict 可能因坏输入/数值失败拒绝一步。不能把最后一次成功传播时刻
  // 当成扫描末端，否则去畸变会夹取未覆盖的点并将其当作有效观测。
  if (propagate_states_.size() < 2 || kf_->GetTime() != scan.end_time) {
    ROS_WARN_THROTTLE(1.0, "Incomplete IMU propagation; skipping LiDAR scan");
    return false;
  }

  const M3 TLI_R = g_lidar_imu.R_;
  const V3 TLI_t = g_lidar_imu.t_;
  const double start_time = measures_.lidar.start_time;
  auto& raw_pc = measures_.lidar.pc;

  std::size_t ptsize = raw_pc->points.size();
  scan_undistort_full_->resize(ptsize);
  std::vector<motion::M3> point_covariances;
  if (uncertainty_motion_) {
    motion::relativeCovariances(knots);
    point_covariances.resize(ptsize);
  }

  tbb::parallel_for(
  tbb::blocked_range<size_t>(0, ptsize),
  [&](const tbb::blocked_range<size_t>& r) {
    for (size_t idx = r.begin(); idx < r.end(); ++idx) {  
      auto& pt_full = scan_undistort_full_->points[idx];
      const auto& pt = raw_pc->points[idx];
      pt_full.intensity = pt.intensity;
      const double query_time = start_time + pt.offset_time;
      const V3 raw(pt.x, pt.y, pt.z);
      const V3 eigen_point = deskewPoint(TLI_R * raw + TLI_t, query_time, propagate_states_, smooth_motion_);
      pt_full.x = eigen_point[0];
      pt_full.y = eigen_point[1];
      pt_full.z = eigen_point[2];
      if (uncertainty_motion_)
        point_covariances[idx] = motion::pointCovariance(eigen_point.cast<double>(), query_time, knots);
    }
  });
  for (const auto& p : *scan_undistort_full_) {
    if (!p.getVector3fMap().allFinite()) {
      scan_undistort_full_->clear();
      return false;
    }
  }
  if (uncertainty_motion_) {
    motion_covariances_.reserve(ptsize);
    for (size_t i=0;i<ptsize;++i) {
      const auto& p = scan_undistort_full_->points[i];
      auto result = motion_covariances_.emplace(std::array<float,3>{p.x,p.y,p.z}, point_covariances[i]);
      // Coincident samples may have different times: sum is a conservative PSD bound.
      if (!result.second) result.first->second += point_covariances[i];
    }
  }
  return true;
}



void SuperLIO::AnalyzeSamplingGeometry() {
  if (!geometry_options_.enable_degeneracy || !geometry_options_.enable_informed_sampling) return;
  geometry::M6 G = geometry::M6::Zero();
  const auto pose = kf_->GetSE3();
  KNNHeapType neighbors;
  ivox_->reset_max_group();
  for (const auto& p : *ds_undistort_) {
    const V3 body(p.x,p.y,p.z), world = pose*body;
    if (!world.allFinite()) continue;
    neighbors.reset(); ivox_->getTopK(world, neighbors);
    std::array<double,4> plane;
    scalar residual;
    if (neighbors.count < 4 || !calc_plane_coeff(neighbors.count, neighbors.points_, plane) ||
        !compute_error(plane, world, body.norm(), residual)) continue;
    const auto J = geometry::poseJacobian(body.cast<double>(), pose.R_.cast<double>(),
                                          geometry::V3(plane[0],plane[1],plane[2]));
    double variance = 0.001;
    if (uncertainty_motion_) {
      const auto it = motion_covariances_.find({p.x,p.y,p.z});
      if (it == motion_covariances_.end() || !it->second.allFinite()) continue;
      const Eigen::Vector3d nb = pose.R_.cast<double>().transpose()*J.tail<3>();
      variance += std::max(0.0, nb.dot(it->second*nb));
    }
    if (!std::isfinite(variance) || variance <= 0) continue;
    if (geometry_options_.plane_huber_delta > 0) {
      const double robust_weight = geometry::huberWeight(residual, geometry_options_.plane_huber_delta);
      if (robust_weight <= 0) continue;
      G += J * (robust_weight/variance) * J.transpose();
    } else {
      G += J * J.transpose()/variance;
    }
  }
  degeneracy_ = geometry_analyzer_.analyze(G, geometry_options_, true);
  sampling_geometry_analyzed_ = true;
}


struct ThreadACC{
  M6d HTVH = M6d::Zero();
  V6d HTVr = V6d::Zero();
  ThreadACC(): HTVH(M6d::Zero()), HTVr(V6d::Zero()) {}
};


// plane_only 用于从同一先验重放平面观测，禁止再次补偿或递归回退；谱保护仍生效。
void SuperLIO::Observe(bool plane_only){
  observation_valid_ = false;
  size_t ptsize = ds_undistort_->size();
  effect_mask_.assign(ptsize, 0);
  effect_knn_mask_.assign(ptsize, 0);
  effect_knn_idxs_.resize(ptsize);
  abcd_vec_.resize(ptsize);
  geometry_stats_ = geometry::Diagnostics{};
  std::vector<geometry::V6> plane_J;
  std::vector<double> plane_error;
  std::vector<double> plane_precisions(ptsize, 1000.0);
  std::vector<motion::M3> sampled_covariances;
  if (uncertainty_motion_) {
    sampled_covariances.reserve(ptsize);
    for (const auto& p : *ds_undistort_) {
      const auto it = motion_covariances_.find({p.x,p.y,p.z});
      if (it == motion_covariances_.end() || !it->second.allFinite()) {
        ROS_ERROR_THROTTLE(1.0, "Missing or invalid motion covariance: rejecting scan update");
        return;
      }
      sampled_covariances.push_back(it->second);
    }
  }
  if (geometry_options_.enable_degeneracy) {
    plane_J.resize(ptsize, geometry::V6::Zero());
    plane_error.resize(ptsize, 0);
  }
  const ESKF prior = *kf_; // includes covariance, gravity, bias and forward propagation state
  const auto prior_pose = kf_->GetSE3();
  bool allow_bump = !plane_only && geometry_options_.enable_bump_measurement && bool(bump_map_);
  bool used_bump = false;
  bool geometry_failed = false; // A failed analysis disables bump for the rest of this frame.
  int observation_calls = 0;
  static std::vector<float> _lengths;
  points_body_v3_.resize(ptsize);
  _lengths.resize(ptsize);

  effect_knn_num_ = ptsize;
  std::iota(effect_knn_idxs_.begin(), effect_knn_idxs_.begin() + ptsize, 0);

  for(size_t i = 0; i < ptsize; ++i){
    const auto& point_body_pcl = ds_undistort_->points[i];
    points_body_v3_[i] = V3(point_body_pcl.x, point_body_pcl.y, point_body_pcl.z);
    _lengths[i] = points_body_v3_[i].norm();
  }

  ivox_->reset_max_group();
  int iter_num = 0;

  auto observation = [&, this](const ESKF::KFState &kf_state, M6 &HTVH, V6 &HTVr) {
    const SE3 pose = kf_state.pose;
    const bool need_converge = kf_state.need_converge;
    const M3d R_transpose = (pose.R_.transpose()).cast<double>();

    tbb::enumerable_thread_specific<ThreadACC> tls_acc;

    tbb::parallel_for(
      tbb::blocked_range<size_t>(0, effect_knn_num_),
      [&](const tbb::blocked_range<size_t>& r) {
        KNNHeapType top_K;
        auto& local_acc = tls_acc.local();
        for (size_t r_s = r.begin(); r_s < r.end(); ++r_s) {
          int idx = effect_knn_idxs_[r_s];
          V3& point_body = points_body_v3_[idx];
          V3 point_world = pose * point_body;

          if(!need_converge){
            top_K.reset();
            ivox_->getTopK(point_world, top_K);
            if(top_K.count < 4){
              effect_mask_[idx] = false;
              effect_knn_mask_[idx] = false;
              continue;
            }
            effect_knn_mask_[idx] = true;
            effect_mask_[idx] = calc_plane_coeff(top_K.count, top_K.points_, abcd_vec_[idx]);
          }

          if(!effect_mask_[idx]) continue;

          auto& abcd = abcd_vec_[idx];
          scalar error;
          effect_mask_[idx] = compute_error(abcd, point_world, _lengths[idx], error);
          if(!effect_mask_[idx]) continue;
          
          {
            V3d normvec(abcd[0], abcd[1], abcd[2]);
            V3d nb = R_transpose * normvec;
            V3d point_body_d = point_body.cast<double>();
            V6d J;
            J.head<3>() = point_body_d.cross(nb);
            J.tail<3>() = normvec;
      
            if (geometry_options_.enable_degeneracy) {
              plane_J[idx] = J;
              plane_error[idx] = error;
            }
            double variance = 0.001;
            if (uncertainty_motion_)
              variance += std::max(0.0, nb.dot(sampled_covariances[idx]*nb));
            const double robust_weight = geometry_options_.plane_huber_delta > 0
                ? geometry::huberWeight(error, geometry_options_.plane_huber_delta) : 1.0;
            if (robust_weight <= 0 || !std::isfinite(variance) || variance <= 0) {
              effect_mask_[idx] = false;
              continue;
            }
            // The same effective precision must be removed if a bump replaces
            // this plane row later in the iteration.
            plane_precisions[idx] = robust_weight/variance;
            local_acc.HTVH += J * plane_precisions[idx] * J.transpose();
            local_acc.HTVr -= J * plane_precisions[idx] * error;
          }
        }
    });

    M6d sum_HTVH = M6d::Zero();
    V6d sum_HTVr = V6d::Zero();
    for(const auto& local_acc : tls_acc){
      sum_HTVH += local_acc.HTVH;
      sum_HTVr += local_acc.HTVr;
    }
    // 平面谱用于选择补偿方向；最终平面/曲面混合系统稍后单独计算谱可靠性。
    if (geometry_options_.enable_degeneracy) {
      using Clock = std::chrono::steady_clock;
      const auto analysis_start = Clock::now();
      degeneracy_ = geometry_analyzer_.analyze(sum_HTVH, geometry_options_, observation_calls == 0 && !sampling_geometry_analyzed_ && !plane_only);
      geometry_failed = geometry_failed || !degeneracy_.valid;
      if (geometry_failed) allow_bump = false;
      geometry_stats_ = geometry::Diagnostics{};
      geometry_stats_.analysis_ms = std::chrono::duration<double, std::milli>(Clock::now()-analysis_start).count();
      const auto bump_start = Clock::now();
      std::vector<geometry::Candidate> candidates;
      double plane_squared = 0;
      for (size_t i = 0; i < ptsize; ++i) {
        if (!effect_mask_[i]) continue;
        ++geometry_stats_.planes;
        plane_squared += plane_error[i]*plane_error[i];
        if (!allow_bump || !degeneracy_.valid || degeneracy_.status == geometry::Status::NORMAL) continue;
        geometry::Surface surface;
        if (!bump_map_->query((pose*points_body_v3_[i]).cast<double>(), surface)) {
          ++geometry_stats_.rejected_map; continue;
        }
        geometry::Candidate c; c.index = i;
        if (geometry::candidate(surface, points_body_v3_[i].cast<double>(), pose.R_.cast<double>(),
                                degeneracy_, geometry_options_, c, geometry_stats_)) {
          if (uncertainty_motion_) {
            const Eigen::Vector3d nb = R_transpose*c.J.tail<3>();
            c.variance += std::max(0.0, nb.dot(sampled_covariances[i]*nb));
          }
          candidates.push_back(c);
        }
      }
      geometry::selectWeakDirectionConstraints(candidates, degeneracy_, geometry_options_);
      // Strategy A: only final accepted indices replace their own plane rows.
      geometry::applyExclusiveBumps(plane_J, plane_error, effect_mask_, candidates, sum_HTVH, sum_HTVr, 1000.0, &plane_precisions);
      double bump_squared = 0;
      for (const auto& c : candidates) {
        plane_squared -= plane_error[c.index]*plane_error[c.index];
        bump_squared += c.residual*c.residual;
        geometry_stats_.avg_mid += c.mid;
        geometry_stats_.avg_gradient += c.gradient;
        geometry_stats_.avg_weak += c.weak_score;
        geometry_stats_.avg_covariance += c.variance;
      }
      geometry_stats_.accepted = candidates.size();
      geometry_stats_.planes -= candidates.size();
      if (!candidates.empty()) {
        used_bump = true;
        geometry_stats_.avg_mid /= candidates.size();
        geometry_stats_.avg_gradient /= candidates.size();
        geometry_stats_.avg_weak /= candidates.size();
        geometry_stats_.avg_covariance /= candidates.size();
        geometry_stats_.bump_rms = std::sqrt(bump_squared/candidates.size());
      }
      geometry_stats_.plane_rms = std::sqrt(std::max(0.0, plane_squared)/std::max(1, geometry_stats_.planes));
      geometry_stats_.bump_ms = std::chrono::duration<double, std::milli>(Clock::now()-bump_start).count();
    }
    // 必须在曲面替换之后加权，避免压制已经恢复的方向；平面回退也不能绕过保护。
    if (geometry_options_.enable_spectral_reliability &&
        !geometry::applySpectralReliability(sum_HTVH, sum_HTVr, geometry_options_, geometry_stats_)) {
      // ESKF rejects zero information and restores the complete propagated prior.
      sum_HTVH.setZero();
      sum_HTVr.setZero();
    }
    ++observation_calls;
    HTVH = sum_HTVH.cast<scalar>();
    HTVr = sum_HTVr.cast<scalar>();

    if(need_converge) return;

    int _effect_knn_num = 0;
    for(size_t i = 0; i < effect_knn_num_; ++i){
      int idx = effect_knn_idxs_[i];
      if(!effect_knn_mask_[idx]) continue;
      effect_knn_idxs_[_effect_knn_num] = idx;
      _effect_knn_num++;
    }

    // LOG(INFO) << "effect_knn_num_: " << effect_knn_num_ << ", _effect_knn_num: " << _effect_knn_num;
    effect_knn_num_ = _effect_knn_num;

    iter_num++;
  };
  const bool update_ok = kf_->UpdateObserve(observation, [&](const ESKF::STATE& dx) {
    return !used_bump || (dx.allFinite() &&
      dx.head<3>().norm() <= geometry_options_.bump_max_rotation_correction &&
      dx.segment<3>(3).norm() <= geometry_options_.bump_max_translation_correction);
  });
  bool abnormal = false;
  if (used_bump) {
    const auto pose = kf_->GetSE3();
    const double translation = (pose.t_-prior_pose.t_).norm();
    const double rotation = SO3(prior_pose.R_.transpose()*pose.R_).log_vee().norm();
    abnormal = !update_ok || !std::isfinite(translation) || !std::isfinite(rotation) || !kf_->GetCov().allFinite() ||
      translation > geometry_options_.bump_max_translation_correction ||
      rotation > geometry_options_.bump_max_rotation_correction;
  }
  const bool changed_sampling = original_sample_ && ds_undistort_ != original_sample_;
  const bool no_final_bump = geometry_stats_.accepted == 0;
  if (!plane_only && geometry_options_.enable_bump_measurement &&
      (abnormal || (no_final_bump && (changed_sampling || used_bump)))) {
    const auto attempted = geometry_stats_;
    // 恢复整帧先验与原采样后重做对应关系，不能混用增强采样的缓存或后验。
    *kf_ = prior;
    if (original_sample_) ds_undistort_ = original_sample_;
    // Rebuild points, masks, correspondences and all iterations from the original prior/sample.
    // plane_only prevents recursion and does not advance frame hysteresis a second time.
    Observe(true);
    geometry_stats_.candidates = attempted.candidates;
    geometry_stats_.rejected_map = attempted.rejected_map;
    geometry_stats_.rejected_mid = attempted.rejected_mid;
    geometry_stats_.rejected_gradient = attempted.rejected_gradient;
    geometry_stats_.rejected_pixel = attempted.rejected_pixel;
    geometry_stats_.rejected_weak = attempted.rejected_weak;
    geometry_stats_.rejected_residual = attempted.rejected_residual;
    geometry_stats_.bump_ms += attempted.bump_ms;
    geometry_stats_.analysis_ms += attempted.analysis_ms;
    geometry_stats_.fallback = true;
    return;
  }
  if (!plane_only && geometry_options_.enable_bump_measurement && no_final_bump &&
      (!degeneracy_.valid || degeneracy_.status != geometry::Status::NORMAL))
    geometry_stats_.fallback = true; // Already used original sampling and only plane rows.

  observation_valid_ = update_ok;
  if (!update_ok) ROS_WARN_THROTTLE(1.0, "LIO observation rejected; retaining IMU prior and skipping map insertion");
}


void SuperLIO::UpdateMap() {
  const size_t ptsize = ds_undistort_->size();
  if (ptsize == 0) return;
  
  last_pose_ = kf_->GetSE3();
  points_world_v3_.resize(ptsize);
  
  const auto R = last_pose_.R_;
  const auto t = last_pose_.t_;
  
  for (size_t i = 0; i < ptsize; ++i) {
    const auto& pt = points_body_v3_[i];
    points_world_v3_[i] = R * pt + t;
  }
  
  ivox_->insert(points_world_v3_);
  UpdateBumpMap();

}


void SuperLIO::Output(){
  auto state = kf_->GetNavState();
  data_wrapper_->pub_odom(state);  

  Eigen::Matrix4f transformation = Eigen::Matrix4f::Identity();
  transformation.block<3, 3>(0, 0) = state.R.R_.cast<float>();
  transformation.block<3, 1>(0, 3) = state.p.cast<float>();

  CloudPtr world_pc(new PointCloudType());
  
  if(g_visual_map){
    static int count = -1;
    count++;
    if(count % g_pub_step != 0){
      return;
    }
    count = 0;
    if(g_visual_dense){
      pcl::transformPointCloud(*scan_undistort_full_, *world_pc, transformation);
      data_wrapper_->pub_cloud_world(world_pc, state.timestamp);
    }else{
      pcl::transformPointCloud(*ds_undistort_, *world_pc, transformation);
      data_wrapper_->pub_cloud_world(world_pc, state.timestamp);
    }
  }
}

void SuperLIO::printTimeRecord(){
  if(!g_time_eva) return;
  time_record_.PrintAll();
}

} // namespace END.
