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
  (this->*state_fn_)();
}


bool SuperLIO::kf_init(){
  static int imu_cout = 0;
  static V3 mean_gyro = V3::Zero();
  static V3 mean_acce = V3::Zero();

  for(auto& imu: measures_.imu){
    imu_cout ++;
    mean_gyro += (imu.gyr - mean_gyro) / imu_cout;
    mean_acce += (imu.acc - mean_acce) / imu_cout;
  }

  /// 100 Hz for 1 second.
  if(imu_cout < 50){
    return false;
  }

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
  options.gyro_var_ = g_imu_ng;
  options.acce_var_ = g_imu_na;
  options.bias_gyro_var_ = g_imu_nbg;
  options.bias_acce_var_ = g_imu_nba;
  options.num_iterations_ = g_kf_max_iterations;
  options.quit_eps_ = g_kf_quit_eps;

  float imu_scale = g_gravity_norm / mean_acce.norm();
  kf_->SetInitialConditions(options, mean_gyro, V3::Zero(), imu_scale, ref_gravity);
  auto state = kf_->GetSysState();
  state.R = SO3(rot);
  state.p = g_odom_robo.t_;        // By default, the robot frame is used as the reference origin.
  state.timestamp = measures_.imu.back().secs;
  kf_->SetX(state);
  sys_init_pose_ = kf_->GetSE3();
  return true;
}


bool SuperLIO::map_init(){
  frame_num_++;

  std::size_t ptsize = measures_.lidar.pc->size();
  points_world_v3_.resize(ptsize);

  const SE3 transform = sys_init_pose_ * g_lidar_imu;

  tbb::parallel_for(
    tbb::blocked_range<size_t>(0, ptsize),
    [&](const tbb::blocked_range<size_t>& r) {
      for (size_t idx = r.begin(); idx < r.end(); ++idx) {
        auto& point_pcl = measures_.lidar.pc->points[idx];
        V3 point_body(point_pcl.x, point_pcl.y, point_pcl.z);
        points_world_v3_[idx] = transform * point_body;
      }
    }
  );

  ivox_->insert(points_world_v3_);
  if (bump_map_) {
    std::vector<geometry::V3> points;
    for (const auto& p : points_world_v3_) points.push_back(p.cast<double>());
    bump_map_->insert(points, transform.t_.cast<double>());
  }
  kf_->SetLastObsTime(measures_.lidar.end_time);

  if(frame_num_ > 3){
    g_flg_map_init = false;
    return true;
  }
  return false;
}


void SuperLIO::stateProcess(){
  frame_num_++;
  if(g_time_eva){
    time_record_.Evaluate([this](){Propagation_Undistort();}, "Undistort");
    time_record_.Evaluate([this]() { DownSample(); }, "DownSample");
    time_record_.Evaluate([this]() { Observe(); }, "Observe");
    if (observation_valid_) time_record_.Evaluate([this]() { UpdateMap(); }, "UpdateMap");
  }else{
    Propagation_Undistort();
    DownSample();
    Observe();
    if (observation_valid_) UpdateMap();
  }
  PublishGeometry();
  Output();
  if (observation_valid_) caceData();
}


void SuperLIO::caceData(){
  if(!g_save_map) return;
  auto state = kf_->GetNavState();
  Eigen::Matrix4f transformation = Eigen::Matrix4f::Identity();
  transformation.block<3, 3>(0, 0) = state.R.R_.cast<float>();
  transformation.block<3, 1>(0, 3) = state.p.cast<float>();

  if(g_if_filter){
    pcl::transformPointCloud(*ds_undistort_, *world_pc_, transformation);
  }else{
    pcl::transformPointCloud(*scan_undistort_full_, *world_pc_, transformation);
  }

  static int scan_wait_num = 0;
  if(!world_pc_->empty()){
    *point_map_ += *world_pc_;
    scan_wait_num++;
  }

  if(g_pcd_save_interval < 0) {
    scan_wait_num = 0;
    return;
  }

  static bool rm_PCD_dir = false;
  if(!rm_PCD_dir){
    rm_PCD_dir = true;
    std::string cmd = "rm -rf " + g_save_map_dir + "/PCD";
    [[maybe_unused]] int res;
    res = system(cmd.c_str());
    cmd = "mkdir -p " + g_save_map_dir + "/PCD";
    res = system(cmd.c_str());
  }

  if (point_map_->size() > 0 && scan_wait_num >= g_pcd_save_interval) {
    pcd_index_++;
    std::string map_name(std::string(g_save_map_dir + "/PCD/scans_") + std::to_string(pcd_index_) +
                               std::string(".pcd"));
    LOG(INFO) << GREEN << " ---> current scan saved to /PCD/scans_" << pcd_index_ << "  size:  " << point_map_->size() << RESET;
    pcl::io::savePCDFileBinary(map_name, *point_map_);
    point_map_->clear();
    scan_wait_num = 0;
  }
}


void SuperLIO::ProcessCaceMap(){
  namespace fs = std::filesystem;

  std::string pcd_folder = g_save_map_dir + "/PCD";
  std::string output_map_name = g_save_map_dir + "/" + g_map_name;

  LOG(INFO) << YELLOW << " ---> Merging PCD fragments in: " << pcd_folder << RESET;

  PointCloudType::Ptr merged_map(new PointCloudType());

  int count = 0;
  std::error_code directory_error;
  fs::directory_iterator fragment_it(pcd_folder, directory_error), end;
  if (directory_error) {
    LOG(WARNING) << "Cannot read PCD directory: " << directory_error.message();
    return;
  }
  for (; fragment_it != end; fragment_it.increment(directory_error)) {
    if (directory_error) break;
    const auto& entry = *fragment_it;
    if (entry.path().extension() == ".pcd" &&
      entry.path().filename().string().find("scans_") != std::string::npos) {
      PointCloudType::Ptr tmp_cloud(new PointCloudType());
      if (pcl::io::loadPCDFile<PointType>(entry.path().string(), *tmp_cloud) == 0) {
        *merged_map += *tmp_cloud;
        count++;
        // LOG(INFO) << GREEN << " ---> Merged: " << entry.path().filename().string() 
        //           << "   size: " << tmp_cloud->size() << RESET;
      } else {
        LOG(WARNING) << RED << " ---> Failed to load: " << entry.path().string() << RESET;
      }
    }
  }

  if (directory_error) {
    LOG(ERROR) << "Cannot iterate PCD directory: " << directory_error.message();
    return;
  }
  if (merged_map->empty()) {
    LOG(WARNING) << "No points in PCD fragments; skipping map save.";
    return;
  }
  LOG(INFO) << YELLOW << " ---> Total merged fragments: " << count << RESET;

  PointCloudType filtered_map;

  if(g_if_filter){
    LOG(INFO) << YELLOW << " ---> Downsampling merged map before final save..." << RESET;
    pcl::VoxelGrid<PointType> voxel_filter;
    voxel_filter.setLeafSize(g_map_ds_size, g_map_ds_size, g_map_ds_size);
    
    voxel_filter.setInputCloud(merged_map);
    voxel_filter.filter(filtered_map);
  }else{
    LOG(INFO) << YELLOW << " ---> Not Downsampling merged map before final save..." << RESET;
    filtered_map = *merged_map;
  }
  
  if (filtered_map.size() > 0) {
    filtered_map.width = filtered_map.size();
    filtered_map.height = 1;
    filtered_map.is_dense = false;
  }

  pcl::io::savePCDFileBinary(output_map_name, filtered_map);

  LOG(INFO) << GREEN << " ---> Final map saved to: " << output_map_name << RESET;
  LOG(INFO) << GREEN << " ---> Final map size: " << filtered_map.size() << RESET;
}


void SuperLIO::saveMap(){
  if(!g_save_map) return;
  if (!point_map_ || (point_map_->empty() && pcd_index_ < 0)) {
    LOG(INFO) << "No map data collected; skipping map save.";
    return;
  }
  std::error_code directory_error;
  std::filesystem::create_directories(
      g_pcd_save_interval > 0 ? g_save_map_dir + "/PCD" : g_save_map_dir,
      directory_error);
  if (directory_error) {
    LOG(ERROR) << "Cannot create map directory: " << directory_error.message();
    return;
  }
  if(g_pcd_save_interval > 0){
    LOG(INFO) << YELLOW << " ---> Saving last cace ... " << RESET;
    if (point_map_->size() > 0) {
      pcd_index_++;
      std::string map_name(std::string(g_save_map_dir + "/PCD/scans_") + std::to_string(pcd_index_) +
                                 std::string(".pcd"));
      LOG(INFO) << GREEN << " ---> current scan saved to /PCD/scans_" << pcd_index_ << "  size:  " << point_map_->size() << RESET;
      pcl::io::savePCDFileBinary(map_name, *point_map_);
      point_map_->clear();
    }
    LOG(INFO) << GREEN << " ---> Save last cace success. " << RESET;
    LOG(INFO) << YELLOW << " ---> Process cace map ... " << RESET;
    ProcessCaceMap();
    return;
  }

  LOG(INFO) << YELLOW << " ---> Saving map..... " << RESET;
  if(!point_map_->empty()){
    std::string map_name = g_save_map_dir + "/" + g_map_name;
    LOG(INFO) << YELLOW << " ---> Save map to: " << map_name << RESET;
    pcl::VoxelGrid<PointType> voxel_fliter;
    PointCloudType latst_map;
    voxel_fliter.setInputCloud(point_map_);
    voxel_fliter.setLeafSize(g_map_ds_size, g_map_ds_size, g_map_ds_size);
    voxel_fliter.filter(latst_map);
    if(latst_map.size() > 0){
      latst_map.width = latst_map.size();
      latst_map.height = 1;
      latst_map.is_dense = false;
    }
    pcl::io::savePCDFileBinary(map_name, latst_map);
    LOG(INFO) << GREEN << " ---> Save map success. File: " << map_name << RESET;
    LOG(INFO) << GREEN << " ---> Map size: " << latst_map.size() << RESET;
  }
}


void SuperLIO::Propagation_Undistort(){
  propagate_states_.clear();
  propagate_states_.emplace_back(kf_->GetDynamicState());
  kf_->SetObsTime(measures_.lidar.end_time);
  for (auto &imu : measures_.imu) {
    if (kf_->Predict(imu) && kf_->GetTime() > propagate_states_.back().time)
      propagate_states_.emplace_back(kf_->GetDynamicState());
  }

  const M3 TLI_R = g_lidar_imu.R_;
  const V3 TLI_t = g_lidar_imu.t_;
  const double start_time = measures_.lidar.start_time;
  auto& raw_pc = measures_.lidar.pc;

  std::size_t ptsize = raw_pc->points.size();
  scan_undistort_full_->resize(ptsize); 

  tbb::parallel_for(
  tbb::blocked_range<size_t>(0, ptsize),
  [&](const tbb::blocked_range<size_t>& r) {
    for (size_t idx = r.begin(); idx < r.end(); ++idx) {  
      auto& pt_full = scan_undistort_full_->points[idx];
      const auto& pt = raw_pc->points[idx];
      pt_full.intensity = pt.intensity;
      const double query_time = start_time + pt.offset_time;
      const V3 raw(pt.x, pt.y, pt.z);
      const V3 eigen_point = deskewPoint(TLI_R * raw + TLI_t, query_time, propagate_states_);
      pt_full.x = eigen_point[0];
      pt_full.y = eigen_point[1];
      pt_full.z = eigen_point[2];
    }
  });
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
    G += 1000 * J * J.transpose();
  }
  degeneracy_ = geometry_analyzer_.analyze(G, geometry_options_, true);
  sampling_geometry_analyzed_ = true;
}


struct ThreadACC{
  M6d HTVH = M6d::Zero();
  V6d HTVr = V6d::Zero();
  ThreadACC(): HTVH(M6d::Zero()), HTVr(V6d::Zero()) {}
};


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
            // Plane reliability is handled by the original correspondence/residual gates.
            // Keep Rp=0.001 in every geometry state; only selected bump rows get their own Rb.
            local_acc.HTVH += J * 1000 * J.transpose();
            local_acc.HTVr -= J * 1000 * error;
          }
        }
    });

    M6d sum_HTVH = M6d::Zero();
    V6d sum_HTVr = V6d::Zero();
    for(const auto& local_acc : tls_acc){
      sum_HTVH += local_acc.HTVH;
      sum_HTVr += local_acc.HTVr;
    }
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
                                degeneracy_, geometry_options_, c, geometry_stats_)) candidates.push_back(c);
      }
      geometry::selectWeakDirectionConstraints(candidates, degeneracy_, geometry_options_);
      // Strategy A: only final accepted indices replace their own plane rows.
      geometry::applyExclusiveBumps(plane_J, plane_error, effect_mask_, candidates, sum_HTVH, sum_HTVr);
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
