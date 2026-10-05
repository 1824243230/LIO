#include "lio/plane_fit.h"

#include "lio/super_lio_reloc.h"

#include <sys/resource.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>

#include <pcl/registration/icp.h>
#include <pcl/registration/ndt.h>
#include <pcl/kdtree/kdtree_flann.h>


using namespace BASIC;

namespace LI2Sup{

inline bool compute_error(
  const std::array<double, 4>& abcd, const V3& point, 
  const float length, scalar& error)
{
  error = abcd[0] * point[0] + abcd[1] * point[1] + abcd[2] * point[2] + abcd[3];
  return length > 81 * error * error;
}


void SuperLIOReLoc::init(){
  ivox_.reset(new OctVoxMapType(OctVoxMapType::Options{g_ivox_resolution, g_ivox_capacity}));
  kf_.reset(new ESKF());
  InitGeometry();
  data_wrapper_->setESKF(kf_);
  
  scan_undistort_full_.reset(new PointCloudType());
  ds_undistort_.reset(new PointCloudType());
  world_pc_.reset(new PointCloudType());
  ds_world_.reset(new PointCloudType());

  point_map_.reset(new PointCloudType());

  init_obs_data_.reset(new PointCloudType());
  
  points_world_v3_.reserve(21000);
  abcd_vec_.resize(20000);
  effect_knn_idxs_.resize(20000);
  voxel_grid_fliter_.setLeafSize(g_voxel_fliter_size);

  LOG(INFO) << GREEN << " ---> [SuperLIO]: initialized." << RESET;

  auto start_time = std::chrono::high_resolution_clock::now();
  if (!SuperLIOReLoc::map_init())
    throw std::runtime_error("Relocation requires a readable, nonempty map");
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  LOG(INFO) << GREEN << " ---> [SuperLIO]: Map init success. Time: " << duration.count() << " ms." << RESET;

  state_fn_ = &SuperLIOReLoc::stateWaitKFInit;
}


bool SuperLIOReLoc::map_init(){
  if (map_loaded_) return true;

  std::string map_name = g_save_map_dir + "/" + g_map_name;
  try {
    if (pcl::io::loadPCDFile<PointType>(map_name, *point_map_) < 0) return false;
  } catch (const std::exception& error) {
    LOG(ERROR) << "Cannot load relocation map " << map_name << ": " << error.what();
    return false;
  }
  // 文件头的 is_dense 可能不准确，逐点检查后再进入地图/配准器。
  auto& points = point_map_->points;
  points.erase(std::remove_if(points.begin(),points.end(),[](const PointType& p) {
    return !p.getVector3fMap().allFinite();
  }),points.end());
  point_map_->width=points.size();point_map_->height=1;point_map_->is_dense=true;
  if (point_map_->empty()) return false;


  VV3 point_map_v3;
  point_map_v3.reserve(point_map_->size());
  for(const auto& point: *point_map_){
    V3 pt(point.x, point.y, point.z);
    point_map_v3.push_back(pt);
  }

  ivox_->insert(point_map_v3);
  if (bump_map_) {
    std::vector<geometry::V3> points;
    for (const auto& p : point_map_v3) points.push_back(p.cast<double>());
    bump_map_->insert(points, geometry::V3::Zero());
  }

  LOG(INFO) << GREEN << " ---> Load map success. File: " << map_name << RESET;
  LOG(INFO) << GREEN << " ---> Map size: " << point_map_->size() << RESET;
  ivox_->printInfo();

  map_loaded_ = true;

  data_wrapper_->set_global_map(point_map_);
  data_wrapper_->set_initial_data(re_init_pose_, flg_get_init_guess_);
  return true;
}



bool SuperLIOReLoc::kf_init(){
  if (measures_.imu.empty()) return false;
  const int need_init_frames = 10;
  /// get init guess from ROS topic.
  if(flg_get_init_guess_){
    imu_initialization_ = ImuInitialization{};
    init_frame_count_ = 0;
    init_obs_data_->clear();
    flg_get_init_guess_ = false;
    return false;
  }

  CloudPtr point_cloud_pcl = CloudPtr(new PointCloudType());
  for(std::size_t i = 0; i < measures_.lidar.pc->size(); i++){
    auto p = measures_.lidar.pc->points[i];
    PointType point;
    point.x = p.x;
    point.y = p.y;
    point.z = p.z;
    point.intensity = p.intensity;
    point_cloud_pcl->points.push_back(point);
  }

  if(init_frame_count_ < need_init_frames){
    *init_obs_data_ += *point_cloud_pcl;
  }
  init_frame_count_++;

  for (const auto& imu : measures_.imu) imu_initialization_.add(imu);
  const auto& mean_gyro = imu_initialization_.mean_gyro;
  const auto& mean_acce = imu_initialization_.mean_acc;

  if(imu_initialization_.count < 20){
    return false;
  }

  if(init_frame_count_ < need_init_frames){
    return false;
  }

  LOG(INFO) << YELLOW << " ---> INIT start... obs_data size: " << init_obs_data_->size() << " target size: " << point_map_->size() << RESET;

  // 零加速度均值无法确定重力方向，继续等待数据，避免归一化产生 NaN。
  if (!mean_acce.allFinite() || !mean_gyro.allFinite() || mean_acce.norm() < 1e-6 ||
      !std::isfinite(g_gravity_norm) || g_gravity_norm <= 0) return false;
  V3 gravity = - mean_acce * g_gravity_norm / mean_acce.norm();
  V3 ref_gravity(0, 0, - g_gravity_norm);
  M3 init_rot = Quat::FromTwoVectors(gravity, ref_gravity).toRotationMatrix();
  V3 n = init_rot.col(0);
  double yaw = atan2(n(1), n(0));

  M3 R_yaw_inv = Eigen::AngleAxis<scalar>(-yaw, V3::UnitZ()).toRotationMatrix(); 
  M3 rot = R_yaw_inv * init_rot;

  M3 init_guess_R_ = re_init_pose_.R_ * rot;
  V3 init_guess_t_ = re_init_pose_.t_;
  M4 init_guess_T = M4::Identity();
  init_guess_T.block<3, 3>(0, 0) = init_guess_R_;
  init_guess_T.block<3, 1>(0, 3) = init_guess_t_;

  pcl::PointCloud<pcl::PointXYZI>::Ptr tmp_src(new pcl::PointCloud<pcl::PointXYZI>());
  pcl::transformPointCloud(*init_obs_data_, *tmp_src, g_lidar_imu.matrix().cast<float>());

  pcl::NormalDistributionsTransform<pcl::PointXYZI, pcl::PointXYZI> ndt;
  ndt.setTransformationEpsilon(1e-4);
  ndt.setEuclideanFitnessEpsilon(1e-4);
  ndt.setMaximumIterations(25);
  ndt.setResolution(1.0);
  ndt.setInputTarget(point_map_);

  pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
  icp.setMaxCorrespondenceDistance(4.0);
  icp.setMaximumIterations(40);
  icp.setTransformationEpsilon(1e-4);
  icp.setEuclideanFitnessEpsilon(1e-4);
  icp.setRANSACIterations(0);
  icp.setInputTarget(point_map_);

  ndt.setInputSource(tmp_src);
  icp.setInputSource(tmp_src);

  pcl::PointCloud<pcl::PointXYZI>::Ptr unused_result(new pcl::PointCloud<pcl::PointXYZI>());
  ndt.align(*unused_result, init_guess_T.matrix().cast<float>());
  icp.align(*unused_result, ndt.getFinalTransformation());

  if (icp.hasConverged() == false || icp.getFitnessScore() > 1.5)
  // if (icp.hasConverged() == false)
  {
    /// reset init state.
    imu_initialization_ = ImuInitialization{};
    init_frame_count_ = 0;
    init_obs_data_->clear();
    LOG(INFO) << RED << " ---> Global ICP Converged Fail! FitnessScore: " << icp.getFitnessScore() << RESET;
    return false;
  } else{
    init_guess_T = icp.getFinalTransformation().cast<scalar>();
    LOG(INFO) << GREEN << " ---> Global ICP Converged Succeed! FitnessScore: " << icp.getFitnessScore() << RESET;
  }

  LOG(INFO) << GREEN << "\n" << init_guess_T << RESET;

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
  // 注册得到的位姿与当前 IMU 一起提交，不能再用旧版的 -1 时间哨兵。
  kf_->SetPoseAtImu(SE3(SO3(init_guess_T.block<3,3>(0,0)),
                       init_guess_T.block<3,1>(0,3)), measures_.imu.back());
  sys_init_pose_ = kf_->GetSE3();

  {
    point_map_->clear();
    point_map_.reset(new PointCloudType());
    init_obs_data_->clear();
    init_obs_data_ = nullptr;
    data_wrapper_->set_initial_data(re_init_pose_, flg_get_init_guess_, true);
  }

  return true;
}


void SuperLIOReLoc::UpdateMap() {
  if(!g_update_map) return;

  if(map_update_delay_ > 0){
    map_update_delay_--;
    std::cout << "\rUpdate map Delay: "
            << 100 - map_update_delay_
            << " %" << std::flush;
    return;
  }

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


void SuperLIOReLoc::Output() {
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



} // namespace END.