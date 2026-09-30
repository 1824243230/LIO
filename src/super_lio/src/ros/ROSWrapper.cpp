
#include "ros/ROSWrapper.h"
#include "lio/scan_validation.h"
#include "lio/measurement_sync.h"
#include "super_lio/CloudPose.h"
#include "super_lio/CloudPose2.h"

#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <stdexcept>

using namespace BASIC;

namespace LI2Sup{

void LoadParamFromRos(ros::NodeHandle& nh){
  // 关键几何参数不能用缺省零值继续运行：错误配置必须在构造地图/滤波器前退出。
  const auto required = [&](const char* key, auto& value) {
    if (!nh.getParam(key, value)) throw std::invalid_argument(std::string("Missing/invalid ROS parameter: ")+key);
  };
  const auto require_finite_array = [](const auto& values, size_t size, const char* key) {
    if (values.size()!=size || !std::all_of(values.begin(),values.end(),
        [](auto v) { return std::isfinite(v); }))
      throw std::invalid_argument(std::string("Invalid extrinsic length/value: ")+key);
  };
  nh.getParam("/lio/map/save_map", g_save_map);
  LOG(INFO) << GREEN << " ---> [Param] map/save_map: " << (g_save_map ? "true" : "false") << RESET;
  nh.getParam("/lio/map/if_filter", g_if_filter);
  nh.getParam("/lio/map/save_map_dir", g_save_map_dir);
  g_save_map_dir = g_root_dir + g_save_map_dir;
  nh.getParam("/lio/map/map_name", g_map_name);
  nh.getParam("/lio/map/ds_size", g_map_ds_size);
  nh.getParam("/lio/map/save_interval", g_pcd_save_interval);

  nh.getParam("/lio/eva/timer", g_time_eva);
  
  // ROS Topic input
  nh.getParam("/lio/ros/lidar_topic",  g_lidar_topic);
  nh.getParam("/lio/ros/imu_topic",    g_imu_topic);

  // sensor cfg
  nh.getParam("/lio/sensor/lidar_type", g_lidar_type);
  double temp_range_dis;
  required("/lio/sensor/blind", temp_range_dis);
  if (!std::isfinite(temp_range_dis) || temp_range_dis < 0) throw std::invalid_argument("Invalid blind range");
  g_blind2 = temp_range_dis * temp_range_dis;
  required("/lio/sensor/maxrange", temp_range_dis);
  if (!std::isfinite(temp_range_dis) || temp_range_dis <= std::sqrt(g_blind2))
    throw std::invalid_argument("maxrange must exceed blind range");
  g_maxrange2 = temp_range_dis * temp_range_dis;
  nh.getParam("/lio/sensor/filter_rate", g_filter_rate);
  g_filter_rate = std::max(1, g_filter_rate);
  nh.getParam("/lio/sensor/enable_downsample", g_enable_downsample);
  nh.getParam("/lio/sensor/voxel_fliter_size", g_voxel_fliter_size);

  nh.getParam("/lio/sensor/gravity_norm", g_gravity_norm);
  nh.getParam("/lio/sensor/imu_type", g_imu_type);
  nh.getParam("/lio/sensor/imu_na",   g_imu_na);
  nh.getParam("/lio/sensor/imu_ng",   g_imu_ng);
  nh.getParam("/lio/sensor/imu_nba",  g_imu_nba);
  nh.getParam("/lio/sensor/imu_nbg",  g_imu_nbg);

  // extrinsic
  std::vector<scalar> extrinsic_lidar_imu, extrinsic_odom_robo;
  required("/lio/extrinsic/lidar_imu", extrinsic_lidar_imu);
  require_finite_array(extrinsic_lidar_imu, 12, "lidar_imu"); // 平移 3 项 + 原有存储约定的旋转 9 项
  V3 __t = V3(extrinsic_lidar_imu[0], 
              extrinsic_lidar_imu[1], 
              extrinsic_lidar_imu[2]);
  M3 __R = M3(extrinsic_lidar_imu.data() + 3);
  g_lidar_imu = SE3(__R, __t);  // lidar in imu frame·
  required("/lio/extrinsic/odom_robo", extrinsic_odom_robo);
  require_finite_array(extrinsic_odom_robo, 6, "odom_robo"); // x,y,z,roll,pitch,yaw，角度单位为度
  __t = V3(extrinsic_odom_robo[0], 
           extrinsic_odom_robo[1], 
           extrinsic_odom_robo[2]);
  auto temp_R = Eigen::AngleAxisd(extrinsic_odom_robo[5]/180 * M_PI, Eigen::Vector3d::UnitZ()) *
                Eigen::AngleAxisd(extrinsic_odom_robo[4]/180 * M_PI, Eigen::Vector3d::UnitY()) *
                Eigen::AngleAxisd(extrinsic_odom_robo[3]/180 * M_PI, Eigen::Vector3d::UnitX());
  g_odom_robo.R_ = temp_R.cast<scalar>();

  /// ATTENTION:
  /// The transpose is intentionally applied here and represents the inverse of R.
  /// Misinterpreting this will lead to incorrect transformations.
  M3 _R = g_odom_robo.R_.transpose();
  g_odom_robo.R_  = _R;
  g_odom_robo = SE3(_R, __t);  // lidar in robot frame

  auto temp_R_yaw = Eigen::AngleAxisd(extrinsic_odom_robo[5]/180 * M_PI, Eigen::Vector3d::UnitZ()).toRotationMatrix();
  g_lidar_robo_yaw = temp_R_yaw.cast<scalar>();

  // hash_map
  int hash_capacity;
  required("/lio/hash_map/hash_capacity", hash_capacity);
  if (hash_capacity <= 0) throw std::invalid_argument("hash_capacity must be positive");
  g_ivox_capacity = hash_capacity;
  nh.getParam("/lio/hash_map/vox_resolution", g_ivox_resolution);
  
  // kf
  nh.getParam("/lio/kf/kf_type", g_kf_type);
  nh.getParam("/lio/kf/kf_max_iterations", g_kf_max_iterations);
  nh.getParam("/lio/kf/kf_align_gravity", g_kf_align_gravity);
  nh.getParam("/lio/kf/kf_quit_eps", g_kf_quit_eps);

  // submaps
  nh.getParam("/lio/submap/submap_resolution", g_submap_resolution);
  nh.getParam("/lio/submap/submap_capacity", g_submap_capacity);

  // visual
  nh.getParam("/lio/output/robot",  g_2_robot);
  nh.getParam("/lio/output/planner", g_planner_enable);
  nh.getParam("/lio/output/plan_env_world",  g_2_plan_env_world);
  nh.getParam("/lio/output/plan_env_body",  g_2_plan_env_body);
  nh.getParam("/lio/output/ml_map",         g_2_ml_map);
  nh.getParam("/lio/output/map",    g_visual_map);
  nh.getParam("/lio/output/dense",  g_visual_dense);
  nh.param("/lio/output/pub_step", g_pub_step, 1);
  g_pub_step = std::max(1, g_pub_step); // 发布计数使用取模，禁止除数为零。

  g_update_map = false;
  nh.getParam("/lio/relocation/update_map", g_update_map);
  std::vector<float> init_pose;
  nh.getParam("/lio/relocation/init_pose", init_pose);
  if(init_pose.size() == 6){
    g_init_px    = init_pose[0];
    g_init_py    = init_pose[1];
    g_init_pz    = init_pose[2];
    g_init_roll  = init_pose[3];
    g_init_pitch = init_pose[4];
    g_init_yaw   = init_pose[5];
  }else{
    g_init_px = 0.0f;
    g_init_py = 0.0f;
    g_init_pz = 0.0f;
    g_init_roll = 0.0f;
    g_init_pitch = 0.0f;
    g_init_yaw = 0.0f;
  }
  
}

std::tuple<float, float, float> getColorFromVelocity(float velocity, float max_velocity = 2.0f) {
  float ratio = std::clamp(velocity / max_velocity, 0.0f, 1.0f);
  float r, g, b;
  if (ratio < 0.5f) {
    float t = ratio / 0.5f;
    r = 0.0f;
    g = t;
    b = 1.0f;
  } else {
    float t = (ratio - 0.5f) / 0.5f;
    r = t;
    g = 1.0f - t;
    b = 1.0f - t;
  }
  return {r, g, b};
}


void livox2pcl(const livox_ros_driver::CustomMsg::ConstPtr& msg, CloudPtr& point_cloud){
  if (!point_cloud) point_cloud.reset(new PointCloudType());
  point_cloud->clear();
  // 本函数为备用无时间点云转换接口；实际 Livox 回调使用带 offset_time 的独立路径。
  // 必须先验证声明长度，避免空包 size-1 下溢及错误长度造成越界。
  if (!msg || msg->point_num != msg->points.size() || msg->points.size() < 2) return;
  point_cloud->reserve(msg->points.size());
  // 直接读不可变原消息中的相邻点，避免并行读取尚未写完的 cloud_full[i-1]。
  // 保留原接口从第 1 点开始、按相邻坐标去重的约定，不使用 vector<bool> 并发写入。
  for (size_t i=1; i<msg->points.size(); ++i) {
    const auto& p=msg->points[i];
    const auto& previous=msg->points[i-1];
    if ((p.tag & 0x30)!=0x10 && (p.tag & 0x30)!=0x00) continue;
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) continue;
    if (std::abs(p.x-previous.x)<=1e-7 && std::abs(p.y-previous.y)<=1e-7 &&
        std::abs(p.z-previous.z)<=1e-7) continue;
    const double distance=double(p.x)*p.x+double(p.y)*p.y+double(p.z)*p.z;
    if (distance<=g_blind2 || distance>=g_maxrange2) continue;
    PointType point;
    point.x=p.x; point.y=p.y; point.z=p.z; point.intensity=p.reflectivity;
    point_cloud->push_back(point);
  }
}


std::string lidarTypeToString(int type) {
  if (type <= 0 || type >= static_cast<int>(LID_TYPE_NAMES.size())) return "UNKNOWN";
  return LID_TYPE_NAMES[type];
}

ROSWrapper::ROSWrapper(){
  ros::SubscribeOptions ops;
  ops.transport_hints = ros::TransportHints().tcpNoDelay();
  
  if(g_lidar_type == LID_TYPE::LIVOX){
    ops.init<livox_ros_driver::CustomMsg>(
      g_lidar_topic, 1000, 
      boost::bind(&ROSWrapper::livoxHandler, this, _1));
  }else{
    ops.init<sensor_msgs::PointCloud2>(
      g_lidar_topic, 1000,
      boost::bind(&ROSWrapper::stdMsgHandler, this, _1));
  }

  LOG(INFO) << GREEN << " ---> Using Lidar type: " << lidarTypeToString(g_lidar_type) << RESET;

  nh_.setCallbackQueue(&self_queue_);

  subLidar_ = nh_.subscribe(ops);
  subIMU_   = nh_.subscribe<sensor_msgs::Imu>(g_imu_topic, 10000,    // 100Hz x 10s
               &ROSWrapper::imuHandler, this, ros::TransportHints().tcpNoDelay());

  /// output
  pub_odom_      = nh_.advertise<nav_msgs::Odometry>("/lio/odom", 100);  /// imu frame -> lidar freq
  pub_path_      = nh_.advertise<nav_msgs::Path>("/lio/path", 1);
  pub_path_robot_ = nh_.advertise<sensor_msgs::PointCloud2>("/lio/path_robot", 1);
  
  msg_path_point_.header.frame_id = "world";
  msg2uav_.header.frame_id = "world";
  path_.header.frame_id = "world";
}


void ROSWrapper::livoxHandler(const livox_ros_driver::CustomMsg::ConstPtr& msg){
  LidarData scan;
  if (!msg || !convertLivoxScan(*msg,g_filter_rate,g_blind2,g_maxrange2,scan)) return;
  lidar_buffer_.push_back(std::move(scan));
}


inline bool validPoint(double x, double y, double z)
{
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
    return false;

  double d2 = x * x + y * y + z * z;
  return (d2 > g_blind2 && d2 < g_maxrange2);
}

void ROSWrapper::stdMsgHandler(const sensor_msgs::PointCloud2::ConstPtr& msg){
  if(msg->data.size() < 10) return;
  
  LidarData lidar_data;
  lidar_data.pc.reset(new pcl::PointCloud<LI2Sup::PointXTZIT>());

  double offset_time = 0.0;
  double dis = 0.0;

  switch (g_lidar_type) {

  case LID_TYPE::HESAI16:
  {
    pcl::PointCloud<hesai_ros::Point> pl_orig;
    pcl::fromROSMsg(*msg, pl_orig);
    if (pl_orig.empty()) return;
    lidar_data.pc->reserve(pl_orig.size() / g_filter_rate + 1);
    const double time_begin = pl_orig.points[0].timestamp;
    lidar_data.start_time = time_begin;
    for(std::size_t i = 0; i < pl_orig.size(); i += g_filter_rate)
    {
      auto& pt = pl_orig.points[i];
      if (!validPoint(pt.x, pt.y, pt.z)) continue;
      offset_time = pt.timestamp - time_begin;
      lidar_data.pc->emplace_back(
          pt.x, pt.y, pt.z, pt.intensity, offset_time);
    }
    lidar_data.end_time = time_begin + offset_time;
    break;
  }
  case LID_TYPE::VEL_NCLT:
  {
    pcl::PointCloud<NCLT::Point> pl_orig;
    pcl::fromROSMsg(*msg, pl_orig);
    lidar_data.pc->reserve(pl_orig.size() / g_filter_rate + 1);
    lidar_data.start_time = msg->header.stamp.toSec();
    
    for(std::size_t i = 0; i < pl_orig.size(); i += g_filter_rate){
      auto& pt = pl_orig.points[i];
      if (!validPoint(pt.x, pt.y, pt.z)) continue;
      offset_time = pt.time * 1e-6;
      lidar_data.pc->emplace_back(
          pt.x, pt.y, pt.z, 1.0, offset_time);
    }
    lidar_data.end_time = lidar_data.start_time + offset_time;
    break;
  }
  case LID_TYPE::VELO16:
  case LID_TYPE::VELO32:
  {
    pcl::PointCloud<velodyne_ros::Point> pl_orig;
    pcl::fromROSMsg(*msg, pl_orig);
    lidar_data.pc->reserve(pl_orig.size() / g_filter_rate + 1);
    lidar_data.start_time = msg->header.stamp.toSec();

    for(std::size_t i = 0; i < pl_orig.size(); i += g_filter_rate){
      auto& pt = pl_orig.points[i];
      if (!validPoint(pt.x, pt.y, pt.z)) continue;
      lidar_data.pc->emplace_back(
          pt.x, pt.y, pt.z, pt.intensity, pt.time);
    }
    break;
  }
  case OUSTER:
  {
    pcl::PointCloud<ouster_ros::Point> pl_orig;
    pcl::fromROSMsg(*msg, pl_orig);
    lidar_data.pc->reserve(pl_orig.size() / g_filter_rate + 1);
    lidar_data.start_time = msg->header.stamp.toSec();

    for(std::size_t i = 0; i < pl_orig.size(); i += g_filter_rate){
      auto& pt = pl_orig.points[i];
      if (!validPoint(pt.x, pt.y, pt.z)) continue;
      offset_time = pt.t * 1e-9;
      lidar_data.pc->emplace_back(
          pt.x, pt.y, pt.z, pt.intensity, offset_time);
    }
    lidar_data.end_time = lidar_data.start_time + offset_time;
    break;
  }
  default:
    return;
  }
  
  if (!finalizeLidarScan(lidar_data)) return;
  lidar_buffer_.push_back(lidar_data);
}



void ROSWrapper::imuHandler(const sensor_msgs::Imu::ConstPtr& msg){
  IMUData data;
  data.secs = msg->header.stamp.toSec();
  data.acc  = V3(msg->linear_acceleration.x,
                 msg->linear_acceleration.y,
                 msg->linear_acceleration.z);
  data.gyr  = V3(msg->angular_velocity.x,
                 msg->angular_velocity.y,
                 msg->angular_velocity.z);
  if (!std::isfinite(data.secs) || !data.acc.allFinite() || !data.gyr.allFinite()) return;

  if (!appendImu(data,imu_buffer_,last_timestamp_imu_)) {
    if (data.secs < last_timestamp_imu_)
      ROS_WARN_THROTTLE(1.0, "Dropping out-of-order IMU; restart the node when replaying an earlier time range");
    return;
  }
  if (!eskf_) return;

  static ros::Publisher pub_imu_odom  = nh_.advertise<nav_msgs::Odometry>("/lio/imu/odom", 10);    /// imu frame -> imu freq
  static ros::Publisher pub_robo_odom = nh_.advertise<nav_msgs::Odometry>("/lio/robo/odom", 10);   /// robot frame -> imu freq
  
  DynamicState imu_state, robo_state;
  if(eskf_->Predict(data, imu_state, robo_state)){
    nav_msgs::Odometry odom_imu, odom_robo;

    {
      odom_imu.pose.pose.position.x = imu_state.p(0);
      odom_imu.pose.pose.position.y = imu_state.p(1);
      odom_imu.pose.pose.position.z = imu_state.p(2);

      Quat q(imu_state.R);
      q.normalize();

      odom_imu.pose.pose.orientation.x = q.x();
      odom_imu.pose.pose.orientation.y = q.y();
      odom_imu.pose.pose.orientation.z = q.z();
      odom_imu.pose.pose.orientation.w = q.w();

      odom_imu.twist.twist.linear.x = imu_state.v(0);
      odom_imu.twist.twist.linear.y = imu_state.v(1);
      odom_imu.twist.twist.linear.z = imu_state.v(2);

      odom_imu.twist.twist.angular.x = imu_state.w(0);
      odom_imu.twist.twist.angular.y = imu_state.w(1);
      odom_imu.twist.twist.angular.z = imu_state.w(2);
    }

    {
      odom_robo.pose.pose.position.x = robo_state.p(0);
      odom_robo.pose.pose.position.y = robo_state.p(1);
      odom_robo.pose.pose.position.z = robo_state.p(2);

      Quat q(robo_state.R);
      q.normalize();

      odom_robo.pose.pose.orientation.x = q.x();
      odom_robo.pose.pose.orientation.y = q.y();
      odom_robo.pose.pose.orientation.z = q.z();
      odom_robo.pose.pose.orientation.w = q.w();
    }

    odom_imu.header.stamp = msg->header.stamp;
    odom_robo.header.stamp = msg->header.stamp;
    odom_imu.header.frame_id = "world";
    odom_robo.header.frame_id = "world";
    pub_imu_odom.publish(odom_imu);
    pub_robo_odom.publish(odom_robo);
  }
}


bool ROSWrapper::sync_measure(MeasureGroup& meas){
  return synchronizeMeasurements(lidar_buffer_,imu_buffer_,last_timestamp_lidar_,meas);
}


void ROSWrapper::pub_odom(const NavState& state){
  nav_msgs::Odometry odom;
  odom.header.frame_id = "world";
  odom.header.stamp = ros::Time().fromSec(state.timestamp);
  odom.pose.pose.position.x = state.p[0];
  odom.pose.pose.position.y = state.p[1];
  odom.pose.pose.position.z = state.p[2];

  V4 temp_q = state.R.coeffs();
  odom.pose.pose.orientation.x = temp_q[0];
  odom.pose.pose.orientation.y = temp_q[1];
  odom.pose.pose.orientation.z = temp_q[2];
  odom.pose.pose.orientation.w = temp_q[3];

  odom.twist.twist.linear.x = state.v[0];
  odom.twist.twist.linear.y = state.v[1];
  odom.twist.twist.linear.z = state.v[2];

  pub_odom_.publish(odom);    // imu frame -> lidar frequency

  V3 robo_position = state.R.R_ * ( - g_odom_robo.R_ * g_odom_robo.t_) + state.p;

  if(g_2_robot){
    static ros::Publisher pub_msg2uav_ = nh_.advertise<geometry_msgs::PoseStamped>("/mavros/vision_pose/pose", 10);
    M3 robo_rotation = state.R.R_ * g_odom_robo.R_;
    msg2uav_.header.stamp = odom.header.stamp;
    msg2uav_.pose.position.x = robo_position[0];
    msg2uav_.pose.position.y = robo_position[1];
    msg2uav_.pose.position.z = robo_position[2];
    Quat robo_quat(robo_rotation);
    msg2uav_.pose.orientation.w = robo_quat.w();
    msg2uav_.pose.orientation.x = robo_quat.x();
    msg2uav_.pose.orientation.y = robo_quat.y();
    msg2uav_.pose.orientation.z = robo_quat.z();
    pub_msg2uav_.publish(msg2uav_);
  }

  // if(1)
  if((last_path_point_ - robo_position).norm() > 0.1)
  {
    /// nav_msgs::Path
    path_.header.stamp = odom.header.stamp;
    geometry_msgs::PoseStamped point;
    point.pose = odom.pose.pose;
    path_.poses.push_back(point);
    pub_path_.publish(path_);

    /// sensor_msgs::PointCloud2
    // pcl::PointXYZRGB point_robot;
    // point_robot.x = robo_position[0];
    // point_robot.y = robo_position[1];
    // point_robot.z = robo_position[2];
    // auto [r, g, b] = getColorFromVelocity(state.v.norm());
    // point_robot.r = static_cast<uint8_t>(r * 255);
    // point_robot.g = static_cast<uint8_t>(g * 255);
    // point_robot.b = static_cast<uint8_t>(b * 255);
    // point_robot.a = 255;
    // pcl::PointCloud<pcl::PointXYZRGB> path_point_;
    // path_point_.push_back(point_robot);
    // pcl::toROSMsg(path_point_, msg_path_point_);
    // msg_path_point_.header.stamp = odom.header.stamp;
    // msg_path_point_.header.frame_id = "world";
    // pub_path_robot_.publish(msg_path_point_);

    last_path_point_ = robo_position;
  }

  // tf::Transform transform;
  // tf::Quaternion q;
  // transform.setOrigin(tf::Vector3(state.p[0], state.p[1], state.p[2]));
  // q.setX(temp_q[0]);
  // q.setY(temp_q[1]);
  // q.setZ(temp_q[2]);
  // q.setW(temp_q[3]);
  // transform.setRotation(q);
  // br_.sendTransform(tf::StampedTransform(transform, odom.header.stamp, "world", "body"));

  // // Visual: for the best field of view in rviz.
  // q.setX(0);
  // q.setY(0);
  // q.setZ(0);
  // q.setW(1);
  // transform.setRotation(q);
  // br_.sendTransform(tf::StampedTransform(transform, odom.header.stamp, "world", "god"));
}



void ROSWrapper::pub_cloud_world(const CloudPtr& pc,double time){
  static ros::Publisher pub_cloud_world_ = nh_.advertise<sensor_msgs::PointCloud2>
                                            ("/lio/cloud_world", 10);
  sensor_msgs::PointCloud2 cloud;
  pcl::toROSMsg(*pc, cloud);
  cloud.header.frame_id = "world";
  cloud.header.stamp = ros::Time().fromSec(time);
  pub_cloud_world_.publish(cloud);
}


void ROSWrapper::pub_cloud2planner(const CloudPtr& pc, double time){
  static ros::Publisher pub_cloud2robot_ = nh_.advertise<sensor_msgs::PointCloud2>
                                            ("/lio/robo/cloud_world", 10);
  sensor_msgs::PointCloud2 cloud;
  pcl::toROSMsg(*pc, cloud);
  cloud.header.frame_id = "world";
  cloud.header.stamp = ros::Time().fromSec(time);
  pub_cloud2robot_.publish(cloud);
}


void ROSWrapper::pub_cloud_body_pose(const CloudPtr& pc, 
  const NavState& state)
{
  static ros::Publisher pub_output2robot_ = nh_.advertise<super_lio::CloudPose>
                                            ("/lio/body/cloud_pose", 10);
  super_lio::CloudPose cloud_pose;
  pcl::toROSMsg(*pc, cloud_pose.cloud);
  cloud_pose.cloud.header.stamp = ros::Time().fromSec(state.timestamp);  
  cloud_pose.pose.position.x = state.p[0];
  cloud_pose.pose.position.y = state.p[1];
  cloud_pose.pose.position.z = state.p[2];
  V4 temp_q = state.R.coeffs();
  cloud_pose.pose.orientation.x = temp_q[0];
  cloud_pose.pose.orientation.y = temp_q[1];
  cloud_pose.pose.orientation.z = temp_q[2];
  cloud_pose.pose.orientation.w = temp_q[3];

  pub_output2robot_.publish(cloud_pose);
}


void ROSWrapper::pub_cloud_world_pose(const CloudPtr& pc, 
   const NavState& state)
{
  static ros::Publisher pub_output2robot_ = nh_.advertise<super_lio::CloudPose>
                                            ("/lio/world/cloud_pose", 10);
  super_lio::CloudPose cloud_pose;
  pcl::toROSMsg(*pc, cloud_pose.cloud);
  cloud_pose.cloud.header.stamp = ros::Time().fromSec(state.timestamp);  
  cloud_pose.pose.position.x = state.p[0];
  cloud_pose.pose.position.y = state.p[1];
  cloud_pose.pose.position.z = state.p[2];
  V4 temp_q = state.R.coeffs();
  cloud_pose.pose.orientation.x = temp_q[0];
  cloud_pose.pose.orientation.y = temp_q[1];
  cloud_pose.pose.orientation.z = temp_q[2];
  cloud_pose.pose.orientation.w = temp_q[3];
  
  pub_output2robot_.publish(cloud_pose);
}


void ROSWrapper::pub_cloud_body_pose( 
      const BASIC::VV3& pc_body,
      const NavState& state)
{
  static ros::Publisher pub_msg_ = nh_.advertise<super_lio::CloudPose2>
                                            ("/lio/dense/cloud_pose", 10);

  super_lio::CloudPose2 cloud_pose;
  cloud_pose.header.stamp = ros::Time().fromSec(state.timestamp);

  cloud_pose.pose.reserve(12);
  cloud_pose.pose.push_back(state.p[0]);
  cloud_pose.pose.push_back(state.p[1]);
  cloud_pose.pose.push_back(state.p[2]);

  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      cloud_pose.pose.push_back(state.R.R_(r, c));

  cloud_pose.cloud_lidar.reserve(pc_body.size() * 3);
  for (const auto& pt : pc_body) {
    cloud_pose.cloud_lidar.push_back(pt[0]);
    cloud_pose.cloud_lidar.push_back(pt[1]);
    cloud_pose.cloud_lidar.push_back(pt[2]);
  }

  pub_msg_.publish(cloud_pose);
}


void ROSWrapper::pub_processing_time(double time, double current_time, double mean_time, double std_time)
{
  static ros::Publisher pub_process_time_ = nh_.advertise<geometry_msgs::PoseStamped>
                                            ("/lio/processing_time", 10);
  geometry_msgs::PoseStamped msg;
  msg.header.stamp = ros::Time().fromSec(time);
  msg.pose.position.x = current_time;
  msg.pose.position.y = mean_time;
  msg.pose.position.z = std_time;
  pub_process_time_.publish(msg);
}


void ROSWrapper::set_global_map(const BASIC::CloudPtr& global_map){
  pcl::toROSMsg(*global_map, global_map_msg_);
  global_map_msg_.header.frame_id = "world";

  static ros::Publisher global_map_pub =
    nh_.advertise<sensor_msgs::PointCloud2>("/lio/global_map", 1, true);

  static ros::Timer global_map_timer =
    nh_.createTimer(
      ros::Duration(1.0),
      [this](const ros::TimerEvent&) {
        static int count = -1;
        static int publish_interval = 1;
        count++;
        if (count % publish_interval != 0) {
          return;
        }
        count = 0;
        publish_interval++;
        if(publish_interval > 10) publish_interval = 10;
        global_map_msg_.header.stamp = ros::Time::now();
        global_map_pub.publish(global_map_msg_);
      });
}

void ROSWrapper::set_initial_data(BASIC::SE3& init_pose, bool& flg_get_init_guess, bool flg_finish_init)
{
  static ros::Subscriber init_pose_sub =
    nh_.subscribe<geometry_msgs::PoseWithCovarianceStamped>(
      "/initialpose", 1,
      [this, &init_pose, &flg_get_init_guess]
      (const geometry_msgs::PoseWithCovarianceStampedConstPtr& msg)
      {
        V3 init_translation;
        init_translation << 
            msg->pose.pose.position.x,
            msg->pose.pose.position.y,
            0.2;

        double x = msg->pose.pose.orientation.x;
        double y = msg->pose.pose.orientation.y;
        double z = msg->pose.pose.orientation.z;
        double w = msg->pose.pose.orientation.w;

        Quat init_rotation(w, x, y, z);

        init_pose = SE3(SO3(init_rotation.toRotationMatrix()), init_translation);

        flg_get_init_guess = true;
        
        LOG(INFO) << YELLOW
                  << " ---> GET Initial guess: "
                  << init_translation.transpose()
                  << " yaw: "
                  << init_rotation.toRotationMatrix()
                          .eulerAngles(0, 1, 2)
                          .transpose()
                  << RESET;
      });

  if (flg_finish_init) {
    init_pose_sub = ros::Subscriber();
  }
}



} // namespace END.
