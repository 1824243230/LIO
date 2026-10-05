#include "lio/super_lio.h"
#include <stdexcept>
#include <iostream>
using namespace LI2Sup;
using namespace BASIC;
void check(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
class Harness : public SuperLIO {
 public:
  void prepare() {
    kf_=std::make_shared<ESKF>();
    scan_undistort_full_.reset(new PointCloudType());
    world_pc_.reset(new PointCloudType());
    point_map_.reset(new PointCloudType());
    PointType p; p.x=1; p.y=2; p.z=3; p.intensity=1;
    scan_undistort_full_->push_back(p);
  }
  void cache() { caceData(); }
  bool retained() const { return point_map_->size()==1 && pcd_index_==-1; }
  bool zero_acceleration_init() {
    measures_.imu.resize(50);
    for(size_t i=0;i<measures_.imu.size();++i) {
      measures_.imu[i].secs=1+i*0.01;
      measures_.imu[i].acc.setZero(); measures_.imu[i].gyr.setZero();
    }
    return kf_init();
  }
};
int main() {
  // 机体朝世界 +Y，但沿世界 +X 运动：ROS 子坐标系中的速度应为 -Y。
  // 同时验证两条实际发布路径共用的消息构造，不需要启动 ROS master。
  NavState nav;
  nav.timestamp=42.125;
  nav.R=SO3(M3(Eigen::AngleAxis<scalar>(M_PI/2, V3::UnitZ())));
  nav.p=V3(4,5,6); nav.v=V3(2,0,0);
  const auto low_rate=makeImuOdometry(nav);
  check(low_rate.header.frame_id=="world" && low_rate.child_frame_id=="body" &&
        low_rate.header.stamp.toSec()==nav.timestamp, "odometry declares both frames and state timestamp");
  check(low_rate.pose.pose.position.x==4 && low_rate.pose.pose.position.y==5 &&
        low_rate.pose.pose.position.z==6, "odometry pose remains in world frame");
  check(std::abs(low_rate.twist.twist.linear.x)<1e-6 &&
        std::abs(low_rate.twist.twist.linear.y+2)<1e-6,
        "world velocity rotates to the declared child frame");
  DynamicState dynamic(nav.timestamp,nav.R.R_,nav.p,nav.v,V3(.1,.2,.3),V3::Zero());
  const auto high_rate=makeImuOdometry(dynamic);
  check(high_rate.header==low_rate.header && high_rate.child_frame_id==low_rate.child_frame_id &&
        high_rate.pose.pose==low_rate.pose.pose && high_rate.twist.twist.linear==low_rate.twist.twist.linear,
        "IMU and LiDAR rate odometry use the same frame convention");
  check(std::abs(high_rate.twist.twist.angular.x-.1)<1e-6 &&
        std::abs(high_rate.twist.twist.angular.y-.2)<1e-6 &&
        std::abs(high_rate.twist.twist.angular.z-.3)<1e-6,
        "body angular velocity is preserved without another rotation");

  g_blind2=0; g_maxrange2=100;
  auto msg=boost::make_shared<livox_ros_driver::CustomMsg>();
  CloudPtr cloud;
  livox2pcl(msg,cloud);
  check(cloud && cloud->empty(), "empty Livox packet is safe with null output");
  msg->point_num=3; msg->points.resize(2);
  livox2pcl(msg,cloud); check(cloud->empty(), "declared size mismatch rejected");
  msg->points.resize(3);
  for(size_t i=0;i<3;++i) { msg->points[i].x=i+1; msg->points[i].tag=0; }
  livox2pcl(msg,cloud);
  check(cloud->size()==2 && (*cloud)[0].x==2 && (*cloud)[1].x==3, "stable point order");
  msg->points[2].x=msg->points[1].x;
  livox2pcl(msg,cloud); check(cloud->size()==1, "adjacent duplicates removed");
  msg->points[1].tag=0x20;
  livox2pcl(msg,cloud); check(cloud->empty(), "tag rejection and raw-neighbor comparison");
  g_save_map=true; g_if_filter=false; g_pcd_save_interval=0;
  Harness h; h.prepare(); h.cache();
  check(h.retained(), "zero save interval retains points for final save without fragments");
  g_gravity_norm=9.81;
  check(!h.zero_acceleration_init(), "zero acceleration cannot initialize gravity");
  std::cout << "input and cache review regression checks passed\n";
}
