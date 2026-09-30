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
