#include "lio/super_lio.h"
#include "lio/measurement_sync.h"
#include "lio/deskew.h"
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace LI2Sup;
using namespace BASIC;
void check(bool ok,const char* message) { if (!ok) throw std::runtime_error(message); }
IMUData sample(double time) { IMUData i; i.secs=time; return i; }
LidarData scan(double start,double end) {
  LidarData l; l.start_time=start;l.end_time=end;l.pc.reset(new pcl::PointCloud<PointXTZIT>());
  for(double t : {start,0.5*(start+end),end}) l.pc->emplace_back(5-(t-1),0,1,1,t-start);
  return l;
}
class Harness : public SuperLIO {
public:
  int outputs=0,updates=0;
  void prepare(bool smooth) {
    smooth_motion_=smooth; uncertainty_motion_=true;
    ESKF::Options o;o.smooth_motion_=smooth;
    kf_=std::make_shared<ESKF>();kf_->SetInitialConditions(o,V3::Zero(),V3::Zero(),1,V3::Zero());
    SysState s;s.timestamp=1;s.v=V3(1,0,0);kf_->SetX(s,sample(1));
    ivox_=std::make_shared<OctVoxMapType>(OctVoxMapType::Options{0.5f,1000});
    scan_undistort_full_.reset(new PointCloudType());ds_undistort_.reset(new PointCloudType());
  }
  void setScan(double start,double end,std::initializer_list<double> imu) {
    measures_.lidar=scan(start,end);measures_.imu.clear();
    for(double t:imu) measures_.imu.push_back(sample(t));
  }
  bool bootstrap() { return map_init(); }
  void run() { stateProcess(); }
  bool propagate() { return Propagation_Undistort(); }
  int frames() const { return frame_num_; }
  double time() const { return kf_->GetTime(); }
  float position() const { return kf_->GetSysState().p.x(); }
  const VV3& world() const { return points_world_v3_; }
  bool valid() const { return observation_valid_; }
  bool emptyCloud() const { return scan_undistort_full_->empty(); }
private:
  void Output() override { ++outputs; }
  void UpdateMap() override { ++updates; }
};
int main() {
  ros::Time::init(); // throttle 日志需要时钟；本测试不启动 ROS 节点或 master。
  g_lidar_imu=SE3();g_time_eva=false;g_save_map=false;
  for(bool smooth : {false,true}) {
    Harness h;h.prepare(smooth);
    for(int i=0;i<4;++i) {
      const double a=1+i*0.1,b=1+(i+1)*0.1;
      h.setScan(a,b,{b});
      check(h.bootstrap()==(i==3),"bootstrap counts completed scans");
      check(h.time()==b && std::abs(h.position()-(b-1))<1e-6,"bootstrap integrates all IMU time");
      for(const auto& p:h.world()) check(std::abs(p.x()-5)<1e-6,"moving bootstrap deskews a fixed world point");
    }
    h.setScan(1.4,1.5,{1.45});h.run();
    check(!h.valid() && h.outputs==0 && h.updates==0 && h.frames()==4 && h.emptyCloud(),
          "incomplete IMU coverage cannot publish or insert scan");
    check(h.time()==1.45,"valid propagation prefix retained");
    h.setScan(1.5,1.6,{1.55,1.6});
    check(h.propagate() && std::abs(h.position()-0.6)<1e-6,"later scan continues from valid prefix without losing time");
    h.setScan(1.5,1.6,{1.6});check(!h.propagate(),"already processed scan rejected");
  }
  // 观测 bookkeeping 不能凭空推进状态，不能造成中间积分区间被跳过。
  ESKF f; f.SetInitialConditions(ESKF::Options(),V3::Zero(),V3::Zero(),1,V3::Zero());
  SysState s;s.timestamp=1;s.v=V3(1,0,0);f.SetX(s,sample(1));
  f.SetLastObsTime(1.08);f.SetObsTime(1.1);
  check(f.Predict(sample(1.1)) && std::abs(f.GetSysState().p.x()-0.1)<1e-6,"nominal IMU endpoint defines integration start");
  for(bool smooth:{false,true}) {
    ESKF::Options o;o.smooth_motion_=smooth;
    ESKF filter;filter.SetInitialConditions(o,V3::Zero(),V3::Zero(),4,V3::Zero());
    s.v.setZero();filter.SetX(s,sample(1));filter.SetObsTime(1.1);filter.init_=true;
    const auto cov=filter.GetCov(); const auto transition=filter.GetTransition();
    IMUData huge=sample(1.1);huge.acc=V3::Constant(std::numeric_limits<float>::max());
    check(!filter.Predict(huge),"overflowing finite input rejected");
    check(filter.GetTime()==1 && (filter.GetCov()-cov).norm()==0 &&
          filter.GetSysState().p.norm()==0 && (filter.GetTransition()-transition).norm()==0,
          "failed propagation is atomic");
    DynamicState imu,robot;
    check(!filter.Predict(huge,imu,robot),"high rate overflow rejected");
    check(filter.Predict(sample(1.1),imu,robot) && imu.v.norm()==0,"high rate recovers without cached invalid endpoint");
    check(filter.Predict(sample(1.1)) && filter.GetTime()==1.1,"main propagation recovers without invalid endpoint");
  }
  std::deque<IMUData> imus;double last_imu=-1;
  check(appendImu(sample(1),imus,last_imu) && appendImu(sample(1.1),imus,last_imu),"queue priming");
  check(!appendImu(sample(1.05),imus,last_imu) && !appendImu(sample(1.1),imus,last_imu) && imus.size()==2 && last_imu==1.1,
        "late/duplicate IMU preserves queued future samples");
  std::deque<LidarData> lidars{scan(1,1.05),scan(1.05,1.1),scan(1.05,1.1),scan(1.1,1.2)};
  double last_scan=-1;MeasureGroup group;
  check(synchronizeMeasurements(lidars,imus,last_scan,group) && group.imu.size()==2 && imus.size()==1 && imus.front().secs==1.1,
        "right bracket borrowed but retained");
  check(synchronizeMeasurements(lidars,imus,last_scan,group) && group.imu.size()==1 && imus.empty(),"exact endpoint consumed once");
  appendImu(sample(1.15),imus,last_imu);
  check(!synchronizeMeasurements(lidars,imus,last_scan,group) && lidars.size()==1,"duplicate scan dropped without consuming IMU");
  check(!synchronizeMeasurements(lidars,imus,last_scan,group) && imus.size()==1,"wait for right bracket retains left samples");
  appendImu(sample(1.25),imus,last_imu);
  check(synchronizeMeasurements(lidars,imus,last_scan,group) && imus.front().secs==1.25,"synchronization resumes on right bracket");
  livox_ros_driver::CustomMsg message;message.point_num=100;message.points.resize(10);
  LidarData converted;
  check(!convertLivoxScan(message,1,0,100,converted),"actual Livox path rejects declared length overflow");
  message.point_num=10;message.header.stamp.fromSec(1);
  for(size_t i=0;i<message.points.size();++i) {
    message.points[i].x=2;message.points[i].tag=0x10;message.points[i].offset_time=(9-i)*1000000;
  }
  check(convertLivoxScan(message,1,0,100,converted) && converted.pc->size()==10 && std::abs(converted.end_time-1.009)<1e-12,
        "actual Livox path uses maximum point time without reordering");
  check(!convertLivoxScan(message,0,0,100,converted),"invalid point stride cannot loop forever");
  std::vector<DynamicState> states(2);states[1].time=1;
  check(!deskewPoint(V3::Ones(),std::numeric_limits<double>::quiet_NaN(),states).allFinite(),"deskew rejects NaN before binary search");
  std::cout<<"deep temporal and input regression checks passed\n";
}
