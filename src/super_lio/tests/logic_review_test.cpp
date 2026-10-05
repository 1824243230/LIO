#include "lio/ESKF.h"
#include "lio/imu_initialization.h"
#include "lio/motion_uncertainty.h"
#include "OctVoxMap/VoxelGridFilter.h"
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace LI2Sup;
using namespace BASIC;
void check(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
int main() {
  VoxelGridClosest<PointType> filter;
  filter.setLeafSize(1);
  CloudPtr input(new PointCloudType()), output;
  auto add=[&](float x,float y,float z) { PointType p; p.x=x;p.y=y;p.z=z;p.intensity=7;input->push_back(p); };
  // 旧拼接键：前两点 15 bit 跨轴相撞，后两点负数转 unsigned 后相撞。
  add(31768,-1000,0); add(-1000,-999,0); add(-1001,0,0); add(-1001,2,3);
  add(31768.2f,-1000,0); // 同一体素应保留更靠近中心的原点。
  add(std::numeric_limits<float>::quiet_NaN(),0,0); add(1e30f,0,0);
  filter.setInputCloud(input); filter.filter(output);
  check(output && output->size()==4,"distinct signed/large voxels must not merge");
  check(output->front().x==31768 && output->front().intensity==7,"select actual center-nearest sample");
  bool threw=false; try { filter.setLeafSize(0); } catch(const std::invalid_argument&) { threw=true; }
  check(threw,"invalid voxel resolution rejected");

  ImuInitialization init,other;
  IMUData a; a.secs=1; a.acc=V3(2,0,0); a.gyr=V3(1,0,0);
  init.add(a); init.add(a);
  IMUData b=a; b.secs=2; b.acc=V3(4,0,0); init.add(b); init.add(a);
  check(init.count==2 && init.mean_acc.x()==3 && other.count==0,"borrowed brackets deduplicated and instances isolated");

  g_gravity_norm=9.81;
  // 建图和重定位现在共用该入口；验证正常时间锚点、偏置保留及首段传播。
  ESKF initialized;
  const V3 bg(0.1f,0.2f,0.3f), ba(0.01f,0.02f,0.03f);
  initialized.SetInitialConditions(ESKF::Options(),bg,ba,1,V3(0,0,-9.81));
  IMUData anchor; anchor.secs=1234.5;anchor.acc=ba+V3(0,0,9.81);anchor.gyr=bg;
  const SE3 initial_pose(SO3::Exp(V3(0,0,0.3f)),V3(1,2,3));
  initialized.SetPoseAtImu(initial_pose,anchor);
  check(initialized.GetTime()==anchor.secs &&
        (initialized.GetSysState().p-initial_pose.t_).norm()==0 &&
        (initialized.GetSysState().R.R_-initial_pose.R_).norm()<1e-7 &&
        (initialized.GetSysState().bg-bg).norm()==0 && (initialized.GetSysState().ba-ba).norm()==0,
        "initial pose handoff uses matching IMU time and preserves initialized biases");
  auto next=anchor; next.secs+=0.01;initialized.SetObsTime(next.secs);
  check(initialized.Predict(next) && (initialized.GetSysState().p-initial_pose.t_).norm()<1e-6 &&
        initialized.GetSysState().v.norm()<1e-6,"first prediction after registration uses actual IMU anchor");
  auto invalid_state=initialized.GetSysState();invalid_state.timestamp=-1;
  bool mismatch=false;
  try { initialized.SetX(invalid_state,anchor); } catch(const std::invalid_argument&) { mismatch=true; }
  check(mismatch && initialized.GetTime()==next.secs,"timestamp validation is retained, not bypassed for relocation");

  ESKF kf; ESKF::Options options;
  kf.SetInitialConditions(options,V3::Zero(),V3::Zero(),1,V3(0,0,-9.81));
  SysState state; state.timestamp=1;
  a.secs=1; a.acc=V3(0,0,9.81); a.gyr.setZero();
  kf.SetX(state,a); kf.init_=true;
  kf.SetObsTime(1.05);
  b=a; b.secs=1.1; b.acc.x()=4;
  check(kf.Predict(b),"scan-end partial prediction");
  DynamicState out,robot;
  check(kf.Predict(b,out,robot),"high-rate prediction before LiDAR update");
  check(kf.UpdateObserve([](const ESKF::KFState&,M6& H,V6& r) { H=M6::Identity();r.setZero(); }),"zero-residual correction");
  const float velocity=kf.GetSysState().v.x();
  IMUData c=b; c.secs=1.15; c.acc.x()=6;
  check(kf.Predict(c,out,robot),"high-rate prediction after correction");
  // 左端必须为 t=1.05 时插值的 ax=2，不能用旧高频末端 t=1.1 的 ax=4。
  check(std::abs(out.v.x()-(velocity+0.4f))<1e-6,"correction resets both IMU time and measurement");
  check(!kf.Predict(c,out,robot),"duplicate high-rate sample rejected");
  c.secs=std::numeric_limits<double>::quiet_NaN();
  check(!kf.Predict(c,out,robot),"nonfinite high-rate time rejected");

  ESKF constant; constant.SetInitialConditions(options,V3::Zero(),V3::Zero(),1,V3(0,0,-9.81));
  a.acc.x()=4; constant.SetX(state,a); constant.SetObsTime(1.1);
  b=a;b.secs=1.1;
  check(constant.Predict(b) && std::abs(constant.GetSysState().v.x()-0.4f)<1e-6,
        "first prediction uses measured initialization anchor rather than zero IMU");
  // 在零先验下，白色加速度测量噪声应直接产生位置方差及位置-速度互协方差。
  ESKF noise_filter; options.acce_var_=0.02;
  noise_filter.SetInitialConditions(options,V3::Zero(),V3::Zero(),1,V3(0,0,-9.81));
  noise_filter.SetX(state,a);noise_filter.SetCov(ESKF::COV::Zero());noise_filter.SetObsTime(1.1);
  check(noise_filter.Predict(b),"noise-only propagation");
  const auto covariance=noise_filter.GetCov();
  check(std::abs(covariance(3,3)-0.25*0.0001*0.02)<1e-10 &&
        std::abs(covariance(3,6)-0.5*0.001*0.02)<1e-10,
        "position and velocity share the same acceleration noise");
  auto changed=constant.GetSysState();changed.timestamp=3;constant.SetX(changed);constant.SetObsTime(3.1);
  b.secs=3.1;
  check(!constant.Predict(b),"state time reset cannot reuse a measurement from another epoch");

  std::vector<motion::Knot> knots(2);knots[1].time=1;
  check(!motion::pointCovariance(Eigen::Vector3d::Ones(),std::numeric_limits<double>::quiet_NaN(),knots).allFinite(),
        "NaN point time rejected before binary search");
  check(!motion::pointCovariance(Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity()),1,knots).allFinite(),
        "invalid point cannot masquerade as zero endpoint covariance");
  std::cout<<"logic review regressions passed\n";
}
