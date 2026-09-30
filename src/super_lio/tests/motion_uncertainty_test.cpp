#include "lio/motion_uncertainty.h"
#include "lio/deskew.h"
#include "lio/super_lio.h"
#include <iostream>
#include <stdexcept>
using namespace LI2Sup;
void check(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
ESKF makeFilter(bool smooth, const SysState& state = SysState()) {
  ESKF::Options o; o.smooth_motion_ = smooth;
  ESKF kf; kf.SetInitialConditions(o, state.bg, state.ba, 1, BASIC::V3::Zero());
  IMUData first; first.secs=0; first.acc=BASIC::V3(2,0,0); first.gyr=BASIC::V3(0,0,3);
  kf.Predict(first); kf.SetX(state); kf.SetObsTime(0.1);
  return kf;
}
class Harness : public SuperLIO {
public:
  void prepare(bool uncertainty, bool enhanced) {
    uncertainty_motion_=uncertainty;
    kf_=std::make_shared<ESKF>();
    ivox_=std::make_shared<OctVoxMapType>(OctVoxMapType::Options{0.5f,1000});
    scan_undistort_full_.reset(new BASIC::PointCloudType());
    ds_undistort_.reset(new BASIC::PointCloudType());
    BASIC::VV3 map;
    for(int x=-16;x<=16;++x) for(int y=-16;y<=16;++y) {
      map.emplace_back(0.07*x,0.07*y,1.0);
      BASIC::PointType p; p.x=0.07*x; p.y=0.07*y; p.z=1.02; p.intensity=7;
      scan_undistort_full_->push_back(p);
      motion_covariances_[{p.x,p.y,p.z}]=0.1*motion::M3::Identity();
    }
    ivox_->insert(map); voxel_grid_fliter_.setLeafSize(0.3f);
    geometry_options_.enable_degeneracy=enhanced;
    geometry_options_.enable_bump_layer=enhanced;
    geometry_options_.enable_bump_measurement=enhanced;
    geometry_options_.enable_informed_sampling=enhanced;
    if(enhanced) bump_map_=std::make_unique<geometry::BumpMap>(geometry_options_,0.5);
    DownSample();
    if (enhanced) ds_undistort_.reset(new BASIC::PointCloudType(*scan_undistort_full_));
  }
  void propagation() {
    smooth_motion_=true; uncertainty_motion_=true;
    kf_=std::make_shared<ESKF>(makeFilter(true));
    scan_undistort_full_.reset(new BASIC::PointCloudType());
    ds_undistort_.reset(new BASIC::PointCloudType());
    measures_.lidar.start_time=0; measures_.lidar.end_time=0.1;
    measures_.lidar.pc.reset(new pcl::PointCloud<PointXTZIT>());
    for(int i=0;i<3;++i) measures_.lidar.pc->emplace_back(3+i,2,1,10+i,0.05*i);
    IMUData imu; imu.secs=0.1; imu.acc=BASIC::V3(2,0,0); imu.gyr=BASIC::V3(0,0,3);
    measures_.imu.push_back(imu);
    Propagation_Undistort();
    check(motion_covariances_.size()==3,"full pipeline creates one covariance per distinct point");
    const auto& first=scan_undistort_full_->front(); const auto& last=scan_undistort_full_->back();
    check(motion_covariances_.at({first.x,first.y,first.z}).trace()>0,"early point has motion uncertainty");
    check(motion_covariances_.at({last.x,last.y,last.z}).norm()==0,"pipeline preserves end zero covariance");
    check(last.intensity==12 && first.intensity==10,"point attributes preserved");
    voxel_grid_fliter_.setLeafSize(0.2); DownSample();
    for(const auto& p:*ds_undistort_) check(motion_covariances_.count({p.x,p.y,p.z})==1,"voxel sample covariance association");
  }
  void run() { Observe(); check(observation_valid_, "weighted Observe succeeds"); }
  double variance() const { return kf_->GetCov()(5,5); }
  double z() const { return kf_->GetSysState().p.z(); }
  bool restored() const { return ds_undistort_==original_sample_; }
};
int main() {
  using namespace motion;
  M3 Iv,Ip;
  integrals(V3::Zero(),0.2,Iv,Ip);
  check((Iv-0.2*M3::Identity()).norm()<1e-14 && (Ip-0.02*M3::Identity()).norm()<1e-14,"zero angular rate limit");
  // Independent quadrature of continuously rotating acceleration.
  const V3 w(1,-2,3), a(2,-1,4); const double dt=0.2;
  integrals(w,dt,Iv,Ip); V3 v=V3::Zero(), p=V3::Zero();
  for(int i=0;i<10000;++i) {
    const double t=(i+0.5)*dt/10000;
    const V3 acc=Eigen::AngleAxisd(w.norm()*t,w.normalized())*a;
    v+=acc*(dt/10000); p+=acc*(dt-t)*(dt/10000);
  }
  check((v-Iv*a).norm()<1e-9 && (p-Ip*a).norm()<1e-9,"integrals match independent quadrature");
  IMUData imu; imu.secs=0.1; imu.acc=BASIC::V3(2,0,0); imu.gyr=BASIC::V3(0,0,3);
  auto kf=makeFilter(true); std::vector<DynamicState> states{kf.GetDynamicState()};
  check(kf.Predict(imu),"smooth predict"); states.push_back(kf.GetDynamicState());
  integrals(V3(0,0,3),0.1,Iv,Ip);
  check((kf.GetSysState().p.cast<double>()-Ip*V3(2,0,0)).norm()<1e-8,"nominal rotating acceleration");
  const BASIC::V3 point(3,2,1);
  check((deskewPoint(point,0.1-1e-9,states,true)-point).norm()<1e-6,"continuous deskew endpoint");
  auto half=makeFilter(true); half.SetObsTime(0.05); check(half.Predict(imu),"clipped endpoint");
  check(half.Predict(imu)==false,"duplicate clipped interval rejected");
  half.SetObsTime(0.1); check(half.Predict(imu),"remaining interval");
  check((half.GetSysState().p-kf.GetSysState().p).norm()<1e-7 &&
        (half.GetSysState().v-kf.GetSysState().v).norm()<1e-7,"nominal propagation partition invariance");
  // Finite-difference the complete nominal state against transition columns.
  const auto F=kf.GetTransition();
  constexpr float eps=0.001f;
  for(int j=0;j<15;++j) {
    SysState plus,minus;
    if(j<3) { BASIC::V3 delta=BASIC::V3::Zero(); delta[j]=eps;
      plus.R=BASIC::SO3::Exp(delta); minus.R=BASIC::SO3::Exp(-delta); }
    else {
      BASIC::V3* xp=j<6?&plus.p:j<9?&plus.v:j<12?&plus.bg:&plus.ba;
      BASIC::V3* xm=j<6?&minus.p:j<9?&minus.v:j<12?&minus.bg:&minus.ba;
      (*xp)[j%3]=eps; (*xm)[j%3]=-eps;
    }
    auto kp=makeFilter(true,plus), km=makeFilter(true,minus); kp.Predict(imu); km.Predict(imu);
    const auto sp=kp.GetSysState(), sm=km.GetSysState();
    Eigen::Matrix<double,9,1> derivative;
    derivative.head<3>()=(kf.GetSysState().R.inverse()*sp.R).log_vee().cast<double>()-
                        (kf.GetSysState().R.inverse()*sm.R).log_vee().cast<double>();
    derivative.segment<3>(3)=(sp.p-sm.p).cast<double>();
    derivative.tail<3>()=(sp.v-sm.v).cast<double>(); derivative/=2*eps;
    check((derivative-F.col(j).head<9>().cast<double>()).norm()<1e-4,"transition matches nominal finite difference");
  }
  std::vector<Knot> knots(3);
  for(int i=0;i<3;++i) { knots[i].time=i*0.1; knots[i].covariance=M18::Identity(); }
  // Shared, completely uncertain absolute pose cancels if no motion/no noise.
  relativeCovariances(knots);
  check(knots[0].relative.norm()==0 && knots[1].relative.norm()==0,"shared pose uncertainty cancels");
  knots[1].transition.block<3,3>(3,6)=0.1*M3::Identity();
  knots[2].transition=knots[1].transition;
  knots[1].noise.topLeftCorner<3,3>()=0.001*M3::Identity();
  knots[2].noise=knots[1].noise;
  for(int i=1;i<3;++i) knots[i].covariance=knots[i].transition*knots[i-1].covariance*knots[i].transition.transpose()+knots[i].noise;
  relativeCovariances(knots);
  check(std::abs(knots[0].relative(3,3)-0.04)<1e-12,"shared velocity uncertainty contributes dt squared");
  check(std::abs(knots[0].relative(0,0)-0.002)<1e-12,"future process noise accumulates");
  check(knots.back().relative.norm()==0 && pointCovariance(V3(2,0,0),0.2,knots).norm()==0,"end covariance exactly zero");
  for(double t : {-0.1,0.0,0.07,0.19}) {
    const M3 C=pointCovariance(V3(3,2,1),t,knots);
    check(C.allFinite() && Eigen::SelfAdjointEigenSolver<M3>(C).eigenvalues().minCoeff()>-1e-12,"interpolated covariance PSD");
  }
  check(pointCovariance(V3(30,0,0),0,knots).trace()>pointCovariance(V3(3,0,0),0,knots).trace(),"range amplifies rotational uncertainty");
  // Compare the O(N) backward recursion with an independently assembled joint covariance.
  std::vector<Knot> joint_knots(4);
  joint_knots[0].covariance=0.02*M18::Identity();
  for(int i=0;i<4;++i) {
    auto& k=joint_knots[i]; k.time=i*0.01;
    k.R=Eigen::AngleAxisd(0.1*i,V3(1,2,3).normalized()).toRotationMatrix();
    k.p=V3(0.1*i,-0.02*i,0.01*i*i);
    if(i) {
      k.transition=kf.GetTransition().cast<double>();
      k.noise=kf.GetProcessNoise().cast<double>();
      k.covariance=k.transition*joint_knots[i-1].covariance*k.transition.transpose()+k.noise;
    }
  }
  relativeCovariances(joint_knots);
  for(int i=0;i<3;++i) {
    const auto& start=joint_knots[i]; const auto& end=joint_knots.back();
    M18 phi=M18::Identity();
    for(int j=i+1;j<4;++j) phi=(joint_knots[j].transition*phi).eval();
    const M18 cross=start.covariance*phi.transpose();
    Eigen::Matrix<double,6,18> A=Eigen::Matrix<double,6,18>::Zero(), B=A;
    A.block<3,3>(0,0)=end.R.transpose()*start.R; A.block<3,3>(3,3)=end.R.transpose();
    B.block<3,3>(0,0)=-M3::Identity(); B.block<3,3>(3,3)=-end.R.transpose();
    B.block<3,3>(3,0)=hat(end.R.transpose()*(start.p-end.p));
    const M6 explicit_joint=A*start.covariance*A.transpose()+B*end.covariance*B.transpose()+
                            A*cross*B.transpose()+B*cross.transpose()*A.transpose();
    check((explicit_joint-start.relative).norm()<1e-12,"recursion equals full joint cross covariance");
  }
  std::vector<geometry::V6> J(1,geometry::V6::Ones());
  std::vector<double> errors{0.2}, precisions{20}; std::vector<unsigned char> valid{1};
  geometry::Candidate candidate; candidate.index=0; candidate.J=geometry::V6::Unit(2); candidate.variance=0.04; candidate.residual=0.1;
  std::vector<geometry::Candidate> candidates{candidate};
  M6 info=20*J[0]*J[0].transpose(); geometry::V6 rhs=-20*J[0]*errors[0];
  geometry::applyExclusiveBumps(J,errors,valid,candidates,info,rhs,1000,&precisions);
  check((info-candidate.J*candidate.J.transpose()/candidate.variance).norm()<1e-10 &&
        (rhs+candidate.J*candidate.residual/candidate.variance).norm()<1e-10,"bump removes actual weighted plane row");
  Harness pipeline; pipeline.propagation();
  Harness base,weighted,fallback;
  base.prepare(false,false); weighted.prepare(true,false); fallback.prepare(true,true);
  base.run(); weighted.run(); fallback.run();
  check(weighted.variance()>base.variance(),"UAMC does not retain overconfident posterior");
  check(fallback.restored() && std::abs(weighted.z()-fallback.z())<1e-6 &&
        std::abs(weighted.variance()-fallback.variance())<1e-6,"fallback retains covariance association and weights");
  std::cout << "motion uncertainty tests passed\n";
}
