#include "lio/plane_fit.h"
#include "lio/deskew.h"
#include "lio/ESKF.h"
#include "lio/scan_validation.h"
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace BASIC;
using namespace LI2Sup;
void check(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }

int main() {
  LidarData scan;
  scan.start_time=10;
  scan.pc.reset(new pcl::PointCloud<PointXTZIT>());
  check(!finalizeLidarScan(scan), "fully filtered scan rejected without back access");
  scan.pc->emplace_back(1,2,3,1,0.1);
  scan.pc->emplace_back(1,2,3,1,0.02);
  scan.pc->emplace_back(1,2,3,1,-0.1);
  check(finalizeLidarScan(scan) && scan.pc->size()==3 && std::abs(scan.end_time-10.1)<1e-7 &&
        std::abs(scan.start_time-9.9)<1e-7 && std::abs(scan.start_time+scan.pc->points[1].offset_time-10.02)<1e-7,
        "signed offsets preserve all absolute timestamps independent of point order");
  std::array<V3, 5> points = {V3(-1,-1,0), V3(-1,1,0), V3(1,-1,0), V3(1,1,0), V3(0,0,0)};
  std::array<double,4> plane;
  check(calc_plane_coeff(5, points, plane), "plane through origin is observable");
  check(std::abs(plane[2]) > 0.999 && std::abs(plane[3]) < 1e-10, "origin plane normal");
  for (auto& point : points) point += V3(10000,20000,30000);
  check(calc_plane_coeff(4, points, plane), "translated plane with four neighbors");
  check(std::abs(std::abs(plane[3])-30000) < 1e-6, "translated plane offset");
  for (int i=0;i<5;++i) points[i]=V3(i,0,0);
  check(!calc_plane_coeff(5, points, plane), "collinear neighbors must not invent a normal");
  points[0].x()=std::numeric_limits<scalar>::quiet_NaN();
  check(!calc_plane_coeff(5, points, plane), "reject nonfinite neighbors");

  std::vector<DynamicState> states;
  for (int i=0;i<=10;++i) states.emplace_back(i*0.1, M3::Identity(), V3(i*0.1,0,0),
                                            V3(1,0,0), V3::Zero(), V3::Zero());
  const V3 point(2,0,0);
  for (double t : {-0.1, 0.0, 0.1, 0.55, 1.0, 1.1}) {
    const V3 actual=deskewPoint(point,t,states);
    check((actual-V3(1+std::clamp(t,0.0,1.0),0,0)).norm()<1e-6,
          "deskew constant velocity including both boundaries");
  }
  check((deskewPoint(point,0,{states.front()})-point).norm()==0, "single state deskew");

  ESKF filter;
  filter.SetInitialConditions(ESKF::Options(),V3::Zero(),V3::Zero(),1.0,V3(0,0,-g_gravity_norm));
  IMUData imu; imu.secs=1.0; imu.acc=V3(1,0,g_gravity_norm);
  filter.SetLastObsTime(1.0); filter.SetObsTime(1.05);
  check(!filter.Predict(imu), "prime IMU history");
  imu.secs=1.1;
  check(filter.Predict(imu), "integrate partial IMU interval");
  check(std::abs(filter.GetTime()-1.05)<1e-12, "state timestamp reaches exact scan end");
  check(std::abs(filter.GetSysState().v.x()-0.05)<1e-6 &&
        std::abs(filter.GetSysState().p.x()-0.00125)<1e-6, "clipped propagation duration");
  filter.SetLastObsTime(1.05); filter.SetObsTime(1.1);
  check(filter.Predict(imu), "retained right bracket integrates next frame remainder");
  check(std::abs(filter.GetSysState().v.x()-0.1)<1e-6 &&
        std::abs(filter.GetSysState().p.x()-0.005)<1e-6, "no lost or double integrated time");
  check(!filter.Predict(imu), "duplicate sample rejected");
  imu.secs=1.09;
  check(!filter.Predict(imu), "out of order sample rejected");

  ESKF weak;
  const auto weak_plane=[](const ESKF::KFState& s,M6& H,V6& b) {
    H.setZero(); b.setZero(); H(5,5)=1000; b[5]=1000*(0.1-s.pose.t_.z());
  };
  check(weak.UpdateObserve(weak_plane), "rank one measurement remains solvable with prior");
  check(std::abs(weak.GetCov()(3,3)-1)<1e-6, "unobserved translation retains prior uncertainty");
  check(weak.GetCov()(5,5)<0.002, "observed translation gains information");
  const ESKF before=weak;
  int calls=0;
  check(!weak.UpdateObserve([&](const ESKF::KFState& s,M6& H,V6& b) {
    weak_plane(s,H,b); b[5]=10;
    if (++calls==2) b[0]=std::numeric_limits<scalar>::quiet_NaN();
  }), "nonfinite second iteration rejected");
  check(calls==2 && (weak.GetSysState().p-before.GetSysState().p).norm()==0 &&
        (weak.GetCov()-before.GetCov()).norm()==0 &&
        (weak.GetGravity()-before.GetGravity()).norm()==0, "later failure rolls back full frame");
  check(!weak.UpdateObserve([](const ESKF::KFState&,M6& H,V6& b) {H.setZero();b.setZero();}),
        "empty observations rejected");
  ESKF invalid;
  auto covariance=invalid.GetCov(); covariance(0,0)=-1; invalid.SetCov(covariance);
  check(!invalid.UpdateObserve(weak_plane), "non-SPD prior rejected");
  std::cout << "LIO numerics passed: plane support, deskew boundaries, IMU bracket, weak prior, atomic rollback\n";
}
