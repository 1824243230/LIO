#include "lio/super_lio.h"
#include <stdexcept>
#include <iostream>
#include <limits>
using namespace LI2Sup;
using namespace BASIC;
void check(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
class Harness : public SuperLIO {
 public:
  void prepare(bool enhanced) {
    kf_=std::make_shared<ESKF>();
    ivox_=std::make_shared<OctVoxMapType>(OctVoxMapType::Options{0.5f,1000});
    scan_undistort_full_.reset(new PointCloudType());
    ds_undistort_.reset(new PointCloudType());
    VV3 history;
    for(int x=-16;x<=16;++x) for(int y=-16;y<=16;++y) {
      history.emplace_back(0.07*x,0.07*y,1.0);
      PointType p; p.x=0.07*x; p.y=0.07*y; p.z=1.02; p.intensity=history.size();
      scan_undistort_full_->push_back(p);
    }
    ivox_->insert(history);
    voxel_grid_fliter_.setLeafSize(0.3f);
    geometry_options_.enable_weak_selection=enhanced;
    geometry_options_.enable_adaptive_covariance=enhanced;
    geometry_options_.enable_bump_layer=enhanced;
    geometry_options_.enable_degeneracy=enhanced;
    geometry_options_.enable_bump_measurement=enhanced;
    geometry_options_.enable_informed_sampling=enhanced;
    if(enhanced) bump_map_=std::make_unique<geometry::BumpMap>(geometry_options_,0.5);
    DownSample();
  }
  void emulate_changed_sampling() {
    ds_undistort_.reset(new PointCloudType(*scan_undistort_full_));
    check(ds_undistort_->size()>original_sample_->size(),"test starts from distinct finer sample");
  }
  void remove_map() { bump_map_.reset(); }
  void fail_analysis() { geometry_options_.rotation_length_scale=std::numeric_limits<double>::quiet_NaN(); }
  void run() { Observe(); }
  auto state() const { return kf_->GetSysState(); }
  auto covariance() const { return kf_->GetCov(); }
  const auto& points() const { return *ds_undistort_; }
  const auto& diagnostics() const { return geometry_stats_; }
  int frame_count() const { return frame_num_; }
  void verify_restoration() const {
    check(ds_undistort_==original_sample_,"original sampling pointer restored");
    check(points_body_v3_.size()==original_sample_->size(),"body points rebuilt for original sample");
    for(size_t i=0;i<points_body_v3_.size();++i) {
      const auto& p=original_sample_->points[i];
      check((points_body_v3_[i]-V3(p.x,p.y,p.z)).norm()==0,"downstream map insertion uses restored points");
    }
  }
};
void verify_case(int failure) {
  Harness original,enhanced;
  original.prepare(false); enhanced.prepare(true);
  enhanced.emulate_changed_sampling();
  if(failure==1) enhanced.remove_map();
  if(failure==2) enhanced.fail_analysis();
  original.run(); enhanced.run();
  enhanced.verify_restoration();
  check(enhanced.diagnostics().fallback && enhanced.diagnostics().accepted==0 &&
        (failure!=0 || enhanced.diagnostics().rejected_map>0),"fallback preserves rejection diagnostics");
  check(enhanced.frame_count()==1 && original.frame_count()==1,"retry advances frame count once");
  check(enhanced.points().size()==original.points().size(),"baseline sample size restored");
  for(size_t i=0;i<original.points().size();++i)
    check((original.points()[i].getVector3fMap()-enhanced.points()[i].getVector3fMap()).norm()==0,
          "baseline sample identity and ordering restored");
  const auto a=original.state(),b=enhanced.state();
  check(a.p.norm()>0.001,"synthetic plane update has nonzero correction");
  check((a.p-b.p).norm()<1e-6 && (a.R.R_-b.R.R_).norm()<1e-6 &&
        (a.v-b.v).norm()<1e-6 && (a.bg-b.bg).norm()<1e-6 && (a.ba-b.ba).norm()<1e-6 &&
        (original.covariance()-enhanced.covariance()).norm()<1e-6,
        "full original-sampling replay equals original SuperLIO update");
}
int main() {
  for(int failure=0;failure<3;++failure) verify_case(failure);
  std::cout << "fallback tests passed: empty map, absent map, invalid analysis vs all-flags-off baseline\n";
}
