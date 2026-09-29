#include "lio/informed_sampling.h"
#include "OctVoxMap/VoxelGridFilter.h"
#include <Eigen/Geometry>
#include <set>
#include <stdexcept>
#include <iostream>
using namespace LI2Sup;
using namespace LI2Sup::geometry;
void check(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
BASIC::CloudPtr sample(const BASIC::CloudPtr& input,double resolution) {
  BASIC::CloudPtr out(new BASIC::PointCloudType());
  VoxelGridClosest<BASIC::PointType> filter;
  filter.setLeafSize(resolution); filter.setInputCloud(input); filter.filter(out); return out;
}
void same(const BASIC::CloudPtr& a,const BASIC::CloudPtr& b) {
  check(a->size()==b->size() && a->header.stamp==b->header.stamp,"baseline size/header unchanged");
  for(size_t i=0;i<a->size();++i)
    check((a->points[i].getVector3fMap()-b->points[i].getVector3fMap()).norm()==0 &&
          a->points[i].intensity==b->points[i].intensity,"baseline point identity/order unchanged");
}
int main() {
  Options o; o.enable_informed_sampling=true; o.bump_history_points=4096;
  o.bump_resolution=0.04; o.bump_min_roughness=0.0001;
  o.mid_threshold=0.001; o.informed_top_k_voxels=1;
  BumpMap map(o,0.5);
  std::vector<V3> history;
  BASIC::CloudPtr input(new BASIC::PointCloudType()); input->header.stamp=123456;
  const M3 R=Eigen::AngleAxisd(0.5,V3::UnitZ()).toRotationMatrix();
  const V3 t(10,0,0);
  for(int region=0;region<2;++region) for(int x=0;x<40;++x) for(int y=0;y<40;++y) {
    double px=0.05+x*0.01, py=0.05+y*0.01;
    V3 world(px+region*0.5,py,0.2+(region?0.022:0.008)*std::sin(40*px)*std::cos(30*py));
    history.push_back(world);
    V3 body=R.transpose()*(world-t);
    BASIC::PointType p; p.x=body.x(); p.y=body.y(); p.z=body.z(); p.intensity=input->size();
    input->push_back(p);
  }
  map.insert(history,V3::Zero()); map.insert(history,V3::Zero());
  check(map.mid(V3(0.7,0.2,0.2))>map.mid(V3(0.2,0.2,0.2)),"synthetic high/low MID regions");
  auto baseline=sample(input,0.3);
  auto output=BASIC::CloudPtr(new BASIC::PointCloudType(*baseline));
  auto original_ptr=output.get();
  DegeneracyResult d; d.valid=true; d.status=Status::NORMAL;
  auto result=informedSample(input,R,t,&map,d,o,output);
  check(!result.applied && output.get()==original_ptr,"NORMAL does not replace original sample"); same(output,baseline);
  d.status=Status::WEAK;
  Options disabled=o; disabled.enable_informed_sampling=false;
  check(!informedSample(input,R,t,&map,d,disabled,output).applied,"feature off preserves baseline"); same(output,baseline);
  check(!informedSample(input,R,t,nullptr,d,o,output).applied,"missing map preserves baseline"); same(output,baseline);
  check(!informedSample(input,R,t,&map,DegeneracyResult{},o,output).applied,"invalid analysis preserves baseline");
  check(!informedSample(input,M3::Identity(),V3::Zero(),&map,d,o,output).applied,
        "queries depend on prior transform"); same(output,baseline);
  result=informedSample(input,R,t,&map,d,o,output);
  check(result.applied && result.selected_voxels==1 && result.fine_points>0 && result.coarse_points>0,
        "WEAK splits fine/coarse with configured Top-1");
  check(output->size()==result.fine_points+result.coarse_points && output->header.stamp==input->header.stamp,
        "combined cloud count and timestamp");
  std::set<float> ids;
  for(const auto& p:*output) check(ids.insert(p.intensity).second,"fine/coarse output has no duplicate input point");
  auto fine=sample(input,o.fine_resolution);
  const auto high_key=map.key(V3(0.7,0.2,0.2));
  size_t expected_fine=0;
  for(const auto& p:*fine) if(map.key(R*V3(p.x,p.y,p.z)+t)==high_key) {
    ++expected_fine; check(ids.count(p.intensity)==1,"selected region retains each fine representative");
  }
  check(expected_fine==result.fine_points,"highest MID voxel receives fine resolution");
  o.informed_top_k_voxels=2; d.status=Status::DEGENERATE;
  result=informedSample(input,R,t,&map,d,o,output);
  check(result.applied && result.selected_voxels==2 && result.fine_points>expected_fine,"configurable Top-2 in DEGENERATE");
  o.mid_threshold=100;
  result=informedSample(input,R,t,&map,d,o,output);
  check(result.applied && result.selected_voxels==0 && result.fine_points==0,"available low-MID regions use coarse sampling");
  same(output,sample(fine,o.coarse_resolution));
  o.coarse_resolution=0.8;
  result=informedSample(input,R,t,&map,d,o,output);
  same(output,sample(fine,o.coarse_resolution));
  check(result.coarse_points==output->size(),"coarse resolution parameter applied");
  std::cout << "sampling tests passed: baseline identity, prior transform, MID, dual resolution, configurable Top-K, fallback\n";
}
