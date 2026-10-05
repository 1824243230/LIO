#include "lio/geometry.h"
#include <Eigen/Geometry>
#include <cassert>
#include <pcl/io/pcd_io.h>
#include "OctVoxMap/OctVoxMap.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>
using namespace LI2Sup::geometry;
void check(bool ok, const char* message) { if(!ok) throw std::runtime_error(message); }
int main() {
  Options o; o.bump_history_points=4096; o.bump_resolution=0.04; o.bump_min_roughness=0.0001;
  BumpMap map(o,0.5);
  V3 probe(0.2,0.2,0.2);
  map.insert({probe},V3::Zero());
  check(map.stats(probe) && map.stats(probe)->count==1 && !map.stats(probe)->plane_valid,
        "sparse geometry statistics");
  check(!map.layer(probe), "sparse voxel does not allocate layer");
  std::vector<V3> flat, ripple;
  V3 sum=V3::Zero(); M3 outer=M3::Zero();
  for(int x=0;x<40;++x) for(int y=0;y<40;++y) {
    double px=0.05+x*0.01, py=0.05+y*0.01;
    flat.emplace_back(px,py,0.2);
    ripple.emplace_back(px,py,0.2+0.015*std::sin(40*px)*std::cos(30*py));
    sum+=ripple.back(); outer+=ripple.back()*ripple.back().transpose();
  }
  BumpMap plane(o,0.5);
  plane.insert(flat,V3::Zero()); plane.insert(flat,V3::Zero());
  check(plane.stats(probe)->plane_valid && !plane.layer(probe), "valid pure plane without image allocation");
  BumpMap surface(o,0.5);
  surface.insert(ripple,V3::Zero());
  const auto* stats=surface.stats(probe);
  check(stats && stats->count==ripple.size() && (stats->sum-sum).norm()<1e-10 &&
        (stats->sum_outer-outer).norm()<1e-10 && (stats->centroid-sum/ripple.size()).norm()<1e-12,
        "stored sufficient statistics match history");
  check(stats->plane_valid && std::abs(stats->normal.norm()-1)<1e-12 &&
        stats->eigenvalues[0]<=stats->eigenvalues[1] && stats->eigenvalues[1]<=stats->eigenvalues[2],
        "valid geometry spectrum and normal");
  check(!surface.layer(probe), "new plane waits for stable normal");
  surface.insert(ripple,V3::Zero());
  const auto* layer=surface.layer(probe);
  check(layer && layer->valid && layer->resolution==o.bump_resolution &&
        layer->image.size()==size_t(layer->width*layer->height) && layer->image.size()==layer->weight.size(),
        "stable surface lazily allocates complete layer");
  check((layer->R_CG*surface.stats(probe)->centroid+layer->t_CG).norm()<1e-10,
        "chart translation matches world-to-local convention");
  double total_weight=0;
  for(double w:layer->weight) total_weight+=w;
  check(std::abs(total_weight-2*ripple.size())<1e-8, "range weights accumulated from real samples");
  // Map preparation has no dependency on degeneracy or measurement activation.
  // The same prepared surface must remain unused while the analyzer reports NORMAL.
  Surface prepared;
  check(surface.query(V3(0.2213,0.2337,0.2),prepared), "prepared map query before activation");
  check(prepared.has_voxel_key && prepared.voxel_key==surface.key(V3(0.2213,0.2337,0.2)),
        "prepared surface carries its map voxel");
  Analyzer analyzer;
  auto normal=analyzer.analyze(M6::Identity(),o,true);
  Options gates=o;
  gates.bump_min_mid=1e-8; gates.bump_min_gradient=1e-8;
  gates.bump_min_pixel_confidence=1e-8; gates.bump_min_weak_score=1e-8;
  gates.weak_direction_score_threshold=0;
  Candidate measurement; Diagnostics diagnostics;
  check(!candidate(prepared,probe,M3::Identity(),normal,gates,measurement,diagnostics),
        "NORMAL rejects measurement even when historical layer exists");
  M6 weak_info=M6::Zero(); weak_info(0,0)=1;
  auto weak=analyzer.analyze(weak_info,o,true);
  check(candidate(prepared,probe,M3::Identity(),weak,gates,measurement,diagnostics),
        "prepared historical layer can activate after degeneracy");
  check(measurement.has_voxel_key && measurement.voxel_key==prepared.voxel_key,
        "map voxel identity reaches accepted candidate");
  check(!candidate(prepared,probe,M3::Identity(),DegeneracyResult{},gates,measurement,diagnostics),
        "invalid analysis never activates prepared layer");

  Options few=o; few.bump_min_points=2000;
  BumpMap insufficient(few,0.5);
  insufficient.insert(ripple,V3::Zero());
  check(insufficient.stats(probe)->count==1600 && !insufficient.layer(probe), "parameterized point-count gate");
  Options rough=o; rough.bump_min_roughness=0.1;
  BumpMap smooth(rough,0.5);
  smooth.insert(ripple,V3::Zero()); smooth.insert(ripple,V3::Zero());
  check(smooth.stats(probe)->plane_valid && !smooth.layer(probe), "parameterized roughness gate");
  std::vector<V3> volume;
  for(int x=0;x<5;++x) for(int y=0;y<5;++y) for(int z=0;z<5;++z)
    volume.emplace_back(0.05+0.08*x,0.05+0.08*y,0.05+0.08*z);
  BumpMap nonplanar(o,0.5);
  nonplanar.insert(volume,V3::Zero()); nonplanar.insert(volume,V3::Zero());
  check(!nonplanar.stats(probe)->plane_valid && !nonplanar.layer(probe), "nonplanar geometry cannot allocate layer");
  Surface invalid_plane_result;
  check(!nonplanar.query(probe,invalid_plane_result), "invalid plane rejects map query");

  Options angle_options=o; angle_options.bump_history_points=ripple.size();
  BumpMap changing(angle_options,0.5);
  changing.insert(ripple,V3::Zero()); changing.insert(ripple,V3::Zero());
  check(changing.layer(probe)!=nullptr, "stable chart before normal change");
  std::vector<V3> tilted;
  const M3 tilt=Eigen::AngleAxisd(20*std::acos(-1.0)/180,V3::UnitX()).toRotationMatrix();
  const V3 center(0.25,0.25,0.2);
  for(const auto& p:ripple) tilted.push_back(center+tilt*(p-center));
  const BumpLayer old_chart = *changing.layer(probe);
  changing.insert(tilted,V3::Zero());
  const auto* moved = changing.layer(probe);
  check(changing.stats(probe)->plane_valid && moved && moved->valid, "normal change reprojects existing image");
  check((moved->R_CG-old_chart.R_CG).norm()>0.1, "reprojection changes chart rotation");
  double expected_weight=0, expected_moment=0, actual_weight=0, actual_moment=0;
  for(int v=0;v<old_chart.height;++v) for(int u=0;u<old_chart.width;++u) {
    const int i=v*old_chart.width+u;
    if(old_chart.weight[i]<=0) continue;
    const V3 local((u-(old_chart.width-1)*0.5)*old_chart.resolution,
                   (v-(old_chart.height-1)*0.5)*old_chart.resolution,old_chart.image[i]);
    const V3 world=old_chart.R_CG.transpose()*(local-old_chart.t_CG);
    expected_weight+=old_chart.weight[i];
    expected_moment+=old_chart.weight[i]*(moved->R_CG*world+moved->t_CG).z();
  }
  for(const auto& p:tilted) {
    const double w=std::min(o.bump_range_weight_max,1.0/(p.norm()+o.bump_range_epsilon));
    expected_weight+=w; expected_moment+=w*(moved->R_CG*p+moved->t_CG).z();
  }
  for(size_t i=0;i<moved->weight.size();++i) {
    actual_weight+=moved->weight[i]; actual_moment+=moved->weight[i]*moved->image[i];
  }
  check(std::abs(actual_weight-expected_weight)<1e-8, "reprojection conserves old confidence without double counting");
  check(std::abs(actual_moment-expected_moment)<1e-8, "reprojection transforms old pixel heights before fusion");

  // Unequal ranges exercise both the configurable cap and inverse-range weighting.
  Options weighted=o; weighted.bump_range_weight_max=0.4; weighted.bump_range_epsilon=2;
  BumpMap incremental(weighted,0.5);
  const std::vector<V3> sensors={V3::Zero(),V3(3,0,0),V3(0,4,0)};
  for(const auto& sensor:sensors) incremental.insert(ripple,sensor);
  const auto* result=incremental.layer(probe);
  check(result && result->valid, "weighted image available");
  std::vector<double> numerator(result->image.size(),0), denominator(result->image.size(),0);
  for(const auto& sensor:sensors) for(const auto& p:ripple) {
    const V3 local=result->R_CG*p+result->t_CG;
    const int u=int(std::round(local.x()/result->resolution+(result->width-1)*0.5));
    const int v=int(std::round(local.y()/result->resolution+(result->height-1)*0.5));
    check(u>=0 && v>=0 && u<result->width && v<result->height, "reference sample in image");
    const double w=std::min(weighted.bump_range_weight_max,1.0/((p-sensor).norm()+weighted.bump_range_epsilon));
    numerator[v*result->width+u]+=w*local.z(); denominator[v*result->width+u]+=w;
  }
  for(size_t i=0;i<denominator.size();++i) {
    check(std::abs(result->weight[i]-denominator[i])<1e-10, "incremental weight equals independent full sum");
    if(denominator[i]>0) check(std::abs(result->image[i]-numerator[i]/denominator[i])<1e-10,
                               "incremental mean equals independent weighted reference");
  }
  const auto weights_before=result->weight;
  incremental.insert({},V3::Zero());
  check(incremental.layer(probe)->weight==weights_before, "empty update does not replay history");
  // A bounded history replacing all ripples by a flat surface must release its old image.
  for(int i=0;i<3;++i) surface.insert(flat,V3::Zero());
  check(surface.stats(probe)->count==size_t(o.bump_history_points) &&
        surface.stats(probe)->plane_valid && !surface.layer(probe), "history replacement releases obsolete layer");
  Options limited=o; limited.bump_capacity=1;
  BumpMap evicted(limited,0.5);
  evicted.insert(ripple,V3::Zero()); evicted.insert(ripple,V3::Zero());
  check(evicted.layer(probe)!=nullptr, "eviction test starts with allocated layer");
  evicted.insert({V3(2,2,2)},V3::Zero());
  check(evicted.size()==1 && !evicted.stats(probe) && !evicted.layer(probe), "eviction removes geometry and image");

  // Identical OctVox insertion/query streams, with a BumpMap populated alongside only one.
  using OctMap=LI2Sup::OctVoxMap<V3,double>;
  OctMap baseline(OctMap::Options{0.5f,100}), hybrid(OctMap::Options{0.5f,100});
  OctMap::Points points;
  std::vector<V3> history;
  for(int x=-2;x<=2;++x) for(int y=-2;y<=2;++y) for(int z=-2;z<=2;++z) {
    V3 p(0.19*x,0.19*y,0.19*z); points.push_back(p); history.push_back(p);
  }
  BumpMap extension(o,0.5);
  for(int frame=0;frame<3;++frame) {
    baseline.insert(points); hybrid.insert(points); extension.insert(history,V3::Zero());
    baseline.reset_max_group(); hybrid.reset_max_group();
    for(const auto& p:points) {
      LI2Sup::KNNHeap<5,V3> a,b; baseline.getTopK(p,a); hybrid.getTopK(p,b);
      check(a.count==b.count && a.count>0, "HKNN count unchanged");
      for(int i=0;i<a.count;++i) check((a.points_[i]-b.points_[i]).norm()==0 && a.dist2_[i]==b.dist2_[i],
                                     "HKNN representatives and distances unchanged");
    }
  }
  std::cout << "map tests passed: statistics, lazy allocation, bounded history, eviction, HKNN invariance\n";
}
