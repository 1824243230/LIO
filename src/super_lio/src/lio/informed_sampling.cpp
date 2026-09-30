#include "lio/informed_sampling.h"
#include "OctVoxMap/VoxelGridFilter.h"
#include <algorithm>
#include <set>
namespace LI2Sup { namespace geometry {
SamplingResult informedSample(const BASIC::CloudPtr& input, const M3& R, const V3& t,
                              const BumpMap* map, const DegeneracyResult& d, const Options& o,
                              BASIC::CloudPtr& output) {
  SamplingResult result;
  if (!o.enable_informed_sampling || !map || !d.valid || d.status == Status::NORMAL ||
      !R.allFinite() || !t.allFinite()) return result;
  BASIC::CloudPtr finite(new BASIC::PointCloudType()), fine(new BASIC::PointCloudType());
  finite->header=input->header;
  for (const auto& p : *input) if (p.getVector3fMap().allFinite()) finite->push_back(p);
  if (finite->empty()) return result;
  VoxelGridClosest<BASIC::PointType> filter;
  filter.setLeafSize(o.fine_resolution);
  filter.setInputCloud(finite);
  filter.filter(fine);
  auto worldPoint = [&](const BASIC::PointType& p) -> V3 { return R*V3(p.x,p.y,p.z)+t; };
  // 仅对已有可靠高度图的区域排名；提高几何丰富区域的点密度，不生成虚拟点。
  std::map<Key,double> regions;
  bool map_available = false;
  for (const auto& p : *fine) {
    const V3 world=worldPoint(p);
    const auto* layer=map->layer(world);
    if (!layer || !layer->valid) continue;
    map_available = true;
    const double mid=map->mid(world);
    if (mid >= o.mid_threshold && mid > 0) regions[map->key(world)]=mid;
  }
  if (!map_available) return result;
  std::vector<std::pair<Key,double>> ranked(regions.begin(),regions.end());
  std::sort(ranked.begin(),ranked.end(),[](const auto& a,const auto& b) {
    return a.second != b.second ? a.second > b.second : a.first < b.first;
  });
  std::set<Key> selected;
  for (size_t i=0;i<std::min(ranked.size(),size_t(o.informed_top_k_voxels));++i)
    selected.insert(ranked[i].first);
  BASIC::CloudPtr coarse_input(new BASIC::PointCloudType()), coarse(new BASIC::PointCloudType());
  BASIC::CloudPtr combined(new BASIC::PointCloudType());
  combined->header=coarse_input->header=input->header;
  for (const auto& p : *fine) {
    const V3 world=worldPoint(p);
    // Only query validated map coordinates; layer() rejects nonfinite or oversized positions.
    if (map->layer(world) && selected.count(map->key(world))) combined->push_back(p);
    else coarse_input->push_back(p);
  }
  result.fine_points=combined->size();
  filter.setLeafSize(o.coarse_resolution);
  filter.setInputCloud(coarse_input);
  filter.filter(coarse);
  *combined += *coarse;
  combined->header=input->header;
  combined->is_dense=true;
  // 原始采样仍由调用者持有；只有成功形成新采样才替换输出。
  output.swap(combined);
  result.applied=true;
  result.selected_voxels=selected.size();
  result.coarse_points=coarse->size();
  return result;
}
} }
