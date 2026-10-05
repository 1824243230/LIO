#include "basic/alias.h"
#include <pcl/io/pcd_io.h>
#include "OctVoxMap/OctVoxMap.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
using namespace BASIC;
using Map=LI2Sup::OctVoxMap<V3,scalar>;
int failures=0;
void check(bool ok,const char* message) { if (!ok) { ++failures; std::cerr<<message<<'\n'; } }
int main() {
  Map single(Map::Options{0.5f,1});
  single.insert(Map::Points{V3(0,0,0)});
  std::vector<float> points;single.getMap(points);
  check(points.size()==3,"capacity=1 must retain the first voxel");
  single.insert(Map::Points{V3(2,0,0)});single.getMap(points);
  check(points.size()==3 && points[0]==2,"capacity=1 evicts only when a second voxel arrives");
  Map lru(Map::Options{0.5f,2});
  lru.insert(Map::Points{V3(0,0,0),V3(2,0,0)});lru.getMap(points);
  check(points.size()==6,"exactly capacity voxels must be retained");
  lru.insert(Map::Points{V3(0,0,0),V3(4,0,0)});lru.getMap(points);
  check(points.size()==6 && points[0]==4 && points[3]==0,"recently refreshed voxel survives LRU eviction");

  // 六个点落在不同子体素，查询的真正第五近点位于后续搜索组。
  // 固定米制停止距离在小分辨率下会过早退出。
  for(float resolution : {0.1f,0.25f,0.5f,1.0f}) {
    const float scale=resolution/0.1f;
    Map::Points cloud{V3(.001,.001,.001),V3(.001,.099,.001),V3(.001,.001,.099),
                      V3(.001,.099,.099),V3(.099,.099,.099),V3(.101,.001,.001)};
    const V3 query=scale*V3(.049,.001,.001);
    for(auto& p:cloud) p*=scale;
    Map map(Map::Options{resolution,100});map.insert(cloud);
    Map::KNNHeapType heap;map.getTopK(query,heap);
    std::vector<float> expected,actual;
    for(const auto& p:cloud) expected.push_back((p-query).squaredNorm());
    std::sort(expected.begin(),expected.end());expected.resize(5);
    for(int i=0;i<heap.count;++i) actual.push_back(heap.dist2_[i]);
    std::sort(actual.begin(),actual.end());
    bool equal=actual.size()==expected.size();
    for(size_t i=0;equal && i<actual.size();++i) equal=std::abs(actual[i]-expected[i])<1e-6*scale*scale;
    if(!equal) std::cerr<<"resolution="<<resolution<<'\n';
    check(equal,"grouped KNN must agree with brute force for contained points at each resolution");
  }
  Map defaults;
  const V3 query(2.01f,0.01f,0.01f);
  defaults.insert(Map::Points{query});
  Map::KNNHeapType heap; defaults.getTopK_VN(query,heap);
  check(heap.count==1 && heap.dist2_[0]==0,"default map resolution must agree between insertion and VN query");
  for(float resolution : {0.0f,-1.0f,std::numeric_limits<float>::infinity(),
                           std::numeric_limits<float>::quiet_NaN()}) {
    bool rejected=false;
    try { Map invalid(Map::Options{resolution,10}); } catch(const std::invalid_argument&) { rejected=true; }
    check(rejected,"invalid map resolution must be rejected before indexing");
  }
  bool rejected=false;
  try { Map invalid(Map::Options{0.5f,0}); } catch(const std::invalid_argument&) { rejected=true; }
  check(rejected,"zero map capacity must be rejected");

  // 清图必须同时清理 resetMap 的短暂禁写状态，否则再次装图也会被跳过。
  Map reset;
  reset.resetMap(std::vector<float>{0,0,0});
  reset.resetMap(std::vector<float>{2,0,0});
  reset.getMap(points);
  check(points.size()==3 && points[0]==2,"consecutive resetMap calls must load the replacement map");
  reset.clear();
  reset.insert(Map::Points{V3(4,0,0)});
  reset.getMap(points);
  check(points.size()==3 && points[0]==4,"clear must remove stale insertion cooldown");

  const float inf=std::numeric_limits<float>::infinity();
  const float nan=std::numeric_limits<float>::quiet_NaN();
  const Map::Points invalid_points{V3(nan,0,0),V3(0,inf,0),V3(0,0,-inf),
                                  V3(1e10f,0,0),V3(-1e10f,0,0)};
  Map bounded(Map::Options{0.5f,1});
  bounded.insert(Map::Points{V3(0,0,0)});
  bounded.insert(invalid_points);
  bounded.getMap(points);
  check(points.size()==3 && points[0]==0,"invalid voxel coordinates must not insert or evict valid map points");
  for(const auto& p:invalid_points) {
    Map::KNNHeapType hknn,vn;
    bounded.getTopK(p,hknn);
    bounded.getTopK_VN(p,vn);
    check(hknn.count==0 && vn.count==0,"invalid query coordinates must not produce neighbors");
  }
  // INT_MAX 不能精确表示成 float；检查前必须转 double，不能把上界舍入到 2^31。
  const float positive_limit=536870912.0f; // 2^29，乘默认子体素倒数 4 后刚好越界。
  const float largest_valid=std::nextafter(positive_limit,0.0f);
  const float negative_limit=-positive_limit; // INT_MIN 端点本身仍合法。
  for(float x:{largest_valid,negative_limit}) {
    Map boundary;
    const V3 point(x,0,0);
    boundary.insert(Map::Points{point});
    Map::KNNHeapType hknn,vn;
    boundary.getTopK(point,hknn);
    boundary.getTopK_VN(point,vn);
    check(hknn.count==1 && hknn.dist2_[0]==0 && vn.count==1 && vn.dist2_[0]==0,
          "representable boundary cells must remain searchable without neighbor offset overflow");
  }
  Map boundary;
  for(float x:{positive_limit,std::nextafter(negative_limit,-inf)}) {
    const V3 point(x,0,0);
    boundary.insert(Map::Points{point});
    Map::KNNHeapType hknn,vn;
    boundary.getTopK(point,hknn);
    boundary.getTopK_VN(point,vn);
    check(hknn.count==0 && vn.count==0,"out of range boundary cells must be rejected consistently");
  }
  boundary.getMap(points);
  check(points.empty(),"out of range boundary insertions must leave the map empty");

  Map tiny(Map::Options{1e-30f,1});
  tiny.insert(Map::Points{V3(1,0,0)});
  tiny.getMap(points);
  check(points.empty(),"finite points with oversized scaled coordinates must be rejected");
  return failures ? 1 : 0;
}
