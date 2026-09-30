#include "lio/super_lio.h"
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <cstdlib>
using namespace LI2Sup;
using namespace BASIC;
void check(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
class Harness : public SuperLIO {
public:
  Harness() {
    kf_=std::make_shared<ESKF>();
    scan_undistort_full_.reset(new PointCloudType());
    ds_undistort_.reset(new PointCloudType());
    point_map_.reset(new PointCloudType());world_pc_.reset(new PointCloudType());
  }
  void feed(float x) {
    scan_undistort_full_->clear();
    PointType p;p.x=x;p.y=0;p.z=1;p.intensity=3;
    scan_undistort_full_->push_back(p);
    *ds_undistort_=*scan_undistort_full_;
    caceData();
  }
  size_t cached() const { return point_map_->size(); }
  int index() const { return pcd_index_; }
  const std::string& directory() const { return fragment_directory_; }
  const std::vector<std::string>& fragments() const { return fragment_paths_; }
};
PointCloudType read(const std::filesystem::path& file) {
  PointCloudType cloud;
  check(pcl::io::loadPCDFile<PointType>(file.string(),cloud)==0,"read saved map");
  return cloud;
}
int main() {
  namespace fs=std::filesystem;
  char temporary[]="/tmp/super_lio_map_io_XXXXXX";
  const char* created=mkdtemp(temporary);check(created,"create isolated test directory");
  const fs::path root=created;
  g_save_map=true;g_if_filter=false;g_map_ds_size=10;g_map_name="map.pcd";
  g_save_map_dir=(root/"memory ; literal").string();g_pcd_save_interval=0;
  Harness memory;memory.feed(0.01f);memory.feed(0.02f);memory.saveMap();
  check(read(fs::path(g_save_map_dir)/g_map_name).size()==2,"if_filter=false must preserve same-voxel points in memory mode");

  g_save_map_dir=(root/"fragments ; literal").string();g_pcd_save_interval=1;
  fs::create_directories(fs::path(g_save_map_dir)/"PCD");
  const fs::path old=fs::path(g_save_map_dir)/"PCD/scans_0.pcd";
  {std::ofstream file(old);file<<"previous run must remain intact";}
  Harness chunks;chunks.feed(1);chunks.feed(2);chunks.saveMap();
  const fs::path final=fs::path(g_save_map_dir)/g_map_name;
  check(fs::exists(old) && read(final).size()==2 && chunks.cached()==0 && chunks.fragments().size()==2,
        "fragment session preserves old data and merges only current run");
  // 另一个实例不能覆盖第一个实例的分片。
  Harness second;second.feed(9);
  check(second.directory()!=chunks.directory(),"independent instances have separate fragment directories");
  fs::rename(chunks.fragments()[0],chunks.fragments()[0]+".held");
  chunks.saveMap();
  const auto preserved=read(final);
  check(preserved.size()==2 && preserved.front().x==1,"missing fragment cannot overwrite final map with partial data");

  g_save_map_dir=(root/"retry").string();
  Harness retry;retry.feed(1);
  const fs::path blocked=fs::path(retry.directory())/"scans_1.pcd";
  fs::create_directory(blocked);
  retry.feed(2); // 临时写入成功但重命名失败：必须保留内存，不能递增分片索引。
  check(retry.cached()==1 && retry.index()==0 && retry.fragments().size()==1,"failed fragment write retains cache and index");
  fs::remove(blocked);
  retry.saveMap();
  check(retry.cached()==0 && retry.index()==1 && read(fs::path(g_save_map_dir)/g_map_name).size()==2,
        "retry saves all points exactly once");

  const fs::path blocker=root/"not-a-directory";
  {std::ofstream file(blocker);file<<"keep";}
  g_save_map_dir=(blocker/"child").string();
  Harness unavailable;unavailable.feed(5);unavailable.saveMap();
  check(unavailable.cached()==1 && unavailable.index()==-1,"directory failure preserves unsaved points");
  fs::remove_all(root); // 仅清理本测试 mkdtemp 创建的目录。
  std::cout<<"map IO success/failure/isolation tests passed\n";
}
