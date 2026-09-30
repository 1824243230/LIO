#include "lio/super_lio.h"
#include <atomic>
#include <chrono>

namespace LI2Sup {
namespace {
namespace fs = std::filesystem;
std::string uniqueSuffix() {
  static std::atomic<unsigned long> sequence{0};
  return std::to_string(std::chrono::system_clock::now().time_since_epoch().count())+
         "-"+std::to_string(sequence++);
}

// 同目录临时文件写完后再改名。写入/重命名失败都不能覆盖已有的完整地图。
bool writeCloud(const fs::path& path, const BASIC::PointCloudType& cloud) {
  const fs::path temporary = path.string()+".tmp-"+uniqueSuffix();
  try {
    if (pcl::io::savePCDFileBinary(temporary.string(),cloud) < 0)
      throw std::runtime_error("PCD writer returned failure");
    std::error_code error;
    fs::rename(temporary,path,error);
    if (error) throw std::runtime_error(error.message());
    return true;
  } catch (const std::exception& error) {
    LOG(ERROR) << "Cannot save map " << path << ": " << error.what();
    std::error_code ignored;
    fs::remove(temporary,ignored); // 仅清理由本次写操作生成的临时文件。
    return false;
  }
}

bool saveFinalMap(const BASIC::CloudPtr& points) {
  if (!points || points->empty()) return false;
  std::error_code error;
  fs::create_directories(g_save_map_dir,error);
  if (error) { LOG(ERROR) << "Cannot create map directory: " << error.message(); return false; }
  BASIC::PointCloudType filtered;
  if (g_if_filter) {
    if (!std::isfinite(g_map_ds_size) || g_map_ds_size <= 0) {
      LOG(ERROR) << "Cannot filter map: invalid map voxel size";
      return false;
    }
    pcl::VoxelGrid<BASIC::PointType> filter;
    filter.setLeafSize(g_map_ds_size,g_map_ds_size,g_map_ds_size);
    filter.setInputCloud(points);
    filter.filter(filtered);
  } else {
    // 非分片与分片路径应遵循相同的 if_filter 语义。
    filtered = *points;
  }
  const fs::path output = fs::path(g_save_map_dir)/g_map_name;
  if (!writeCloud(output,filtered)) return false;
  LOG(INFO) << "Map saved to " << output << ", points: " << filtered.size();
  return true;
}
}

void SuperLIO::caceData() {
  if (!g_save_map || !point_map_) return;
  const auto state = kf_->GetNavState();
  Eigen::Matrix4f transform = Eigen::Matrix4f::Identity();
  transform.block<3,3>(0,0) = state.R.R_;
  transform.block<3,1>(0,3) = state.p;
  pcl::transformPointCloud(g_if_filter ? *ds_undistort_ : *scan_undistort_full_, *world_pc_, transform);
  if (!world_pc_->empty()) {
    *point_map_ += *world_pc_;
    ++scans_since_flush_;
  }
  if (g_pcd_save_interval > 0 && scans_since_flush_ >= g_pcd_save_interval)
    FlushMapFragment();
}

bool SuperLIO::FlushMapFragment() {
  if (!point_map_ || point_map_->empty()) return true;
  namespace fs = std::filesystem;
  // 每个估计器实例只管理自己的分片，不执行 shell 命令，也不删除旧会话。
  if (fragment_directory_.empty()) {
    const fs::path base = fs::path(g_save_map_dir)/"PCD";
    std::error_code error;
    fs::create_directories(base,error);
    if (error) { LOG(ERROR) << "Cannot create fragment directory: " << error.message(); return false; }
    for (int attempt=0;attempt<32 && fragment_directory_.empty();++attempt) {
      const fs::path directory = base/("run-"+uniqueSuffix());
      if (fs::create_directory(directory,error)) fragment_directory_=directory.string();
      else if (error) { LOG(ERROR) << "Cannot create fragment session: " << error.message(); return false; }
    }
    if (fragment_directory_.empty()) return false;
  }
  const std::string file = (fs::path(fragment_directory_)/("scans_"+std::to_string(pcd_index_+1)+".pcd")).string();
  if (!writeCloud(file,*point_map_)) return false;
  // 成功落盘后才提交索引并清空缓存；失败时保留全部点，供后续扫描或退出重试。
  fragment_paths_.push_back(file);
  ++pcd_index_;
  point_map_->clear();
  scans_since_flush_=0;
  return true;
}

void SuperLIO::ProcessCaceMap() {
  BASIC::CloudPtr merged(new BASIC::PointCloudType());
  // 只合并本实例成功写入的分片，旧运行文件不能混进当前地图。
  for (const auto& path : fragment_paths_) {
    BASIC::PointCloudType fragment;
    try {
      if (pcl::io::loadPCDFile<BASIC::PointType>(path,fragment) < 0)
        throw std::runtime_error("PCD reader returned failure");
    } catch (const std::exception& error) {
      LOG(ERROR) << "Cannot load fragment " << path << ": " << error.what()
                 << "; keeping existing final map";
      return; // 缺少任一分片时禁止发布一个貌似完整的最终地图。
    }
    *merged += fragment;
  }
  saveFinalMap(merged);
}

void SuperLIO::saveMap() {
  if (!g_save_map || !point_map_) return;
  if (point_map_->empty() && fragment_paths_.empty()) return;
  if (g_pcd_save_interval > 0) {
    if (FlushMapFragment()) ProcessCaceMap();
  } else {
    saveFinalMap(point_map_);
  }
}
}
