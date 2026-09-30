
#ifndef VOXEL_GRID_CLOSEST_H
#define VOXEL_GRID_CLOSEST_H


#include <pcl/point_cloud.h>
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "tsl/robin_hood.h"


namespace LI2Sup {


template<typename PointType>
class VoxelGridClosest {
private:
  using Point = PointType;
  using PointCloud = pcl::PointCloud<Point>;
  using CloudPtr = typename PointCloud::Ptr;

  CloudPtr cloud_;
  float voxel_size_ = 0.5f;
  float inv_voxel_size_ = 2.0f;
  // 保存完整有符号坐标作为键。哈希碰撞由相等比较解决，不能把坐标
  // 直接按 15 bit 拼接：负坐标和超过 32767 的坐标会覆盖相邻轴。
  using Key = std::array<int,3>;
  struct KeyHash {
    size_t operator()(const Key& key) const {
      size_t h = 0;
      for (int x : key) h ^= std::hash<int>{}(x)+0x9e3779b9+(h<<6)+(h>>2);
      return h;
    }
  };
  robin_hood::unordered_flat_map<Key, std::size_t, KeyHash> voxel_map_;

  std::vector<Point, Eigen::aligned_allocator<Point>> points_;
  std::vector<double> dist2_;


public:
  VoxelGridClosest() {
    dist2_.reserve(10000);
    points_.reserve(10000);
    voxel_map_.reserve(10000);
  }

  void setLeafSize(float lx) {
    if (!std::isfinite(lx) || lx <= 0 || !std::isfinite(1.0f/lx))
      throw std::invalid_argument("voxel size must be finite and positive");
    voxel_size_ = lx;
    inv_voxel_size_ = 1.0f / lx;
  }

  void setInputCloud(const CloudPtr& cloud) {
    cloud_ = cloud;
  }

  void filter(CloudPtr& output) {
    voxel_map_.clear();
    dist2_.clear();
    points_.clear();

    if (!output) output.reset(new PointCloud());
    if (!cloud_) { output->clear(); return; }
    for (const auto& pt : cloud_->points) {
      const Eigen::Vector3d pf = pt.getVector3fMap().template cast<double>();
      if (!pf.allFinite()) continue;
      const Eigen::Vector3d cell = (pf * double(inv_voxel_size_)).array().round();
      // 转换为整型前检查范围，避免极端输入触发未定义行为。
      if ((cell.array() < std::numeric_limits<int>::min()).any() ||
          (cell.array() > std::numeric_limits<int>::max()).any()) continue;
      const Key key{int(cell.x()), int(cell.y()), int(cell.z())};
      const Eigen::Vector3d center = double(voxel_size_) * cell;
      const double d2 = (pf-center).squaredNorm();

      auto it = voxel_map_.find(key);
      if (it == voxel_map_.end()) {
        voxel_map_.emplace(key, points_.size());
        points_.push_back(pt);
        dist2_.push_back(d2);
      } else if (d2 < dist2_[it->second]) {
        points_[it->second] = pt;
        dist2_[it->second] = d2;
      }
    }

    output->points.swap(points_);
    output->width = output->points.size();
    output->height = 1;
    output->is_dense = true;
    output->header = cloud_->header;
  }
};

}
#endif // VOXEL_GRID_CLOSEST_H