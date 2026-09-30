#pragma once
#include "common/ds.h"
#include "lio/scan_validation.h"

namespace LI2Sup {
inline bool appendImu(const IMUData& sample, std::deque<IMUData>& buffer, double& last_time) {
  // 迟到/重复包不应清空已经按时收到的队列，更不能把旧 IMU 与新地图混合。
  // 数据集时间整体回绕需要重启整个估计器，仅清队列不足以重置地图与状态。
  if (!std::isfinite(sample.secs) || sample.secs <= last_time ||
      !sample.acc.allFinite() || !sample.gyr.allFinite()) return false;
  buffer.push_back(sample);
  last_time = sample.secs;
  return true;
}

inline bool synchronizeMeasurements(std::deque<LidarData>& lidar, std::deque<IMUData>& imu,
                                    double& last_scan_end, MeasureGroup& output) {
  if (lidar.empty() || imu.empty()) return false;
  const auto& scan = lidar.front();
  if (!std::isfinite(scan.start_time) || !std::isfinite(scan.end_time) ||
      scan.start_time > scan.end_time || scan.end_time <= last_scan_end) {
    lidar.pop_front();
    return false;
  }
  // 未收到右端测量时不能消耗左侧样本，也不把半个扫描写入 output。
  if (imu.back().secs < scan.end_time) return false;
  output.lidar = scan;
  output.imu.clear();
  while (!imu.empty() && imu.front().secs <= scan.end_time) {
    output.imu.push_back(imu.front());
    imu.pop_front();
  }
  // 右侧 IMU 只借用不弹出：同一个测量仍要参与下一帧的剩余区间积分。
  if (!imu.empty() && (output.imu.empty() || output.imu.back().secs < scan.end_time))
    output.imu.push_back(imu.front());
  last_scan_end = scan.end_time;
  lidar.pop_front();
  return true;
}

// 与实际 Livox 回调共用转换函数，防止只测试备用转换器而漏检生产路径。
template<class Message>
bool convertLivoxScan(const Message& msg, int step, double blind2, double maxrange2, LidarData& scan) {
  if (step <= 0 || msg.point_num < 10 || msg.point_num != msg.points.size()) return false;
  scan.pc.reset(new pcl::PointCloud<PointXTZIT>());
  scan.pc->reserve(msg.points.size()/step+1);
  scan.start_time = msg.header.stamp.toSec();
  for (size_t i=0;i<msg.points.size();i+=step) {
    const auto& p = msg.points[i];
    const auto tag = p.tag & 0x30;
    if (tag != 0x10 && tag != 0) continue;
    const double d2 = double(p.x)*p.x+double(p.y)*p.y+double(p.z)*p.z;
    if (!std::isfinite(d2) || d2 <= blind2 || d2 >= maxrange2) continue;
    scan.pc->emplace_back(p.x,p.y,p.z,p.reflectivity,p.offset_time*1e-9);
  }
  return finalizeLidarScan(scan);
}
}
