#pragma once
#include "common/ds.h"
#include <algorithm>
#include <cmath>

namespace LI2Sup {
inline bool finalizeLidarScan(LidarData& scan) {
  if (!scan.pc || !std::isfinite(scan.start_time)) return false;
  auto& points = scan.pc->points;
  points.erase(std::remove_if(points.begin(), points.end(), [](const auto& p) {
    return !std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) ||
           !std::isfinite(p.offset_time);
  }), points.end());
  if (points.empty()) return false;
  double first = points.front().offset_time, last = first;
  for (const auto& point : points) {
    first = std::min(first, double(point.offset_time));
    last = std::max(last, double(point.offset_time));
  }
  // Some Velodyne bags stamp the scan at its end: negative point offsets
  // are valid. Rebase while preserving every point's absolute timestamp.
  scan.end_time = scan.start_time + last;
  scan.start_time += first;
  for (auto& point : points) point.offset_time -= first;
  scan.pc->width = points.size();
  scan.pc->height = 1;
  return std::isfinite(scan.end_time);
}
}
