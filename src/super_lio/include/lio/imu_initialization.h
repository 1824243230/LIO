#pragma once
#include "common/ds.h"
#include <limits>

namespace LI2Sup {
// 同步器会借用扫描结束后的 IMU，并在下一帧再次返回它。
// 统计必须按时间去重，且属于估计器实例，不能使用跨实例的 static 均值。
struct ImuInitialization {
  size_t count = 0;
  double last_time = -std::numeric_limits<double>::infinity();
  BASIC::V3 mean_gyro = BASIC::V3::Zero();
  BASIC::V3 mean_acc = BASIC::V3::Zero();
  void add(const IMUData& imu) {
    if (!std::isfinite(imu.secs) || imu.secs <= last_time ||
        !imu.acc.allFinite() || !imu.gyr.allFinite()) return;
    last_time = imu.secs;
    ++count;
    mean_gyro += (imu.gyr-mean_gyro)/count;
    mean_acc += (imu.acc-mean_acc)/count;
  }
};
}
