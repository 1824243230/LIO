#pragma once
#include "common/ds.h"
#include <cmath>
#include <deque>
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

// Optional initialization from a recent quiet interval. Running means over all
// startup IMU samples are unsafe when the bag starts while the sensor rotates.
// Use a time window so moving samples age out before estimating the initial bias.
class StationaryImuWindow {
public:
  struct Estimate {
    bool ready = false;
    size_t count = 0;
    double duration = 0.0;
    BASIC::V3 mean_gyro = BASIC::V3::Zero();
    BASIC::V3 mean_acc = BASIC::V3::Zero();
    BASIC::V3 std_gyro = BASIC::V3::Zero();
    BASIC::V3 std_acc = BASIC::V3::Zero();
  };

  // Hawkins has 200 Hz IMU data. One second gives about 200 samples, while
  // min_samples still permits this option on sensors down to roughly 50 Hz.
  static constexpr double kWindowSeconds = 1.0;
  static constexpr double kMinDurationSeconds = 0.95;
  static constexpr size_t kMinSamples = 50;
  static constexpr double kMaxGyroMeanNorm = 0.02; // rad/s
  static constexpr double kMaxGyroStd = 0.01;      // rad/s, each axis
  static constexpr double kMaxAccStd = 0.2;        // m/s^2, each axis

  void add(const IMUData& imu) {
    // The synchronization layer can lend the right-hand bracket twice.
    if (!std::isfinite(imu.secs) || imu.secs <= last_time_ ||
        !imu.acc.allFinite() || !imu.gyr.allFinite()) return;
    last_time_ = imu.secs;
    samples_.push_back(imu);
    while (!samples_.empty() && imu.secs - samples_.front().secs > kWindowSeconds)
      samples_.pop_front();
  }

  Estimate estimate() const {
    Estimate result;
    result.count = samples_.size();
    if (samples_.size() < kMinSamples) return result;
    result.duration = samples_.back().secs - samples_.front().secs;
    if (result.duration < kMinDurationSeconds) return result;

    Eigen::Vector3d mean_gyro = Eigen::Vector3d::Zero();
    Eigen::Vector3d mean_acc = Eigen::Vector3d::Zero();
    for (const auto& sample : samples_) {
      mean_gyro += sample.gyr.cast<double>();
      mean_acc += sample.acc.cast<double>();
    }
    mean_gyro /= samples_.size();
    mean_acc /= samples_.size();
    Eigen::Vector3d var_gyro = Eigen::Vector3d::Zero();
    Eigen::Vector3d var_acc = Eigen::Vector3d::Zero();
    for (const auto& sample : samples_) {
      var_gyro += (sample.gyr.cast<double>() - mean_gyro).cwiseAbs2();
      var_acc += (sample.acc.cast<double>() - mean_acc).cwiseAbs2();
    }
    const Eigen::Vector3d std_gyro = (var_gyro / samples_.size()).cwiseSqrt();
    const Eigen::Vector3d std_acc = (var_acc / samples_.size()).cwiseSqrt();
    result.mean_gyro = mean_gyro.cast<BASIC::scalar>();
    result.mean_acc = mean_acc.cast<BASIC::scalar>();
    result.std_gyro = std_gyro.cast<BASIC::scalar>();
    result.std_acc = std_acc.cast<BASIC::scalar>();
    result.ready = mean_gyro.allFinite() && mean_acc.allFinite() &&
                   std_gyro.allFinite() && std_acc.allFinite() &&
                   mean_gyro.norm() <= kMaxGyroMeanNorm &&
                   std_gyro.maxCoeff() <= kMaxGyroStd &&
                   std_acc.maxCoeff() <= kMaxAccStd &&
                   mean_acc.norm() > 1e-6;
    return result;
  }

private:
  double last_time_ = -std::numeric_limits<double>::infinity();
  std::deque<IMUData> samples_;
};
}
