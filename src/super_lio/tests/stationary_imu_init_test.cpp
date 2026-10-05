#include "lio/imu_initialization.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace LI2Sup;
using namespace BASIC;

namespace {
void check(bool ok, const char* message) {
  if (!ok) throw std::runtime_error(message);
}

IMUData sample(double time, const V3& gyro, const V3& acc) {
  IMUData imu;
  imu.secs = time;
  imu.gyr = gyro;
  imu.acc = acc;
  return imu;
}
}

int main() {
  StationaryImuWindow window;
  const V3 gravity(0.53f, 0.23f, 9.80f);
  // A low-variance constant spin is still motion and must not become gyro bias.
  for (int i = 0; i <= 200; ++i) {
    auto imu = sample(i * 0.005, V3(0.0f, 0.0f, -0.75f), gravity);
    window.add(imu);
    window.add(imu); // borrowed right-hand IMU sample from the next scan
  }
  auto estimate = window.estimate();
  check(!estimate.ready && estimate.count == 201, "spin and duplicate samples rejected");

  // The moving samples must leave the window before a quiet mean is accepted.
  for (int i = 201; i <= 399; ++i) {
    const double time = i * 0.005;
    const float wave = std::sin(float(i) * 0.2f);
    window.add(sample(time, V3(0.001f * wave, -0.002f + 0.003f * wave,
                                0.001f * wave),
                      gravity + V3(0.05f * wave, 0.01f * wave, 0.005f * wave)));
  }
  check(!window.estimate().ready, "mixed spin and quiet samples rejected");
  window.add(sample(2.005, V3(0.0f, -0.002f, 0.0f), gravity));
  estimate = window.estimate();
  check(estimate.ready && estimate.count >= 199 && estimate.duration >= 0.95,
        "one-second quiet window accepted");
  check(std::abs(estimate.mean_gyro.z()) < 0.001f &&
        std::abs(estimate.mean_gyro.y() + 0.002f) < 0.001f &&
        (estimate.mean_acc - gravity).norm() < 0.01f,
        "quiet-window means exclude startup motion");

  StationaryImuWindow vibration;
  for (int i = 0; i <= 200; ++i) {
    const float wave = i % 2 ? 1.0f : -1.0f;
    vibration.add(sample(i * 0.005, V3::Zero(),
                         gravity + V3(0.3f * wave, 0.0f, 0.0f)));
  }
  check(!vibration.estimate().ready, "high accelerometer variation rejected");

  StationaryImuWindow short_window;
  for (int i = 0; i < 50; ++i)
    short_window.add(sample(i * 0.005, V3::Zero(), gravity));
  check(!short_window.estimate().ready, "a short burst cannot pass the duration gate");

  std::cout << "stationary IMU initialization passed\n";
}
