#pragma once
#include "common/ds.h"
#include "lio/motion_uncertainty.h"
#include <algorithm>
#include <vector>
#include <limits>

namespace LI2Sup {
// States must be strictly increasing. Clamp boundary points to the supported
// interval; never dereference end() or divide by a repeated timestamp.
inline BASIC::V3 deskewPoint(const BASIC::V3& point_imu, double time,
                            const std::vector<DynamicState>& states, bool smooth = false) {
  if (!std::isfinite(time) || !point_imu.allFinite())
    return BASIC::V3::Constant(std::numeric_limits<BASIC::scalar>::quiet_NaN());
  if (states.empty()) return point_imu;
  const auto& end = states.back();
  if (states.size() == 1 || time >= end.time) return point_imu;
  BASIC::M3 rotation;
  BASIC::V3 position;
  if (time <= states.front().time) {
    rotation = states.front().R;
    position = states.front().p;
  } else {
    // 找到点时间右侧的状态；边界已提前处理，因此 next 和前驱均可解引用。
    const auto next = std::upper_bound(states.begin(), states.end(), time,
        [](double t, const DynamicState& state) { return t < state.time; });
    const auto& head = *std::prev(next);
    const double tau = time - head.time;
    const double fraction = tau / (next->time - head.time);
    rotation = BASIC::Quat(head.R).slerp(fraction, BASIC::Quat(next->R)).toRotationMatrix();
    position = head.p + head.v * tau + 0.5 * next->a * tau * tau;
    if (smooth) {
      Eigen::Matrix3d Iv, Ip;
      motion::integrals(next->w.cast<double>(), tau, Iv, Ip);
      rotation = head.R * BASIC::SO3::Exp(next->w, tau).R_;
      position = head.p + head.v*tau +
          (head.R.cast<double>()*Ip*next->specific_force.cast<double>()).cast<BASIC::scalar>() +
          0.5*next->gravity*tau*tau;
    }
  }
  // 先从采样时刻 IMU 系变换到世界系，再变回扫描结束时刻 IMU 系。
  return end.R.transpose() * (rotation * point_imu + position - end.p);
}
}
