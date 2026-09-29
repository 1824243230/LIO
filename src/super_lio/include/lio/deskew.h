#pragma once
#include "common/ds.h"
#include <algorithm>
#include <vector>

namespace LI2Sup {
// States must be strictly increasing. Clamp boundary points to the supported
// interval; never dereference end() or divide by a repeated timestamp.
inline BASIC::V3 deskewPoint(const BASIC::V3& point_imu, double time,
                            const std::vector<DynamicState>& states) {
  if (states.empty()) return point_imu;
  const auto& end = states.back();
  if (states.size() == 1 || time >= end.time) return point_imu;
  BASIC::M3 rotation;
  BASIC::V3 position;
  if (time <= states.front().time) {
    rotation = states.front().R;
    position = states.front().p;
  } else {
    const auto next = std::upper_bound(states.begin(), states.end(), time,
        [](double t, const DynamicState& state) { return t < state.time; });
    const auto& head = *std::prev(next);
    const double tau = time - head.time;
    const double fraction = tau / (next->time - head.time);
    rotation = BASIC::Quat(head.R).slerp(fraction, BASIC::Quat(next->R)).toRotationMatrix();
    position = head.p + head.v * tau + 0.5 * next->a * tau * tau;
  }
  return end.R.transpose() * (rotation * point_imu + position - end.p);
}
}
