#pragma once
#include "basic/alias.h"
#include <Eigen/Eigenvalues>
#include <array>

namespace LI2Sup {
// Centered orthogonal fit also represents planes through the world origin.
inline bool calc_plane_coeff(int count, const std::array<BASIC::V3, 5>& points,
                             std::array<double, 4>& abcd) {
  if (count < 4 || count > 5) return false;
  Eigen::Vector3d center = Eigen::Vector3d::Zero();
  for (int i = 0; i < count; ++i) {
    if (!points[i].allFinite()) return false;
    center += points[i].cast<double>();
  }
  center /= count;
  Eigen::Matrix3d scatter = Eigen::Matrix3d::Zero();
  for (int i = 0; i < count; ++i) {
    const Eigen::Vector3d delta = points[i].cast<double>() - center;
    scatter.noalias() += delta * delta.transpose();
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(scatter);
  if (solver.info() != Eigen::Success || !solver.eigenvalues().allFinite()) return false;
  const auto values = solver.eigenvalues();
  // A line or a point has no identifiable plane normal.
  if (values[1] <= std::max(1e-10, 1e-4 * values[2])) return false;
  const Eigen::Vector3d normal = solver.eigenvectors().col(0);
  const double offset = -normal.dot(center);
  for (int i = 0; i < count; ++i)
    if (std::abs(normal.dot(points[i].cast<double>()) + offset) > 0.1) return false;
  abcd = {normal.x(), normal.y(), normal.z(), offset};
  return true;
}
}
