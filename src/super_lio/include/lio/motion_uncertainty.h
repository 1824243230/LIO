#pragma once
#include <Eigen/Dense>
#include <array>
#include <unordered_map>
#include <vector>

namespace LI2Sup::motion {
using V3 = Eigen::Vector3d;
using M3 = Eigen::Matrix3d;
using M6 = Eigen::Matrix<double,6,6>;
using M18 = Eigen::Matrix<double,18,18>;
M3 hat(const V3& v);
// Exact integrals of Exp(w*t) for constant body angular rate and specific force.
void integrals(const V3& w, double dt, M3& velocity, M3& position);
struct Knot {
  double time = 0;
  M3 R = M3::Identity();
  V3 p = V3::Zero();
  M18 covariance = M18::Zero();
  M18 transition = M18::Identity(); // from previous knot
  M18 noise = M18::Zero();         // independent process covariance of that step
  M6 relative = M6::Zero();       // [left rotation, additive translation], scan-end frame
};
void relativeCovariances(std::vector<Knot>& knots);
M3 pointCovariance(const V3& deskewed, double time, const std::vector<Knot>& knots);
struct PointHash {
  size_t operator()(const std::array<float,3>& p) const {
    size_t h = 0;
    for (float x : p) h ^= std::hash<float>{}(x) + 0x9e3779b9 + (h<<6) + (h>>2);
    return h;
  }
};
// Sampling copies actual points. Exact XYZ keys preserve association through both
// voxel selection and informed sampling, without repurposing intensity/padding.
using PointCovariances = std::unordered_map<std::array<float,3>, M3, PointHash>;
}
