#pragma once
#include <Eigen/Core>
#include <array>
#include <map>
#include <memory>
#include <vector>
#include <cstdint>

namespace LI2Sup { namespace geometry {
using V3 = Eigen::Vector3d;
using M3 = Eigen::Matrix3d;
using V6 = Eigen::Matrix<double, 6, 1>;
using M6 = Eigen::Matrix<double, 6, 6>;
using Key = std::array<int, 3>;
struct Options {
#define GOPT(type, name, value) type name = value;
#include "lio/geometry_options.def"
#undef GOPT
  bool valid(double voxel_size) const;
};
enum class Status { NORMAL, WEAK, DEGENERATE };
struct DegeneracyResult {
  V6 eigenvalues = V6::Zero();
  M6 eigenvectors = M6::Identity(); // eigenvectors in scaled pose coordinates
  std::array<bool, 6> weak{};
  double ratio = 0;
  int weak_dim = 0;
  Status status = Status::NORMAL;
  bool valid = false;
};
class Analyzer {
 public:
  DegeneracyResult analyze(const M6& information, const Options& o, bool advance);
 private:
  Status status_ = Status::NORMAL;
};
struct Surface {
  double residual = 0, mid = 0, gradient = 0, confidence = 0;
  V3 normal = V3::Zero(); // d residual / d world point; includes height gradient
};
// Statistics describe the bounded retained history, not an unbounded lifetime count.
struct VoxelGeometryStats {
  V3 sum = V3::Zero();
  M3 sum_outer = M3::Zero();
  uint32_t count = 0;
  V3 centroid = V3::Zero();
  V3 normal = V3::Zero();
  V3 eigenvalues = V3::Zero(); // ascending covariance spectrum; normal is minimum mode
  bool plane_valid = false;
};
struct BumpLayer {
  double resolution = 0;
  int width = 0, height = 0;
  std::vector<double> image, weight;
  M3 R_CG = M3::Identity();
  V3 t_CG = V3::Zero(); // pC = R_CG * pG + t_CG
  double mid = 0;
  bool valid = false;
  bool observedPixel(size_t index) const;
  void updateMid(); // Unweighted mean over finite, observed pixels only.
  bool query(const V3& world_point, Surface& out) const;
};
class BumpMap {
 public:
  BumpMap(const Options& options, double voxel_size) : o_(options), size_(voxel_size) {}
  void insert(const std::vector<V3>& points, const V3& sensor);
  bool query(const V3& point, Surface& out) const;
  double mid(const V3& point) const;
  Key key(const V3& p) const;
  size_t size() const { return cells_.size(); }
  // Read-only inspection. Returned pointers must not be retained across insert().
  const VoxelGeometryStats* stats(const V3& point) const;
  const BumpLayer* layer(const V3& point) const;
 private:
  struct Sample { V3 p; double w; };
  struct Cell {
    std::vector<Sample> samples;
    std::vector<Sample> pending; // New observations since the last image update.
    size_t cursor = 0;
    uint64_t touched = 0;
    V3 previous_normal = V3::Zero();
    VoxelGeometryStats geom;
    std::unique_ptr<BumpLayer> layer;
  };
  void rebuild(Cell& c);
  Options o_;
  double size_;
  uint64_t clock_ = 0;
  std::map<Key, Cell> cells_;
  std::map<uint64_t, Key> ages_;
};
V6 poseJacobian(const V3& body, const M3& rotation, const V3& normal);
struct Candidate {
  size_t index = 0;
  V6 J = V6::Zero();
  double residual = 0, quality = 0, variance = 0, weak_score = 0;
  double mid = 0, gradient = 0;
  double w_mid = 0, w_grad = 0, w_pixel = 0, w_weak = 0;
  double local_quality = 0, robust_weight = 0;
  double degeneracy_severity = 0, variance_factor = 1;
};
// Strict quality threshold followed by deterministic Top-K of unique point indices.
void selectWeakDirectionConstraints(std::vector<Candidate>& candidates,
                                    const DegeneracyResult& degeneracy, const Options& options);
// Replaces only final selected rows. Filters invalid/duplicate candidates in place.
// A rejected candidate leaves its original plane row untouched.
void applyExclusiveBumps(const std::vector<V6>& plane_jacobians,
                         const std::vector<double>& plane_residuals,
                         const std::vector<unsigned char>& plane_valid,
                         std::vector<Candidate>& selected, M6& information, V6& rhs,
                         double plane_precision = 1000.0,
                         const std::vector<double>* plane_precisions = nullptr);
struct Diagnostics {
  // 最终观测系统中降权方向的数量和最小权重，非轨迹精度指标。
  int spectral_attenuated = 0;
  double spectral_min_gain = 1;
  bool spectral_valid = false;
  int candidates = 0, accepted = 0, planes = 0;
  int rejected_map = 0, rejected_mid = 0, rejected_gradient = 0;
  int rejected_pixel = 0, rejected_weak = 0, rejected_residual = 0;
  double avg_mid = 0, avg_gradient = 0, avg_weak = 0, avg_covariance = 0;
  double plane_rms = 0, bump_rms = 0, analysis_ms = 0, bump_ms = 0;
  bool fallback = false;
};
// Attenuate unreliable modes of the FINAL measurement system, including bumps.
// Atomic on failure. Never modifies the IMU prior or adds measurement information.
bool applySpectralReliability(M6& information, V6& rhs, const Options& options,
                              Diagnostics& diagnostics);
// Huber IRLS weight; invalid inputs return zero and must be rejected.
double huberWeight(double residual, double delta);
bool degeneracySeverity(double ratio, const Options& options, double& severity);
// Base formula uses variance_factor=1; optional severity scales the numerator only.
// Returns false for invalid inputs; callers must reject the measurement.
bool adaptiveBumpVariance(double quality, const Options& options, double& variance,
                          double variance_factor = 1.0);
bool candidate(const Surface& s, const V3& body, const M3& rotation,
               const DegeneracyResult& d, const Options& o, Candidate& c, Diagnostics& stats);
} }
