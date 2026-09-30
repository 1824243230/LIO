#include "lio/geometry.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <set>

namespace LI2Sup { namespace geometry {
namespace { double unit(double x) { return std::clamp(x, 0.0, 1.0); } }
// 统一校验所有开关组合、迟滞阈值及方差边界，避免部分模块带着非法配置运行。
bool Options::valid(double size) const {
#define GOPT(type, name, value) if (!std::isfinite(static_cast<double>(name))) return false;
#include "lio/geometry_options.def"
#undef GOPT
  return size > 0 && std::isfinite(size) && rotation_length_scale > 0 &&
    spectral_full_ratio > 0 && spectral_full_ratio <= 1 &&
    (!enable_spectral_reliability || enable_degeneracy) &&
    degeneracy_ratio_enter_degenerate > 0 &&
    degeneracy_ratio_enter_degenerate < degeneracy_ratio_exit_degenerate &&
    degeneracy_ratio_exit_degenerate < degeneracy_ratio_enter_weak &&
    degeneracy_ratio_enter_weak < degeneracy_ratio_exit_weak &&
    degeneracy_ratio_exit_weak <= 1 && degeneracy_eigen_ratio_threshold > 0 &&
    degeneracy_eigen_ratio_threshold <= 1 && bump_capacity > 0 &&
    bump_min_points >= 6 && bump_history_points >= bump_min_points &&
    bump_resolution > 0 && size / bump_resolution <= 64 &&
    bump_min_roughness >= 0 && bump_plane_ratio_max > 0 && bump_plane_ratio_max < 1 &&
    bump_plane_reproject_angle_deg > 0 && bump_plane_reproject_angle_deg < 90 &&
    bump_range_weight_max > 0 && bump_range_epsilon > 0 &&
    fine_resolution > 0 && coarse_resolution >= fine_resolution &&
    informed_top_k_voxels > 0 && max_bump_constraints > 0 && mid_threshold >= 0 &&
    bump_min_mid >= 0 && bump_min_gradient > 0 && bump_min_pixel_confidence > 0 &&
    bump_min_weak_score > 0 && weak_direction_score_threshold >= 0 &&
    weak_direction_score_threshold <= 1 && bump_mid_scale > 0 && bump_gradient_scale > 0 &&
    bump_pixel_scale > 0 && bump_weak_scale > 0 && bump_residual_max > 0 && bump_huber_delta > 0 &&
    sigma_b_min > 0 && sigma_b_base >= sigma_b_min && sigma_b_max >= sigma_b_base &&
    bump_degenerate_variance_factor > 0 && bump_degenerate_variance_factor <= 1 &&
    bump_max_translation_correction > 0 && bump_max_rotation_correction > 0;
}
bool applySpectralReliability(M6& information, V6& rhs, const Options& o,
                              Diagnostics& diagnostics) {
  diagnostics.spectral_valid = false;
  diagnostics.spectral_attenuated = 0;
  diagnostics.spectral_min_gain = 1;
  if (!information.allFinite() || !rhs.allFinite() ||
      !std::isfinite(o.rotation_length_scale) || o.rotation_length_scale <= 0 ||
      !std::isfinite(o.spectral_full_ratio) || o.spectral_full_ratio <= 0 ||
      o.spectral_full_ratio > 1) return false;
  // 用特征长度把旋转和平移放到可比较的尺度：x = S z，z = [ℓδθ, δp]。
  // 信息矩阵作合同变换，右端项也必须同步变换，不能只缩放特征值。
  M6 S = M6::Identity(), inverse_S = M6::Identity();
  S.topLeftCorner<3,3>() /= o.rotation_length_scale;
  inverse_S.topLeftCorner<3,3>() *= o.rotation_length_scale;
  const M6 A = S * (0.5 * information + 0.5 * information.transpose()) * S;
  const V6 b = S * rhs;
  if (!A.allFinite() || !b.allFinite()) return false;
  Eigen::SelfAdjointEigenSolver<M6> eig(A);
  if (eig.info() != Eigen::Success || !eig.eigenvalues().allFinite() ||
      !eig.eigenvectors().allFinite()) return false;
  const double largest = eig.eigenvalues().maxCoeff();
  if (largest <= 0 || eig.eigenvalues().minCoeff() < -1e-8 * largest) return false;
  V6 gains, eigenvalues;
  int attenuated = 0;
  for (int k = 0; k < 6; ++k) {
    const double ratio = std::max(0.0, eig.eigenvalues()[k]) / largest;
    // Treat numerical null modes as null; do not preserve unsupported rhs components.
    const double t = ratio <= 1e-12 ? 0 : unit(ratio / o.spectral_full_ratio);
    gains[k] = t*t*(3-2*t); // C1 smooth transition, exactly one in reliable modes.
    eigenvalues[k] = gains[k] * std::max(0.0, eig.eigenvalues()[k]);
    attenuated += gains[k] < 1;
  }
  // 仅在有弱方向时重建；强约束系统保持原数值，避免无意义的分解重构误差。
  if (attenuated) {
    const M6 next_A = inverse_S * eig.eigenvectors() * eigenvalues.asDiagonal() *
                      eig.eigenvectors().transpose() * inverse_S;
    // 同一 gi 同时作用于二次项和一次项，保证状态更新与后验置信度相匹配。
    const V6 next_b = inverse_S * eig.eigenvectors() * gains.asDiagonal() *
                      eig.eigenvectors().transpose() * b;
    if (!next_A.allFinite() || !next_b.allFinite()) return false;
    // 全部结果验证通过后再提交；失败时调用者仍持有原始观测方程。
    information = 0.5 * next_A + 0.5 * next_A.transpose();
    rhs = next_b;
  }
  diagnostics.spectral_valid = true;
  diagnostics.spectral_attenuated = attenuated;
  diagnostics.spectral_min_gain = gains.minCoeff();
  return true;
}
DegeneracyResult Analyzer::analyze(const M6& information, const Options& o, bool advance) {
  DegeneracyResult d;
  if (!information.allFinite()) return d;
  M6 scale = M6::Identity();
  scale.topLeftCorner<3,3>() /= o.rotation_length_scale;
  const M6 G = scale * (0.5 * (information + information.transpose())) * scale;
  if (!G.allFinite()) return d; // Also reject overflow during symmetrization/scaling.
  Eigen::SelfAdjointEigenSolver<M6> solver(G);
  if (solver.info() != Eigen::Success || !solver.eigenvalues().allFinite() ||
      !solver.eigenvectors().allFinite() || solver.eigenvalues()[5] <= 1e-12 ||
      solver.eigenvalues()[0] < -1e-8 * solver.eigenvalues()[5]) return d;
  const double denominator = solver.eigenvalues()[5] + 1e-12;
  for (int k = 0; k < 6; ++k) {
    d.eigenvalues[k] = std::max(0.0, solver.eigenvalues()[5-k]);
    d.eigenvectors.col(k) = solver.eigenvectors().col(5-k);
    d.weak[k] = d.eigenvalues[k] / denominator < o.degeneracy_eigen_ratio_threshold;
    d.weak_dim += d.weak[k];
  }
  d.ratio = d.eigenvalues[5] / denominator;
  Status next = status_;
  if (d.ratio < o.degeneracy_ratio_enter_degenerate) next = Status::DEGENERATE;
  else if (status_ == Status::DEGENERATE && d.ratio < o.degeneracy_ratio_exit_degenerate) next = Status::DEGENERATE;
  else if (d.ratio < o.degeneracy_ratio_enter_weak) next = Status::WEAK;
  else if (status_ != Status::NORMAL && d.ratio < o.degeneracy_ratio_exit_weak) next = Status::WEAK;
  else next = Status::NORMAL;
  // 每帧只推进一次迟滞状态；同帧迭代和回退重算只更新谱，不能重复推进。
  if (advance) status_ = next;
  d.status = status_;
  d.valid = true;
  return d;
}
bool BumpLayer::observedPixel(size_t index) const {
  return index < image.size() && index < weight.size() &&
    std::isfinite(image[index]) && std::isfinite(weight[index]) && weight[index] > 0;
}
void BumpLayer::updateMid() {
  mid = 0; valid = false;
  if (width <= 0 || height <= 0 || image.size() != size_t(width)*size_t(height) ||
      image.size() != weight.size()) return;
  size_t observed = 0;
  for (size_t i=0; i<image.size(); ++i) if (observedPixel(i)) {
    ++observed;
    mid += (std::abs(image[i])-mid)/observed;
  }
  valid = observed > 0 && std::isfinite(mid);
}
Key BumpMap::key(const V3& p) const {
  return {int(std::floor(p.x()/size_)), int(std::floor(p.y()/size_)), int(std::floor(p.z()/size_))};
}
void BumpMap::insert(const std::vector<V3>& points, const V3& sensor) {
  if (!sensor.allFinite()) return;
  std::set<Key> touched;
  for (const auto& p : points) {
    if (!p.allFinite() || p.cwiseAbs().maxCoeff()/size_ > 1e8) continue;
    const Key k = key(p);
    auto it = cells_.find(k);
    if (it == cells_.end()) {
      if (cells_.size() >= size_t(o_.bump_capacity)) {
        const Key oldest = ages_.begin()->second;
        touched.erase(oldest);
        cells_.erase(oldest);
        ages_.erase(ages_.begin());
      }
      it = cells_.try_emplace(k).first;
    }
    Cell& c = it->second;
    ages_.erase(c.touched);
    c.touched = ++clock_;
    ages_.emplace(c.touched, k);
    Sample s{p, std::min(o_.bump_range_weight_max, 1.0 / ((p-sensor).norm()+o_.bump_range_epsilon))};
    if (c.samples.size() < size_t(o_.bump_history_points)) c.samples.push_back(s);
    else { c.samples[c.cursor] = s; c.cursor = (c.cursor+1) % c.samples.size(); }
    if (c.layer) c.pending.push_back(s);
    touched.insert(k);
  }
  for (const auto& k : touched) {
    auto& c = cells_.at(k);
    rebuild(c);
    std::vector<Sample>().swap(c.pending);
  }
}
void BumpMap::rebuild(Cell& c) {
  c.geom = VoxelGeometryStats{};
  auto& geom = c.geom;
  geom.count = static_cast<uint32_t>(c.samples.size());
  for (const auto& s : c.samples) {
    geom.sum += s.p;
    geom.sum_outer += s.p*s.p.transpose();
  }
  if (geom.count) geom.centroid = geom.sum / geom.count;
  if (c.samples.size() < size_t(o_.bump_min_points)) {
    c.layer.reset(); c.previous_normal.setZero(); return;
  }
  // Centered second pass avoids cancellation of E[pp^T] - E[p]E[p]^T far from origin.
  M3 cov = M3::Zero();
  for (const auto& s : c.samples) { const V3 x = s.p-geom.centroid; cov += x*x.transpose(); }
  cov /= geom.count;
  if (!cov.allFinite()) { c.layer.reset(); c.previous_normal.setZero(); return; }
  Eigen::SelfAdjointEigenSolver<M3> eig(cov);
  if (eig.info() != Eigen::Success || !eig.eigenvalues().allFinite() || !eig.eigenvectors().allFinite()) {
    c.layer.reset(); c.previous_normal.setZero(); return;
  }
  geom.eigenvalues = eig.eigenvalues().cwiseMax(0.0);
  V3 n = eig.eigenvectors().col(0);
  if (n.dot(c.previous_normal) < 0) n = -n;
  geom.normal = n;
  geom.plane_valid = geom.eigenvalues[1] >= 1e-8 &&
    geom.eigenvalues[0]/geom.eigenvalues[1] <= o_.bump_plane_ratio_max;
  // A pure plane remains valid geometry, but does not allocate a height image.
  if (!geom.plane_valid || std::sqrt(geom.eigenvalues[0]) < o_.bump_min_roughness) {
    c.layer.reset(); c.previous_normal.setZero(); return;
  }
  const double cos_limit = std::cos(o_.bump_plane_reproject_angle_deg * std::acos(-1.0)/180);
  const bool stable = c.previous_normal.squaredNorm() > 0.5 && n.dot(c.previous_normal) >= cos_limit;
  c.previous_normal = n;
  if (!c.layer && !stable) return;
  const bool initializing = !c.layer;
  const bool reproject = c.layer && n.dot(c.layer->R_CG.row(2).transpose()) < cos_limit;
  // Weighted incremental mean: algebraically (W*I + w*h)/(W+w), without large W*I.
  auto accumulate = [](BumpLayer& l, const V3& local, double w) {
    if (!local.allFinite() || !std::isfinite(w) || w <= 0) return;
    const double uf = std::round(local.x()/l.resolution+(l.width-1)*0.5);
    const double vf = std::round(local.y()/l.resolution+(l.height-1)*0.5);
    if (uf < 0 || vf < 0 || uf >= l.width || vf >= l.height) return;
    const size_t i = size_t(vf)*l.width+size_t(uf);
    const double total = l.weight[i]+w;
    if (!std::isfinite(total)) return;
    const double height = l.image[i]+(w/total)*(local.z()-l.image[i]);
    if (!std::isfinite(height)) return;
    l.image[i] = height;
    l.weight[i] = total;
  };
  if (initializing || reproject) {
    auto next = std::make_unique<BumpLayer>();
    const V3 x = n.unitOrthogonal();
    next->R_CG.row(0) = x.transpose();
    next->R_CG.row(1) = n.cross(x).transpose();
    next->R_CG.row(2) = n.transpose();
    next->t_CG = -next->R_CG * geom.centroid;
    next->resolution = o_.bump_resolution;
    next->width = int(std::ceil(2 * std::sqrt(3.0) * size_/o_.bump_resolution)) + 3;
    next->height = next->width;
    next->image.assign(next->width*next->height, 0);
    next->weight.assign(next->width*next->height, 0);
    if (reproject) {
      const auto& old = *c.layer;
      for (int v=0; v<old.height; ++v) for (int u=0; u<old.width; ++u) {
        const size_t i = size_t(v)*old.width+u;
        if (!old.observedPixel(i)) continue;
        const V3 local((u-(old.width-1)*0.5)*old.resolution,
                       (v-(old.height-1)*0.5)*old.resolution, old.image[i]);
        const V3 world = old.R_CG.transpose()*(local-old.t_CG);
        // Each old pixel transfers its confidence once; collisions use weighted fusion.
        accumulate(*next, next->R_CG*world+next->t_CG, old.weight[i]);
      }
    }
    c.layer = std::move(next);
  }
  auto& l = *c.layer;
  // Initial allocation seeds retained history once. Existing layers consume only new points,
  // including after reprojection; replaying history here would double-count measurements.
  const auto& updates = initializing ? c.samples : c.pending;
  for (const auto& s : updates) accumulate(l, l.R_CG*s.p+l.t_CG, s.w);
  l.updateMid();
}

const VoxelGeometryStats* BumpMap::stats(const V3& point) const {
  if (!point.allFinite() || point.cwiseAbs().maxCoeff()/size_ > 1e8) return nullptr;
  auto it = cells_.find(key(point));
  return it == cells_.end() ? nullptr : &it->second.geom;
}
const BumpLayer* BumpMap::layer(const V3& point) const {
  if (!point.allFinite() || point.cwiseAbs().maxCoeff()/size_ > 1e8) return nullptr;
  auto it = cells_.find(key(point));
  return it == cells_.end() ? nullptr : it->second.layer.get();
}
double BumpMap::mid(const V3& p) const {
  if (!p.allFinite() || p.cwiseAbs().maxCoeff()/size_ > 1e8) return 0;
  auto it = cells_.find(key(p));
  return it != cells_.end() && it->second.geom.plane_valid && it->second.layer &&
    it->second.layer->valid ? it->second.layer->mid : 0;
}
bool BumpMap::query(const V3& p, Surface& out) const {
  out = Surface{};
  if (!p.allFinite() || p.cwiseAbs().maxCoeff()/size_ > 1e8) return false;
  auto it = cells_.find(key(p));
  if (it == cells_.end() || !it->second.geom.plane_valid || !it->second.layer) return false;
  return it->second.layer->query(p, out);
}
bool BumpLayer::query(const V3& p, Surface& out) const {
  out = Surface{};
  if (!valid || !p.allFinite() || !R_CG.allFinite() || !t_CG.allFinite() ||
      !std::isfinite(resolution) || resolution <= 0 || !std::isfinite(mid) || mid < 0 ||
      width < 2 || height < 2 || image.size() != size_t(width)*size_t(height) ||
      weight.size() != image.size()) return false;
  const V3 local = R_CG*p+t_CG;
  double u = local.x()/resolution+(width-1)*0.5;
  double v = local.y()/resolution+(height-1)*0.5;
  // Check finite coordinates and full bilinear support before any float-to-int conversion.
  if (!local.allFinite() || !std::isfinite(u) || !std::isfinite(v) ||
      u < 0 || v < 0 || u >= width-1 || v >= height-1) return false;
  const int x = int(std::floor(u)), y = int(std::floor(v));
  const size_t base = size_t(y)*width+x;
  const size_t ids[4] = {base, base+1, base+width, base+width+1};
  double confidence = weight[ids[0]];
  for (size_t i : ids) { if (!observedPixel(i)) return false; confidence = std::min(confidence, weight[i]); }
  u -= x; v -= y;
  const double a = image[ids[0]], b = image[ids[1]], c = image[ids[2]], d = image[ids[3]];
  const double h = (1-u)*(1-v)*a+u*(1-v)*b+(1-u)*v*c+u*v*d;
  const double gx = ((1-v)*(b-a)+v*(d-c))/resolution;
  const double gy = ((1-u)*(c-a)+u*(d-b))/resolution;
  Surface result;
  result.residual = local.z()-h;
  result.normal = R_CG.transpose()*V3(-gx,-gy,1);
  result.mid = mid; result.gradient = std::hypot(gx,gy); result.confidence = confidence;
  if (!result.normal.allFinite() || !std::isfinite(result.residual) ||
      !std::isfinite(result.gradient) || !std::isfinite(result.confidence)) return false;
  out = result;
  return true;
}
V6 poseJacobian(const V3& body, const M3& rotation, const V3& normal) {
  V6 J; J.head<3>() = body.cross(rotation.transpose()*normal); J.tail<3>() = normal; return J;
}
void selectWeakDirectionConstraints(std::vector<Candidate>& candidates,
                                    const DegeneracyResult& d, const Options& o) {
  if (!d.valid || d.status == Status::NORMAL || o.max_bump_constraints <= 0 ||
      !std::isfinite(o.weak_direction_score_threshold)) { candidates.clear(); return; }
  candidates.erase(std::remove_if(candidates.begin(), candidates.end(), [&](const Candidate& c) {
    return !std::isfinite(c.quality) || c.quality <= o.weak_direction_score_threshold || c.quality > 1 ||
      !c.J.allFinite() || !std::isfinite(c.residual) || !std::isfinite(c.variance) || c.variance <= 0 ||
      !std::isfinite(c.weak_score) || (o.enable_weak_selection && c.weak_score < o.bump_min_weak_score);
  }), candidates.end());
  std::stable_sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b) {
    return a.quality != b.quality ? a.quality > b.quality : a.index < b.index;
  });
  std::set<size_t> indices;
  size_t count = 0;
  for (const auto& c : candidates) {
    if (!indices.insert(c.index).second) continue;
    candidates[count++] = c;
    if (count == size_t(o.max_bump_constraints)) break;
  }
  candidates.resize(count);
}
// 一个点只能贡献平面或曲面中的一条观测，先减旧行再加新行，防止重复计入信息。
void applyExclusiveBumps(const std::vector<V6>& plane_jacobians,
                         const std::vector<double>& plane_residuals,
                         const std::vector<unsigned char>& plane_valid,
                         std::vector<Candidate>& selected, M6& information, V6& rhs,
                         double plane_precision, const std::vector<double>* plane_precisions) {
  if (!information.allFinite() || !rhs.allFinite() || !std::isfinite(plane_precision) || plane_precision <= 0) {
    selected.clear(); return;
  }
  std::set<size_t> replaced;
  size_t accepted = 0;
  for (const auto& c : selected) {
    const size_t i = c.index;
    if (i >= plane_jacobians.size() || i >= plane_residuals.size() || i >= plane_valid.size() ||
        !plane_valid[i] || replaced.count(i) || !c.J.allFinite() || !std::isfinite(c.residual) ||
        !std::isfinite(c.variance) || c.variance <= 0 || !plane_jacobians[i].allFinite() ||
        !std::isfinite(plane_residuals[i])) continue;
    if (plane_precisions && (i >= plane_precisions->size() ||
        !std::isfinite((*plane_precisions)[i]) || (*plane_precisions)[i] <= 0)) continue;
    const double precision = plane_precisions ? (*plane_precisions)[i] : plane_precision;
    const auto& Jp = plane_jacobians[i];
    const M6 next_information = information - precision*Jp*Jp.transpose() + c.J*c.J.transpose()/c.variance;
    const V6 next_rhs = rhs + precision*Jp*plane_residuals[i] - c.J*c.residual/c.variance;
    if (!next_information.allFinite() || !next_rhs.allFinite()) continue;
    information = next_information;
    rhs = next_rhs;
    replaced.insert(i);
    selected[accepted++] = c;
  }
  selected.resize(accepted);
}
double huberWeight(double residual, double delta) {
  if (!std::isfinite(residual) || !std::isfinite(delta) || delta <= 0) return 0;
  const double magnitude = std::abs(residual);
  return magnitude <= delta ? 1.0 : delta/magnitude;
}
bool degeneracySeverity(double ratio, const Options& o, double& severity) {
  severity = 0;
  const double normal = o.degeneracy_ratio_exit_weak;
  const double degenerate = o.degeneracy_ratio_enter_degenerate;
  if (!std::isfinite(ratio) || ratio < 0 || ratio > 1 || !std::isfinite(normal) ||
      !std::isfinite(degenerate) || degenerate < 0 || normal <= degenerate || normal > 1) return false;
  severity = unit((normal-ratio)/(normal-degenerate+1e-12));
  return true;
}
bool adaptiveBumpVariance(double quality, const Options& o, double& variance, double factor) {
  variance = 0;
  if (!std::isfinite(quality) || quality < 0 || quality > 1 ||
      !std::isfinite(factor) || factor <= 0 || factor > 1 ||
      !std::isfinite(o.sigma_b_base) || !std::isfinite(o.sigma_b_min) || !std::isfinite(o.sigma_b_max) ||
      o.sigma_b_min <= 0 || o.sigma_b_base < o.sigma_b_min || o.sigma_b_max < o.sigma_b_base) return false;
  const double base = o.sigma_b_base*o.sigma_b_base;
  const double minimum = o.sigma_b_min*o.sigma_b_min;
  const double maximum = o.sigma_b_max*o.sigma_b_max;
  if (!std::isfinite(base) || !std::isfinite(maximum) || minimum <= 0) return false;
  variance = std::clamp(base*factor/(quality+1e-12), minimum, maximum);
  return std::isfinite(variance) && variance > 0;
}
bool candidate(const Surface& s, const V3& body, const M3& rotation,
               const DegeneracyResult& d, const Options& o, Candidate& c, Diagnostics& stats) {
  ++stats.candidates;
  const size_t index = c.index;
  c = Candidate{}; c.index = index;
  if (!body.allFinite() || !rotation.allFinite()) return false;
  if (!d.valid || d.status == Status::NORMAL || !s.normal.allFinite() ||
      !std::isfinite(s.mid) || !std::isfinite(s.gradient) || !std::isfinite(s.confidence)) return false;
  if (s.mid < o.bump_min_mid) { ++stats.rejected_mid; return false; }
  if (s.gradient < o.bump_min_gradient) { ++stats.rejected_gradient; return false; }
  if (s.confidence < o.bump_min_pixel_confidence) { ++stats.rejected_pixel; return false; }
  if (!std::isfinite(s.residual) || std::abs(s.residual) >= o.bump_residual_max) { ++stats.rejected_residual; return false; }
  // 曲面法向含高度梯度；归一化会改变残差导数，因此保留其实际幅值。
  c.J = poseJacobian(body, rotation, s.normal);
  if (!c.J.allFinite()) return false;
  V6 scaled = c.J; scaled.head<3>() /= o.rotation_length_scale;
  if (!scaled.allFinite() || !d.eigenvectors.allFinite() || !std::isfinite(d.ratio)) return false;
  c.weak_score = 0;
  for (int k = 0; k < 6; ++k) if (d.weak[k]) c.weak_score += std::pow(scaled.dot(d.eigenvectors.col(k)), 2);
  if (!std::isfinite(c.weak_score)) { ++stats.rejected_weak; return false; }
  if (o.enable_weak_selection && c.weak_score < o.bump_min_weak_score) { ++stats.rejected_weak; return false; }
  c.w_mid = unit(s.mid/o.bump_mid_scale);
  c.w_grad = unit(s.gradient/o.bump_gradient_scale);
  c.w_pixel = unit(s.confidence/o.bump_pixel_scale);
  c.w_weak = o.enable_weak_selection ? unit(c.weak_score/o.bump_weak_scale) : 1.0;
  c.local_quality = unit(c.w_mid*c.w_grad*c.w_pixel*c.w_weak);
  c.robust_weight = huberWeight(s.residual, o.bump_huber_delta);
  if (c.robust_weight <= 0) return false;
  c.quality = unit(c.local_quality*c.robust_weight);
  if (!std::isfinite(c.quality)) return false;
  if (c.quality <= o.weak_direction_score_threshold) { ++stats.rejected_weak; return false; }
  // Only reliable candidates reach global severity modulation; it cannot rescue bad local quality.
  if (!degeneracySeverity(d.ratio, o, c.degeneracy_severity)) return false;
  c.variance_factor = (1-c.degeneracy_severity) + c.degeneracy_severity*o.bump_degenerate_variance_factor;
  if (o.enable_adaptive_covariance) {
    if (!adaptiveBumpVariance(c.quality, o, c.variance, c.variance_factor)) return false;
  } else {
    c.variance = o.sigma_b_base*o.sigma_b_base/c.robust_weight;
  }
  c.residual = s.residual; c.mid = s.mid; c.gradient = s.gradient;
  return std::isfinite(c.variance) && c.variance > 0;
}
} }
