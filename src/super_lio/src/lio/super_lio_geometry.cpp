#include "lio/super_lio.h"
#include "lio/informed_sampling.h"
#include <std_msgs/Float64.h>
#include <std_msgs/Float64MultiArray.h>
#include <algorithm>
#include <set>

namespace LI2Sup {
namespace {
const std::vector<std::string> metrics = {
  "degeneracy_ratio", "weak_dim", "degeneracy_status", "bump_candidate_count",
  "bump_accepted_count", "bump_avg_mid", "bump_avg_gradient", "bump_avg_weak_score",
  "bump_avg_covariance", "plane_measurement_count", "bump_update_time", "degeneracy_analysis_time",
  "plane_residual_rms", "bump_residual_rms", "bump_fallback", "rejected_by_map", "rejected_by_mid",
  "rejected_by_gradient", "rejected_by_pixel", "rejected_by_weak_score", "rejected_by_residual",
  "candidate_bump_count", "accepted_bump_count", "weak_dimension_count", "average_weak_score", "average_bump_covariance"};
}
void SuperLIO::InitGeometry() {
  ros::NodeHandle nh;
#define GOPT(type, name, value) nh.param("/lio/geometry/" #name, geometry_options_.name, geometry_options_.name);
#include "lio/geometry_options.def"
#undef GOPT
  if (!geometry_options_.valid(g_ivox_resolution)) {
    ROS_ERROR("Invalid /lio/geometry parameters: disabling geometry extensions");
    geometry_options_ = geometry::Options{};
  }
  if (geometry_options_.enable_bump_layer)
    bump_map_ = std::make_unique<geometry::BumpMap>(geometry_options_, g_ivox_resolution);
  if (geometry_options_.enable_degeneracy) {
    for (const auto& name : metrics)
      geometry_publishers_.push_back(nh.advertise<std_msgs::Float64>("/super_lio/"+name, 10));
    geometry_publishers_.push_back(nh.advertise<std_msgs::Float64MultiArray>("/super_lio/weak_direction_vectors", 10));
    std::string csv;
    nh.param<std::string>("/lio/geometry/csv_path", csv, "");
    if (!csv.empty()) {
      geometry_csv_.open(csv);
      if (!geometry_csv_) ROS_ERROR_STREAM("Cannot open geometry CSV: " << csv);
      else {
        geometry_csv_ << "timestamp";
        for (const auto& name : metrics) geometry_csv_ << ',' << name;
        for (int k=0;k<6;++k) geometry_csv_ << ",lambda" << k;
        for (int k=0;k<6;++k) {
          geometry_csv_ << ",weak_" << k;
          for (int j=0;j<6;++j) geometry_csv_ << ",eigenvector_" << k << "_" << j;
        }
        geometry_csv_ << '\n';
      }
    }
  }
}
void SuperLIO::UpdateBumpMap() {
  if (!bump_map_) return;
  const auto pose = kf_->GetSE3();
  std::vector<geometry::V3> points;
  points.reserve(scan_undistort_full_->size());
  // Use undistorted full-resolution history so ordinary downsampling cannot erase microgeometry.
  for (const auto& p : *scan_undistort_full_) {
    BASIC::V3 body(p.x,p.y,p.z);
    if (body.allFinite()) points.push_back((pose*body).cast<double>());
  }
  bump_map_->insert(points, (pose * g_lidar_imu.t_).cast<double>());
}
void SuperLIO::DownSample() {
  voxel_grid_fliter_.setInputCloud(scan_undistort_full_);
  voxel_grid_fliter_.filter(ds_undistort_);
  original_sample_ = ds_undistort_;
  sampling_geometry_analyzed_ = false;
  const auto pose = kf_->GetSE3();
  // Fine sampling selects original points. If none of the full scan can query
  // a supported historical surface at the prior pose, retain the baseline
  // sample and avoid an extra KNN pass, fine/coarse resampling and replay.
  // Observe still analyzes geometry and attempts measurements at each iterate.
  // Sampling-only ablations deliberately retain their existing behavior.
  if (geometry_options_.enable_informed_sampling && geometry_options_.enable_bump_measurement) {
    bool supported_surface = false;
    if (bump_map_) for (const auto& p : *scan_undistort_full_) {
      const BASIC::V3 body(p.x,p.y,p.z);
      if (!body.allFinite()) continue;
      geometry::Surface surface;
      if (bump_map_->query((pose*body).cast<double>(), surface)) {
        supported_surface = true;
        break;
      }
    }
    if (!supported_surface) return;
  }
  AnalyzeSamplingGeometry();
  geometry::informedSample(scan_undistort_full_, pose.R_.cast<double>(), pose.t_.cast<double>(),
                           bump_map_.get(), degeneracy_, geometry_options_, ds_undistort_);
}
void SuperLIO::PublishGeometry() {
  if (geometry_publishers_.empty()) return;
  const auto& s = geometry_stats_;
  const std::vector<double> values = {
    degeneracy_.ratio, double(degeneracy_.weak_dim), degeneracy_.valid ? double(int(degeneracy_.status)) : -1.0,
    double(s.candidates), double(s.accepted), s.avg_mid, s.avg_gradient, s.avg_weak, s.avg_covariance,
    double(s.planes), s.bump_ms, s.analysis_ms, s.plane_rms, s.bump_rms, double(s.fallback),
    double(s.rejected_map), double(s.rejected_mid), double(s.rejected_gradient), double(s.rejected_pixel),
    double(s.rejected_weak), double(s.rejected_residual), double(s.candidates), double(s.accepted),
    double(degeneracy_.weak_dim), s.avg_weak, s.avg_covariance};
  for (size_t i=0;i<values.size();++i) { std_msgs::Float64 m; m.data=values[i]; geometry_publishers_[i].publish(m); }
  std_msgs::Float64MultiArray vectors;
  vectors.layout.dim.resize(2);
  vectors.layout.dim[0].label="weak_direction";
  vectors.layout.dim[0].size=degeneracy_.valid ? degeneracy_.weak_dim : 0;
  vectors.layout.dim[0].stride=vectors.layout.dim[0].size*6;
  vectors.layout.dim[1].label="rotation_body_xyz_translation_world_xyz_scaled";
  vectors.layout.dim[1].size=6;
  vectors.layout.dim[1].stride=6;
  if (degeneracy_.valid) for (int k=0;k<6;++k) if (degeneracy_.weak[k])
    for (int j=0;j<6;++j) vectors.data.push_back(degeneracy_.eigenvectors(j,k));
  geometry_publishers_.back().publish(vectors);
  if (geometry_csv_) {
    geometry_csv_.precision(16);
    geometry_csv_ << measures_.lidar.end_time;
    for (double v : values) geometry_csv_ << ',' << v;
    for (int k=0;k<6;++k) geometry_csv_ << ',' << degeneracy_.eigenvalues[k];
    for (int k=0;k<6;++k) {
      geometry_csv_ << ',' << int(degeneracy_.valid && degeneracy_.weak[k]);
      for (int j=0;j<6;++j) geometry_csv_ << ',' << degeneracy_.eigenvectors(j,k);
    }
    geometry_csv_ << '\n';
  }
}
}
