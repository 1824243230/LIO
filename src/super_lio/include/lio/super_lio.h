

#ifndef SUPER_LIO_H_
#define SUPER_LIO_H_

#include <queue>
#include <vector>
#include <iostream>
#include <cassert>
#include <filesystem>

#include <pcl/io/pcd_io.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/voxel_grid.h>

#include "basic/alias.h"
#include "common/ds.h"
#include "common/timer.h"
#include "params.h"
#include "ESKF.h"
#include "lio/geometry.h"
#include "lio/motion_uncertainty.h"
#include "lio/imu_initialization.h"
#include <fstream>
#include "OctVoxMap/OctVoxMap.hpp"
#include "OctVoxMap/VoxelGridFilter.h"
#include "ros/ROSWrapper.h"

namespace LI2Sup{

class SuperLIO{
public:
  SuperLIO(){};
  ~SuperLIO(){};

  void setROSWrapper(const ROSWrapper::Ptr& wrapper){
    data_wrapper_ = wrapper;
  }
  virtual void init();
  void process();
  void saveMap();
  void printTimeRecord();

protected:
  void stateWaitKFInit();
  void stateWaitMapInit();
  void stateProcess();
  virtual bool kf_init();
  virtual bool map_init();
  bool Propagation_Undistort();
  void DownSample();
  void AnalyzeSamplingGeometry();
  void Observe(bool plane_only = false);
  void InitGeometry();
  void UpdateBumpMap();
  void PublishGeometry();
  virtual void UpdateMap();
  virtual void Output();
  void caceData();
  void ProcessCaceMap();
  bool FlushMapFragment();

  using StateFn = void (SuperLIO::*)();
  using OctVoxMapType = OctVoxMap<BASIC::V3, BASIC::scalar>;
  using KNNHeapType = KNNHeap<5, BASIC::V3>;
  StateFn state_fn_;
  ESKF::Ptr kf_;
  OctVoxMapType::Ptr ivox_;
  VoxelGridClosest<BASIC::PointType> voxel_grid_fliter_;
  ROSWrapper::Ptr data_wrapper_;
  MeasureGroup measures_;
  
  ImuInitialization imu_initialization_;
  bool flg_init_ = false;
  bool flg_first_scan_ = true;
  std::vector<DynamicState> propagate_states_;
  bool smooth_motion_ = false;
  bool uncertainty_motion_ = false;
  motion::PointCovariances motion_covariances_;
  BASIC::CloudPtr scan_undistort_full_;
  BASIC::CloudPtr ds_undistort_;
  BASIC::CloudPtr original_sample_; // Original center-based cloud retained for full frame rollback.
  BASIC::CloudPtr point_map_, world_pc_, ds_world_;
  int frame_num_ = 0;
  bool observation_valid_ = true;
  bool sampling_geometry_analyzed_ = false;
  BASIC::SE3 sys_init_pose_;
  BASIC::SE3 last_pose_;

  std::size_t effect_knn_num_ = 0;
  BASIC::VV3 points_world_v3_, points_body_v3_;
  std::vector<unsigned char> effect_mask_, effect_knn_mask_;
  std::vector<int> effect_knn_idxs_;
  std::vector<std::pair<BASIC::M6, BASIC::V6>> H_R_;
  std::vector<std::array<double, 4>> abcd_vec_;
  int pcd_index_ = -1;
  int scans_since_flush_ = 0;
  std::string fragment_directory_;
  std::vector<std::string> fragment_paths_;

  geometry::Options geometry_options_;
  geometry::Analyzer geometry_analyzer_;
  geometry::DegeneracyResult degeneracy_;
  geometry::Diagnostics geometry_stats_;
  std::unique_ptr<geometry::BumpMap> bump_map_;
  std::vector<ros::Publisher> geometry_publishers_;
  std::ofstream geometry_csv_;
  Timer time_record_;
};

} // namespace END.

#endif
