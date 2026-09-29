#pragma once
#include "basic/alias.h"
#include "lio/geometry.h"
namespace LI2Sup { namespace geometry {
struct SamplingResult {
  bool applied = false;
  size_t selected_voxels = 0, fine_points = 0, coarse_points = 0;
};
// output already contains the original center-based sample; inactive/failure paths leave it intact.
SamplingResult informedSample(const BASIC::CloudPtr& input, const M3& prior_rotation,
                              const V3& prior_translation, const BumpMap* map,
                              const DegeneracyResult& degeneracy, const Options& options,
                              BASIC::CloudPtr& output);
} }
