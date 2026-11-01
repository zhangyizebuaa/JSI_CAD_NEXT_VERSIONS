#ifndef JSI_CAD_AABB_DISTANCE_BOUNDS_HPP
#define JSI_CAD_AABB_DISTANCE_BOUNDS_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>

struct AabbDistanceBoundsSq {
  double lower;
  double upper;
};

inline AabbDistanceBoundsSq aabb_distance_bounds_sq(const float *center1, const float *radius1, size_t i,
                                                     const float *center2, const float *radius2, size_t j) {
  AabbDistanceBoundsSq bounds{0.0, 0.0};
  for (int axis = 0; axis < 3; ++axis) {
    const size_t a = i * 3 + static_cast<size_t>(axis);
    const size_t b = j * 3 + static_cast<size_t>(axis);
    const double delta = std::fabs(static_cast<double>(center1[a]) - center2[b]);
    const double radii = static_cast<double>(radius1[a]) + radius2[b];
    const double gap = std::max(0.0, delta - radii);
    bounds.lower += gap * gap;
    const double span = delta + radii;
    bounds.upper += span * span;
  }
  return bounds;
}

#endif
