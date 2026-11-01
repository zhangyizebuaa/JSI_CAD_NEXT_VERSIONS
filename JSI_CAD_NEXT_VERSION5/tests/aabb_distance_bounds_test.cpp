#include "distance_analysis/aabb_distance_bounds.hpp"

#include <cassert>
#include <cmath>

int main() {
  const float c1[] = {0, 0, 0, 10, 0, 0};
  const float r1[] = {1, 1, 1, 1, 1, 1};
  const float c2[] = {2, 0, 0, 20, 0, 0};
  const float r2[] = {1, 1, 1, 1, 1, 1};

  const auto touching = aabb_distance_bounds_sq(c1, r1, 0, c2, r2, 0);
  assert(touching.lower == 0.0);
  assert(touching.upper == 24.0);

  const auto separated = aabb_distance_bounds_sq(c1, r1, 1, c2, r2, 1);
  assert(separated.lower == 64.0);
  assert(separated.upper == 152.0);
  assert(separated.lower > touching.upper);
  assert(!(touching.lower > touching.upper));

  const auto equal_boxes = aabb_distance_bounds_sq(c1, r1, 0, c1, r1, 0);
  assert(equal_boxes.lower == 0.0);
  assert(equal_boxes.upper == 12.0);
}
