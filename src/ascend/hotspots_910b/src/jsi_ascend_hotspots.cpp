#include "jsi_ascend_hotspots.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace {

inline bool AabbOverlap(int i, int j, const float *center1, const float *radius1, const float *center2,
                        const float *radius2) {
  float x_center_1 = center1[i * 3];
  float y_center_1 = center1[i * 3 + 1];
  float z_center_1 = center1[i * 3 + 2];
  float x_half_size_1 = radius1[i * 3];
  float y_half_size_1 = radius1[i * 3 + 1];
  float z_half_size_1 = radius1[i * 3 + 2];

  float x_center_2 = center2[j * 3];
  float y_center_2 = center2[j * 3 + 1];
  float z_center_2 = center2[j * 3 + 2];
  float x_half_size_2 = radius2[j * 3];
  float y_half_size_2 = radius2[j * 3 + 1];
  float z_half_size_2 = radius2[j * 3 + 2];

  return std::fabs(x_center_1 - x_center_2) < (x_half_size_1 + x_half_size_2) &&
         std::fabs(y_center_1 - y_center_2) < (y_half_size_1 + y_half_size_2) &&
         std::fabs(z_center_1 - z_center_2) < (z_half_size_1 + z_half_size_2);
}

}  // namespace

extern "C" int jsi_ascend_distance_topk_mask(const float *center1, size_t cnt1, const float *center2, size_t cnt2,
                                             int topk, int *mask1, int *mask2, float *kth_distance) {
  if (!center1 || !center2 || !mask1 || !mask2 || !kth_distance || cnt1 == 0 || cnt2 == 0 || topk <= 0) return -1;

  std::memset(mask1, 0, cnt1 * sizeof(int));
  std::memset(mask2, 0, cnt2 * sizeof(int));

  const size_t total = cnt1 * cnt2;
  std::vector<float> distance(total, 0.0f);

  for (long long idx = 0; idx < static_cast<long long>(total); ++idx) {
    size_t i = static_cast<size_t>(idx) / cnt2;
    size_t j = static_cast<size_t>(idx) % cnt2;
    float dx = std::fabs(center1[i * 3] - center2[j * 3]);
    float dy = std::fabs(center1[i * 3 + 1] - center2[j * 3 + 1]);
    float dz = std::fabs(center1[i * 3 + 2] - center2[j * 3 + 2]);
    distance[idx] = std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  size_t k = std::min(static_cast<size_t>(topk), total) - 1;
  std::vector<float> scratch = distance;
  std::nth_element(scratch.begin(), scratch.begin() + static_cast<long long>(k), scratch.end());
  *kth_distance = scratch[k];

  for (long long i = 0; i < static_cast<long long>(cnt1); ++i) {
    int hit = 0;
    const size_t base = static_cast<size_t>(i) * cnt2;
    for (size_t j = 0; j < cnt2; ++j) {
      if (distance[base + j] <= *kth_distance) {
        hit = 1;
        break;
      }
    }
    mask1[i] = hit;
  }

  for (long long j = 0; j < static_cast<long long>(cnt2); ++j) {
    int hit = 0;
    for (size_t i = 0; i < cnt1; ++i) {
      if (distance[i * cnt2 + static_cast<size_t>(j)] <= *kth_distance) {
        hit = 1;
        break;
      }
    }
    mask2[j] = hit;
  }

  return 0;
}

extern "C" int jsi_ascend_intersection_mask(const float *center1, const float *radius1, int cnt1, const float *center2,
                                            const float *radius2, int cnt2, int *mask1, int *mask2) {
  if (!center1 || !radius1 || !center2 || !radius2 || !mask1 || !mask2 || cnt1 <= 0 || cnt2 <= 0) return -1;

  std::memset(mask1, 0, static_cast<size_t>(cnt1) * sizeof(int));
  std::memset(mask2, 0, static_cast<size_t>(cnt2) * sizeof(int));

  for (int i = 0; i < cnt1; ++i) {
    int hit = 0;
    for (int j = 0; j < cnt2; ++j) {
      if (AabbOverlap(i, j, center1, radius1, center2, radius2)) {
        hit = 1;
        break;
      }
    }
    mask1[i] = hit;
  }

  for (int j = 0; j < cnt2; ++j) {
    int hit = 0;
    for (int i = 0; i < cnt1; ++i) {
      if (AabbOverlap(i, j, center1, radius1, center2, radius2)) {
        hit = 1;
        break;
      }
    }
    mask2[j] = hit;
  }

  return 0;
}
