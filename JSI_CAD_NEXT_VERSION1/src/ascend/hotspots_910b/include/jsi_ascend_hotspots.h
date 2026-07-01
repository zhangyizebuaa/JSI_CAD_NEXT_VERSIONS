#ifndef _JSI_ASCEND_HOTSPOTS_H
#define _JSI_ASCEND_HOTSPOTS_H

#include <cstddef>

#ifdef __cplusplus
extern "C" {
#endif

int jsi_ascend_distance_topk_mask(const float *center1, size_t cnt1, const float *center2, size_t cnt2, int topk,
                                  int *mask1, int *mask2, float *kth_distance);

int jsi_ascend_intersection_mask(const float *center1, const float *radius1, int cnt1, const float *center2,
                                 const float *radius2, int cnt2, int *mask1, int *mask2);

#ifdef __cplusplus
}
#endif

#endif
