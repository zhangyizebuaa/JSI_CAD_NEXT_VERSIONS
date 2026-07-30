#ifndef _ASCEND_HOTSPOTS_HPP
#define _ASCEND_HOTSPOTS_HPP

#include <cstddef>

bool ascendHotspotsReady();

int ascendDistanceTopkMask(const float *center1, size_t cnt1, const float *center2, size_t cnt2, int topk, int *mask1,
                           int *mask2, float *kth_distance);

int ascendIntersectionMask(const float *center1, const float *radius1, int cnt1, const float *center2, const float *radius2,
                          int cnt2, int *mask1, int *mask2);

#endif
