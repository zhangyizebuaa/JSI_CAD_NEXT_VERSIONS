#include "ascend_hotspots.hpp"

#include <mutex>

#if defined(USE_ASCEND_ACL) && defined(__linux__)
#include <dlfcn.h>

using DistanceFn = int (*)(const float *, size_t, const float *, size_t, int, int *, int *, float *);
using IntersectionFn = int (*)(const float *, const float *, int, const float *, const float *, int, int *, int *);

static void *g_lib = nullptr;
static DistanceFn g_distance_fn = nullptr;
static IntersectionFn g_intersection_fn = nullptr;
static std::once_flag g_once;

static void LoadHotspotLib() {
  g_lib = dlopen("libjsi_ascend_hotspots.so", RTLD_LAZY | RTLD_LOCAL);
  if (!g_lib) return;
  g_distance_fn = reinterpret_cast<DistanceFn>(dlsym(g_lib, "jsi_ascend_distance_topk_mask"));
  g_intersection_fn = reinterpret_cast<IntersectionFn>(dlsym(g_lib, "jsi_ascend_intersection_mask"));
}

#endif

bool ascendHotspotsReady() {
#if defined(USE_ASCEND_ACL) && defined(__linux__)
  std::call_once(g_once, LoadHotspotLib);
  return g_distance_fn != nullptr && g_intersection_fn != nullptr;
#else
  return false;
#endif
}

int ascendDistanceTopkMask(const float *center1, size_t cnt1, const float *center2, size_t cnt2, int topk, int *mask1,
                           int *mask2, float *kth_distance) {
#if defined(USE_ASCEND_ACL) && defined(__linux__)
  std::call_once(g_once, LoadHotspotLib);
  if (!g_distance_fn) return -1;
  return g_distance_fn(center1, cnt1, center2, cnt2, topk, mask1, mask2, kth_distance);
#else
  (void)center1;
  (void)cnt1;
  (void)center2;
  (void)cnt2;
  (void)topk;
  (void)mask1;
  (void)mask2;
  (void)kth_distance;
  return -1;
#endif
}

int ascendIntersectionMask(const float *center1, const float *radius1, int cnt1, const float *center2, const float *radius2,
                           int cnt2, int *mask1, int *mask2) {
#if defined(USE_ASCEND_ACL) && defined(__linux__)
  std::call_once(g_once, LoadHotspotLib);
  if (!g_intersection_fn) return -1;
  return g_intersection_fn(center1, radius1, cnt1, center2, radius2, cnt2, mask1, mask2);
#else
  (void)center1;
  (void)radius1;
  (void)cnt1;
  (void)center2;
  (void)radius2;
  (void)cnt2;
  (void)mask1;
  (void)mask2;
  return -1;
#endif
}
