#ifndef _DISTANCE_ANALYSIS_CUH
#define _DISTANCE_ANALYSIS_CUH

#include "common/common.cuh"

struct CheckDistanceTask {
    size_t aabb_cnt1;
    float *d_aabb_layer1_center;
    float *d_aabb_layer1_radius;
    float *d_aabb_layer1_diagonal;
    size_t aabb_cnt2;
    float *d_aabb_layer2_center;
    float *d_aabb_layer2_radius;
    float *d_aabb_layer2_diagonal;
};

struct AABBResult {
    int naabb_width;
    float *d_aabb_center;
    float *d_aabb_radius;
    float *d_aabb_diagonal;
    float *d_aabb_k;      // per-patch curvature discount (0 if CAD_DIST_USE_K=0)

    float ustart;
    float vstart;
    float ustop;
    float vstop;
    int nuv;

    int eval_task_side;
};

struct CheckDistanceResult {
    size_t naabbs_1;
    int *d_mask_1;
    size_t naabbs_2;
    int *d_mask_2;
    float h_kth_distance;
};

#endif