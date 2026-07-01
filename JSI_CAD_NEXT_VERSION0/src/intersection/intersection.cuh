#ifndef _INTERSECTION_CUH
#define _INTERSECTION_CUH

struct CheckIntersectTask {
    int aabb_cnt1;
    float *d_aabb_layer1_center;
    float *d_aabb_layer1_radius;
    int aabb_cnt2;
    float *d_aabb_layer2_center;
    float *d_aabb_layer2_radius;
};

struct AABBResult {
    int naabb_width;
    float *d_aabb_center;
    float *d_aabb_radius;
    float *d_aabb_diagonal;

    float ustart;
    float vstart;
    float ustop;
    float vstop;
    int nuv;

    int eval_task_side;
};

struct CheckedIntersectionResult {
    int naabbs_1;
    int *d_mask_1;
    int naabbs_2;
    int *d_mask_2;
};

#endif