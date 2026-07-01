#include "intersection.cuh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <utility>

#include <vector>

#include "common/common.cuh"
#include "common/profiler.hpp"
#include "evaluation/evaluation.cuh"
#include "jsi_cad.hpp"

static int global_nuv = 16;

static inline bool AabbOverlap(int i, int j, const float *center1, const float *radius1, const float *center2,
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

static void gather_idx_kernel_host(int *data_out, size_t size, const int *pos, const int *pre) {
  for (size_t i = 0; i < size; ++i) {
    if (pre[i]) data_out[pos[i] - 1] = static_cast<int>(i);
  }
}

static void construct_aabbs_host(int un, int vn, const float *points, float *center, float *radius, float *diagonal) {
  for (int i = 0; i < un - 1; ++i) {
    for (int j = 0; j < vn - 1; ++j) {
      int id = i * (vn - 1) + j;
      int p0 = i * vn + j;
      int p1 = i * vn + j + 1;
      int p2 = (i + 1) * vn + j;
      int p3 = (i + 1) * vn + j + 1;
      float xmin =
          std::min(points[p0 * ONE_UV_SIZE], std::min(points[p1 * ONE_UV_SIZE],
                                                       std::min(points[p2 * ONE_UV_SIZE], points[p3 * ONE_UV_SIZE])));
      float ymin = std::min(points[p0 * ONE_UV_SIZE + 1],
                            std::min(points[p1 * ONE_UV_SIZE + 1],
                                     std::min(points[p2 * ONE_UV_SIZE + 1], points[p3 * ONE_UV_SIZE + 1])));
      float zmin = std::min(points[p0 * ONE_UV_SIZE + 2],
                            std::min(points[p1 * ONE_UV_SIZE + 2],
                                     std::min(points[p2 * ONE_UV_SIZE + 2], points[p3 * ONE_UV_SIZE + 2])));
      float xmax =
          std::max(points[p0 * ONE_UV_SIZE], std::max(points[p1 * ONE_UV_SIZE],
                                                       std::max(points[p2 * ONE_UV_SIZE], points[p3 * ONE_UV_SIZE])));
      float ymax = std::max(points[p0 * ONE_UV_SIZE + 1],
                            std::max(points[p1 * ONE_UV_SIZE + 1],
                                     std::max(points[p2 * ONE_UV_SIZE + 1], points[p3 * ONE_UV_SIZE + 1])));
      float zmax = std::max(points[p0 * ONE_UV_SIZE + 2],
                            std::max(points[p1 * ONE_UV_SIZE + 2],
                                     std::max(points[p2 * ONE_UV_SIZE + 2], points[p3 * ONE_UV_SIZE + 2])));
      float x_span = xmax - xmin;
      float y_span = ymax - ymin;
      float z_span = zmax - zmin;
      center[id * 3] = (xmin + xmax) / 2.0f;
      center[id * 3 + 1] = (ymin + ymax) / 2.0f;
      center[id * 3 + 2] = (zmin + zmax) / 2.0f;
      radius[id * 3] = x_span / 2.0f;
      radius[id * 3 + 1] = y_span / 2.0f;
      radius[id * 3 + 2] = z_span / 2.0f;
      diagonal[id] = std::sqrt(x_span * x_span + y_span * y_span + z_span * z_span);
    }
  }
}

static void para_check_aabbs_host(int aabb_cnt1, const float *center1, const float *radius1, int aabb_cnt2,
                                  const float *center2, const float *radius2, int *mask1, int *mask2) {
  for (int i = 0; i < aabb_cnt1; ++i) {
    int hit = 0;
    for (int j = 0; j < aabb_cnt2; ++j) {
      if (AabbOverlap(i, j, center1, radius1, center2, radius2)) {
        hit = 1;
        break;
      }
    }
    mask1[i] = hit;
  }

  for (int j = 0; j < aabb_cnt2; ++j) {
    int hit = 0;
    for (int i = 0; i < aabb_cnt1; ++i) {
      if (AabbOverlap(i, j, center1, radius1, center2, radius2)) {
        hit = 1;
        break;
      }
    }
    mask2[j] = hit;
  }
}

static AABBResult run_evaluation_and_construction_task(EvalAndConstructTask &task) {
  CAD_PROFILE_SCOPE("run_evaluation_and_construction_task_intersection");
  int nuv = task.nuv;
  float ustep = (task.ustop - task.ustart) / nuv;
  float vstep = (task.vstop - task.vstart) / nuv;
  int naabb_width = nuv - 1;
  float *results = static_cast<float *>(allocate_from_workspace(ONE_UV_SIZE * nuv * nuv * sizeof(float)));
  float *center = static_cast<float *>(allocate_from_workspace(3 * naabb_width * naabb_width * sizeof(float)));
  float *radius = static_cast<float *>(allocate_from_workspace(3 * naabb_width * naabb_width * sizeof(float)));
  float *diagonal = static_cast<float *>(allocate_from_workspace(naabb_width * naabb_width * sizeof(float)));
  EG_spline2dDeriv_irrational(task.d_ivec, task.d_data, ustep, vstep, task.ustart, task.vstart, nuv, results);
  construct_aabbs_host(nuv, nuv, results, center, radius, diagonal);
  AABBResult res;
  res.naabb_width = naabb_width;
  res.d_aabb_center = center;
  res.d_aabb_radius = radius;
  res.d_aabb_diagonal = diagonal;
  res.eval_task_side = task.eval_task_side;
  res.ustart = task.ustart;
  res.vstart = task.vstart;
  res.ustop = task.ustop;
  res.vstop = task.vstop;
  res.nuv = task.nuv;
  return res;
}

static CheckedIntersectionResult check_intersection_of_aabbs(CheckIntersectTask t) {
  CAD_PROFILE_SCOPE("check_intersection_of_aabbs");
  int *mask1 = static_cast<int *>(allocate_from_workspace(t.aabb_cnt1 * sizeof(int)));
  int *mask2 = static_cast<int *>(allocate_from_workspace(t.aabb_cnt2 * sizeof(int)));
  cadMemset(mask1, 0, t.aabb_cnt1 * sizeof(int));
  cadMemset(mask2, 0, t.aabb_cnt2 * sizeof(int));
  // FORMER: always `para_check_aabbs` kernel semantics (no Ascend shortcut).
  para_check_aabbs_host(t.aabb_cnt1, t.d_aabb_layer1_center, t.d_aabb_layer1_radius, t.aabb_cnt2,
                        t.d_aabb_layer2_center, t.d_aabb_layer2_radius, mask1, mask2);
  CheckedIntersectionResult res;
  res.d_mask_1 = mask1;
  res.d_mask_2 = mask2;
  res.naabbs_1 = t.aabb_cnt1;
  res.naabbs_2 = t.aabb_cnt2;
  return res;
}

static std::pair<int, int> find_aabb(int id, std::vector<AABBResult> &aabbs) {
  int total = 0;
  for (int i = 0; i < static_cast<int>(aabbs.size()); i++) {
    auto &aabb = aabbs[i];
    total += aabb.naabb_width * aabb.naabb_width;
    if (total > id) return {i, id - (total - aabb.naabb_width * aabb.naabb_width)};
  }
  return {-1, -1};
}

void intersection(EvalAndConstructTask t1, EvalAndConstructTask t2) {
  CAD_PROFILE_SCOPE("intersection");
  std::queue<EvalAndConstructTask> eval_queue;
  std::vector<AABBResult> aabb_1_queue;
  std::vector<AABBResult> aabb_2_queue;
  eval_queue.push(t1);
  eval_queue.push(t2);
  int round = 0;

  while (true) {
    ++round;
    printf("Round: %d\n", round);
    int naabbs_1 = 0;
    int naabbs_2 = 0;
    while (!eval_queue.empty()) {
      auto task = eval_queue.front();
      eval_queue.pop();
      auto result = run_evaluation_and_construction_task(task);
      if (result.eval_task_side == 1) {
        naabbs_1 += result.naabb_width * result.naabb_width;
        aabb_1_queue.push_back(result);
      } else {
        naabbs_2 += result.naabb_width * result.naabb_width;
        aabb_2_queue.push_back(result);
      }
    }

    float *center_1 = static_cast<float *>(allocate_from_workspace(3 * naabbs_1 * sizeof(float)));
    float *radius_1 = static_cast<float *>(allocate_from_workspace(3 * naabbs_1 * sizeof(float)));
    float *center_2 = static_cast<float *>(allocate_from_workspace(3 * naabbs_2 * sizeof(float)));
    float *radius_2 = static_cast<float *>(allocate_from_workspace(3 * naabbs_2 * sizeof(float)));
    int pos = 0;
    for (auto &aabb : aabb_1_queue) {
      int elem_cnt = aabb.naabb_width * aabb.naabb_width;
      cadMemcpy(center_1 + 3 * pos, aabb.d_aabb_center, elem_cnt * 3 * sizeof(float), CadMemcpyKind::DeviceToDevice);
      cadMemcpy(radius_1 + 3 * pos, aabb.d_aabb_radius, elem_cnt * 3 * sizeof(float), CadMemcpyKind::DeviceToDevice);
      pos += elem_cnt;
    }
    pos = 0;
    for (auto &aabb : aabb_2_queue) {
      int elem_cnt = aabb.naabb_width * aabb.naabb_width;
      cadMemcpy(center_2 + 3 * pos, aabb.d_aabb_center, elem_cnt * 3 * sizeof(float), CadMemcpyKind::DeviceToDevice);
      cadMemcpy(radius_2 + 3 * pos, aabb.d_aabb_radius, elem_cnt * 3 * sizeof(float), CadMemcpyKind::DeviceToDevice);
      pos += elem_cnt;
    }

    CheckIntersectTask check_intersection;
    check_intersection.aabb_cnt1 = naabbs_1;
    check_intersection.aabb_cnt2 = naabbs_2;
    check_intersection.d_aabb_layer1_center = center_1;
    check_intersection.d_aabb_layer2_center = center_2;
    check_intersection.d_aabb_layer1_radius = radius_1;
    check_intersection.d_aabb_layer2_radius = radius_2;
    auto checked_result = check_intersection_of_aabbs(check_intersection);

    int *prefix_sum_1 = static_cast<int *>(allocate_from_workspace(sizeof(int32_t) * checked_result.naabbs_1));
    int *prefix_sum_2 = static_cast<int *>(allocate_from_workspace(sizeof(int32_t) * checked_result.naabbs_2));
    int count_intersections_1 = PrefixSumAdd(checked_result.d_mask_1, prefix_sum_1, checked_result.naabbs_1);
    int count_intersections_2 = PrefixSumAdd(checked_result.d_mask_2, prefix_sum_2, checked_result.naabbs_2);
    int *intersection_ids_1 = static_cast<int *>(allocate_from_workspace(sizeof(int32_t) * count_intersections_1));
    int *intersection_ids_2 = static_cast<int *>(allocate_from_workspace(sizeof(int32_t) * count_intersections_2));
    auto *h_intersection_ids_1 = static_cast<int *>(std::malloc(sizeof(int32_t) * count_intersections_1));
    auto *h_intersection_ids_2 = static_cast<int *>(std::malloc(sizeof(int32_t) * count_intersections_2));
    gather_idx_kernel_host(intersection_ids_1, checked_result.naabbs_1, prefix_sum_1, checked_result.d_mask_1);
    gather_idx_kernel_host(intersection_ids_2, checked_result.naabbs_2, prefix_sum_2, checked_result.d_mask_2);
    cadMemcpy(h_intersection_ids_1, intersection_ids_1, sizeof(int32_t) * count_intersections_1, CadMemcpyKind::DeviceToHost);
    cadMemcpy(h_intersection_ids_2, intersection_ids_2, sizeof(int32_t) * count_intersections_2, CadMemcpyKind::DeviceToHost);

    printf("Round %d: count_intersections_1 = %d, count_intersections_2 = %d\n", round, count_intersections_1,
           count_intersections_2);
    if (round == 6) break;
    // FORMER `intersection.cu` adaptive step (matches legacy CUDA sources verbatim).
    if (count_intersections_1 + count_intersections_1 > 512 && global_nuv > 4) global_nuv /= 2;

    for (int i = 0; i < count_intersections_1; i++) {
      auto res = find_aabb(h_intersection_ids_1[i], aabb_1_queue);
      if (res.first < 0) {
        printf("ERROR! idx == -1\n");
        std::exit(-1);
      }
      auto &aabb = aabb_1_queue[res.first];
      float ustep = (aabb.ustop - aabb.ustart) / aabb.nuv;
      float vstep = (aabb.vstop - aabb.vstart) / aabb.nuv;
      int x = res.second % aabb.naabb_width;
      int y = res.second / aabb.naabb_width;
      EvalAndConstructTask t;
      t.d_data = t1.d_data;
      t.d_ivec = t1.d_ivec;
      t.d_ndata = t1.d_ndata;
      t.eval_task_side = 1;
      t.nuv = global_nuv;
      t.ustart = aabb.ustart + ustep * x;
      t.ustop = t.ustart + ustep;
      t.vstart = aabb.vstart + vstep * y;
      t.vstop = t.vstart + vstep;
      eval_queue.push(t);
    }

    for (int i = 0; i < count_intersections_2; i++) {
      auto res = find_aabb(h_intersection_ids_2[i], aabb_2_queue);
      if (res.first < 0) {
        printf("ERROR! idx == -1\n");
        std::exit(-1);
      }
      auto &aabb = aabb_2_queue[res.first];
      float ustep = (aabb.ustop - aabb.ustart) / aabb.nuv;
      float vstep = (aabb.vstop - aabb.vstart) / aabb.nuv;
      int x = res.second % aabb.naabb_width;
      int y = res.second / aabb.naabb_width;
      EvalAndConstructTask t;
      t.d_data = t2.d_data;
      t.d_ivec = t2.d_ivec;
      t.d_ndata = t2.d_ndata;
      t.eval_task_side = 2;
      t.nuv = global_nuv;
      t.ustart = aabb.ustart + ustep * x;
      t.ustop = t.ustart + ustep;
      t.vstart = aabb.vstart + vstep * y;
      t.vstop = t.vstart + vstep;
      eval_queue.push(t);
    }

    std::free(h_intersection_ids_1);
    std::free(h_intersection_ids_2);
    aabb_1_queue.clear();
    aabb_2_queue.clear();
  }
}
