#include "intersection.cuh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <utility>

#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

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
  EG_spline2dDeriv_irrational_threshold(task.d_ivec, task.d_data, ustep, vstep, task.ustart, task.vstart, nuv,
                                        results, 512);
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

struct IntersectPairRef {
  int slot;
  int i;
  int j;
};

struct IntersectChain {
  std::vector<AABBResult> side1;
  std::vector<AABBResult> side2;
};

static int intersect_collect_pairs(const IntersectChain &ch, std::vector<IntersectPairRef> &pairs) {
  pairs.clear();
  if (ch.side1.empty() || ch.side1.size() != ch.side2.size()) return 0;
#ifdef _OPENMP
  std::vector<std::vector<IntersectPairRef>> slot_pairs(ch.side1.size());
#pragma omp parallel for schedule(dynamic) if (ch.side1.size() > 4)
  for (long long slot_ll = 0; slot_ll < static_cast<long long>(ch.side1.size()); ++slot_ll) {
    const size_t slot = static_cast<size_t>(slot_ll);
    auto &local_pairs = slot_pairs[slot];
#else
  for (size_t slot = 0; slot < ch.side1.size(); ++slot) {
    auto &local_pairs = pairs;
#endif
    const auto &a1 = ch.side1[slot];
    const auto &a2 = ch.side2[slot];
    const int n1 = a1.naabb_width * a1.naabb_width;
    const int n2 = a2.naabb_width * a2.naabb_width;
    for (int i = 0; i < n1; ++i) {
      for (int j = 0; j < n2; ++j) {
        if (AabbOverlap(i, j, a1.d_aabb_center, a1.d_aabb_radius, a2.d_aabb_center, a2.d_aabb_radius)) {
          local_pairs.push_back({static_cast<int>(slot), i, j});
        }
      }
    }
  }
#ifdef _OPENMP
  size_t total = 0;
  for (const auto &local : slot_pairs) {
    total += local.size();
  }
  pairs.reserve(total);
  for (auto &local : slot_pairs) {
    pairs.insert(pairs.end(), local.begin(), local.end());
  }
#endif
  return static_cast<int>(pairs.size());
}

static EvalAndConstructTask make_intersect_refine_task(const EvalAndConstructTask &tmpl, const AABBResult &aabb,
                                                       int bbox_id, int side, int nuv) {
  const float ustep = (aabb.ustop - aabb.ustart) / aabb.nuv;
  const float vstep = (aabb.vstop - aabb.vstart) / aabb.nuv;
  const int x = bbox_id % aabb.naabb_width;
  const int y = bbox_id / aabb.naabb_width;
  EvalAndConstructTask t = tmpl;
  t.eval_task_side = side;
  t.nuv = nuv;
  t.ustart = aabb.ustart + ustep * x;
  t.ustop = t.ustart + ustep;
  t.vstart = aabb.vstart + vstep * y;
  t.vstop = t.vstart + vstep;
  return t;
}

static void run_intersect_refine_batch(const std::vector<EvalAndConstructTask> &tasks, std::vector<AABBResult> &out) {
  out.resize(tasks.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if (tasks.size() > 4)
#endif
  for (size_t k = 0; k < tasks.size(); ++k) {
    EvalAndConstructTask t = tasks[k];
    out[k] = run_evaluation_and_construction_task(t);
  }
}

void intersection(EvalAndConstructTask t1, EvalAndConstructTask t2) {
  CAD_PROFILE_SCOPE("intersection");
  IntersectChain chain;
  chain.side1.push_back(run_evaluation_and_construction_task(t1));
  chain.side2.push_back(run_evaluation_and_construction_task(t2));

  int round = 0;
  while (true) {
    ++round;
    printf("Round: %d (pair refine, aligned_slots=%zu)\n", round, chain.side1.size());
    std::vector<IntersectPairRef> pairs;
    const int npairs = intersect_collect_pairs(chain, pairs);
    printf("Round %d: overlap pairs = %d\n", round, npairs);
    if (round == 6) break;
    if (npairs == 0) break;
    if (npairs > 512 && global_nuv > 4) global_nuv /= 2;

    IntersectChain next;
    std::vector<EvalAndConstructTask> tasks1;
    std::vector<EvalAndConstructTask> tasks2;
    tasks1.reserve(pairs.size());
    tasks2.reserve(pairs.size());
    for (size_t k = 0; k < pairs.size(); ++k) {
      const auto &p = pairs[k];
      if (p.slot < 0 || static_cast<size_t>(p.slot) >= chain.side1.size()) {
        printf("ERROR! bad intersect pair slot=%d\n", p.slot);
        std::exit(-1);
      }
      tasks1.push_back(make_intersect_refine_task(t1, chain.side1[static_cast<size_t>(p.slot)], p.i, 1, global_nuv));
      tasks2.push_back(make_intersect_refine_task(t2, chain.side2[static_cast<size_t>(p.slot)], p.j, 2, global_nuv));
    }
    run_intersect_refine_batch(tasks1, next.side1);
    run_intersect_refine_batch(tasks2, next.side2);
    chain = std::move(next);
  }
}
