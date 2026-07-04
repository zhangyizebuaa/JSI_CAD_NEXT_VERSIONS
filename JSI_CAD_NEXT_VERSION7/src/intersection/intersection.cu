#include "intersection.cuh"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <unordered_map>
#include <utility>

#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(JSI_CAD_EXPLICIT_SVE) && defined(__ARM_FEATURE_SVE)
#include <arm_sve.h>
#define JSI_CAD_HAS_SVE 1
#else
#define JSI_CAD_HAS_SVE 0
#endif

#include "common/common.cuh"
#include "common/profiler.hpp"
#include "evaluation/evaluation.cuh"
#include "jsi_cad.hpp"

#ifndef JSI_CAD_PRESERVE_PAIR_RELATION
#define JSI_CAD_PRESERVE_PAIR_RELATION 0
#endif

static int global_nuv = 16;
static int initial_k_prints = 0;
static float initial_u_span[3] = {1.0f, 1.0f, 1.0f};
static float initial_v_span[3] = {1.0f, 1.0f, 1.0f};

static bool use_compute_k_padding() {
  const char *value = std::getenv("CAD_INTERSECT_USE_K");
  return value == nullptr || std::strcmp(value, "0") != 0;
}

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

#if JSI_CAD_HAS_SVE
static inline int AabbOverlapHitJsSVE(int i, int j, int j_end, const float *center1, const float *radius1,
                                      const float *center2, const float *radius2, int *hit_js) {
  const svbool_t pg = svwhilelt_b32(j, j_end);
  const svfloat32_t cx1 = svdup_f32(center1[i * 3]);
  const svfloat32_t cy1 = svdup_f32(center1[i * 3 + 1]);
  const svfloat32_t cz1 = svdup_f32(center1[i * 3 + 2]);
  const svfloat32_t rx1 = svdup_f32(radius1[i * 3]);
  const svfloat32_t ry1 = svdup_f32(radius1[i * 3 + 1]);
  const svfloat32_t rz1 = svdup_f32(radius1[i * 3 + 2]);

  const svfloat32x3_t c2 = svld3_f32(pg, center2 + static_cast<size_t>(j) * 3);
  const svfloat32x3_t r2 = svld3_f32(pg, radius2 + static_cast<size_t>(j) * 3);
  const svfloat32_t cx2 = svget3_f32(c2, 0);
  const svfloat32_t cy2 = svget3_f32(c2, 1);
  const svfloat32_t cz2 = svget3_f32(c2, 2);
  const svfloat32_t rx2 = svget3_f32(r2, 0);
  const svfloat32_t ry2 = svget3_f32(r2, 1);
  const svfloat32_t rz2 = svget3_f32(r2, 2);

  svbool_t hit =
      svcmplt_f32(pg, svabs_f32_x(pg, svsub_f32_x(pg, cx1, cx2)), svadd_f32_x(pg, rx1, rx2));
  hit = svand_b_z(pg, hit,
                  svcmplt_f32(pg, svabs_f32_x(pg, svsub_f32_x(pg, cy1, cy2)),
                              svadd_f32_x(pg, ry1, ry2)));
  hit = svand_b_z(pg, hit,
                  svcmplt_f32(pg, svabs_f32_x(pg, svsub_f32_x(pg, cz1, cz2)),
                              svadd_f32_x(pg, rz1, rz2)));
  if (!svptest_any(pg, hit)) return 0;
  const int hits = static_cast<int>(svcntp_b32(pg, hit));
  svst1_s32(svwhilelt_b32(0, hits), hit_js, svcompact_s32(hit, svindex_s32(j, 1)));
  return hits;
}
#endif

static void gather_idx_kernel_host(int *data_out, size_t size, const int *pos, const int *pre) {
  for (size_t i = 0; i < size; ++i) {
    if (pre[i]) data_out[pos[i] - 1] = static_cast<int>(i);
  }
}

static void construct_aabbs_host(int un, int vn, const float *points, float curvature_padding, float *center,
                                 float *radius, float *diagonal) {
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
      radius[id * 3] = x_span / 2.0f + curvature_padding;
      radius[id * 3 + 1] = y_span / 2.0f + curvature_padding;
      radius[id * 3 + 2] = z_span / 2.0f + curvature_padding;
      x_span += 2.0f * curvature_padding;
      y_span += 2.0f * curvature_padding;
      z_span += 2.0f * curvature_padding;
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
  float curvature_padding = 0.0f;
  if (use_compute_k_padding()) {
    CAD_PROFILE_SCOPE("compute_K_intersection");
    std::vector<float> second_deriv_max(static_cast<size_t>(nuv) * nuv * 3);
    float raw_k = 0.0f;
    compute_K(results, second_deriv_max.data(), &raw_k, nuv);
    const int side = (task.eval_task_side == 2) ? 2 : 1;
    const float u_scale = std::fabs(task.ustop - task.ustart) / initial_u_span[side];
    const float v_scale = std::fabs(task.vstop - task.vstart) / initial_v_span[side];
    const float patch_scale = std::max(u_scale, v_scale);
    curvature_padding = raw_k * patch_scale * patch_scale;
    if (!std::isfinite(curvature_padding) || curvature_padding < 0.0f) curvature_padding = 0.0f;
    if (initial_k_prints < 2) {
      printf("Compute_K initial side%d raw=%g scale=%g padding=%g\n", task.eval_task_side, raw_k, patch_scale,
             curvature_padding);
      std::fflush(stdout);
      ++initial_k_prints;
    }
  }
  construct_aabbs_host(nuv, nuv, results, curvature_padding, center, radius, diagonal);
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
#if JSI_CAD_HAS_SVE
      const int vl = static_cast<int>(svcntw());
      int hit_js[64];
      for (int j = 0; j < n2; j += vl) {
        const int j_end = std::min(j + vl, n2);
        const int hits =
            AabbOverlapHitJsSVE(i, j, j_end, a1.d_aabb_center, a1.d_aabb_radius, a2.d_aabb_center, a2.d_aabb_radius,
                                hit_js);
        for (int h = 0; h < hits; ++h) {
          local_pairs.push_back({static_cast<int>(slot), i, hit_js[h]});
        }
      }
#else
      for (int j = 0; j < n2; ++j) {
        if (AabbOverlap(i, j, a1.d_aabb_center, a1.d_aabb_radius, a2.d_aabb_center, a2.d_aabb_radius)) {
          local_pairs.push_back({static_cast<int>(slot), i, j});
        }
      }
#endif
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
    EvalAndConstructTask task = tasks[k];
    out[k] = run_evaluation_and_construction_task(task);
  }
}

void intersection(EvalAndConstructTask t1, EvalAndConstructTask t2) {
  CAD_PROFILE_SCOPE("intersection");
  initial_k_prints = 0;
  initial_u_span[1] = std::max(std::fabs(t1.ustop - t1.ustart), 1.0e-12f);
  initial_v_span[1] = std::max(std::fabs(t1.vstop - t1.vstart), 1.0e-12f);
  initial_u_span[2] = std::max(std::fabs(t2.ustop - t2.ustart), 1.0e-12f);
  initial_v_span[2] = std::max(std::fabs(t2.vstop - t2.vstart), 1.0e-12f);
  IntersectChain chain;
  chain.side1.push_back(run_evaluation_and_construction_task(t1));
  chain.side2.push_back(run_evaluation_and_construction_task(t2));
  printf("Compute_K AABB padding: %s\n", use_compute_k_padding() ? "on" : "off");
  std::fflush(stdout);

  int round = 0;
  while (true) {
    ++round;
    printf("Round: %d (pair refine, aligned_slots=%zu)\n", round, chain.side1.size());
    std::vector<IntersectPairRef> pairs;
    const int npairs = intersect_collect_pairs(chain, pairs);
    printf("Round %d: overlap pairs = %d\n", round, npairs);
    std::fflush(stdout);
    if (round == 6) break;
    if (npairs == 0) break;
    if (npairs > 512 && global_nuv > 4) global_nuv /= 2;

    struct IntersectEndpoint {
      int slot;
      int bbox;
    };
    struct IntersectRefineJob {
      IntersectEndpoint side1;
      IntersectEndpoint side2;
    };
    std::vector<IntersectRefineJob> jobs;
    if constexpr (JSI_CAD_PRESERVE_PAIR_RELATION) {
      jobs.reserve(pairs.size());
      for (const auto &p : pairs) {
        jobs.push_back({{p.slot, p.i}, {p.slot, p.j}});
      }
    } else {
      std::vector<IntersectEndpoint> endpoints1;
      std::vector<IntersectEndpoint> endpoints2;
      auto add_unique = [](std::vector<IntersectEndpoint> &endpoints, int slot, int bbox) {
        for (const auto &endpoint : endpoints) {
          if (endpoint.slot == slot && endpoint.bbox == bbox) return;
        }
        endpoints.push_back({slot, bbox});
      };
      for (const auto &p : pairs) {
        add_unique(endpoints1, p.slot, p.i);
        add_unique(endpoints2, p.slot, p.j);
      }
      jobs.reserve(endpoints1.size() * endpoints2.size());
      for (const auto &endpoint1 : endpoints1) {
        for (const auto &endpoint2 : endpoints2) {
          jobs.push_back({endpoint1, endpoint2});
        }
      }
    }

    struct IntersectReuseJob {
      size_t side1_unit;
      size_t side2_unit;
    };
    std::vector<EvalAndConstructTask> tasks1;
    std::vector<EvalAndConstructTask> tasks2;
    std::vector<IntersectReuseJob> reuse_jobs;
    std::unordered_map<uint64_t, size_t> task_index1;
    std::unordered_map<uint64_t, size_t> task_index2;
    tasks1.reserve(jobs.size());
    tasks2.reserve(jobs.size());
    reuse_jobs.reserve(jobs.size());
    task_index1.reserve(jobs.size());
    task_index2.reserve(jobs.size());
    auto find_or_add_task = [](std::unordered_map<uint64_t, size_t> &task_index,
                               std::vector<EvalAndConstructTask> &tasks, const IntersectEndpoint &endpoint,
                               const EvalAndConstructTask &task) {
      const uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(endpoint.slot)) << 32) |
                           static_cast<uint32_t>(endpoint.bbox);
      const auto inserted = task_index.emplace(key, tasks.size());
      if (inserted.second) tasks.push_back(task);
      return inserted.first->second;
    };

    for (const auto &job : jobs) {
      if (job.side1.slot < 0 || job.side2.slot < 0 ||
          static_cast<size_t>(job.side1.slot) >= chain.side1.size() ||
          static_cast<size_t>(job.side2.slot) >= chain.side2.size()) {
        printf("ERROR! bad intersect endpoint slots=%d,%d\n", job.side1.slot, job.side2.slot);
        std::exit(-1);
      }
      const auto &a1 = chain.side1[static_cast<size_t>(job.side1.slot)];
      const auto &a2 = chain.side2[static_cast<size_t>(job.side2.slot)];
      const size_t side1_unit =
          find_or_add_task(task_index1, tasks1, job.side1,
                           make_intersect_refine_task(t1, a1, job.side1.bbox, 1, global_nuv));
      const size_t side2_unit =
          find_or_add_task(task_index2, tasks2, job.side2,
                           make_intersect_refine_task(t2, a2, job.side2.bbox, 2, global_nuv));
      reuse_jobs.push_back({side1_unit, side2_unit});
    }

    std::vector<AABBResult> unique_res1;
    std::vector<AABBResult> unique_res2;
    run_intersect_refine_batch(tasks1, unique_res1);
    run_intersect_refine_batch(tasks2, unique_res2);

    IntersectChain next;
    next.side1.resize(reuse_jobs.size());
    next.side2.resize(reuse_jobs.size());
    for (size_t k = 0; k < reuse_jobs.size(); ++k) {
      next.side1[k] = unique_res1[reuse_jobs[k].side1_unit];
      next.side2[k] = unique_res2[reuse_jobs[k].side2_unit];
    }
    chain = std::move(next);
  }
}
