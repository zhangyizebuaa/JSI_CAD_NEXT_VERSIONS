#include "distance_analysis.cuh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
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

#include "evaluation/evaluation.cuh"
#include "jsi_cad.hpp"
#include "common/profiler.hpp"

static int global_nuv = 64;
static int next_nuv = 16;
static int k_of_topk = 8;
static int max_layers = 4;
static bool use_time = false;

static void gather_idx_kernel_host(int *data_out, size_t size, int max_num, const int *pos, const int *pre) {
  for (size_t i = 0; i < size; ++i) {
    if (pre[i] && pre[i] < max_num) data_out[pos[i] - 1] = static_cast<int>(i);
  }
}

static void construct_aabbs_distance_host(int un, int vn, const float *points, float *center, float *radius,
                                          float *diagonal) {
  // Matches FORMER `construct_aabbs_distance` kernel (including k margin); k==0 matches legacy behavior.
  const float k = 0.0f;
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
      float x_span = xmax - xmin + 2.0f * k;
      float y_span = ymax - ymin + 2.0f * k;
      float z_span = zmax - zmin + 2.0f * k;
      center[id * 3] = (xmin - k + xmax + k) / 2.0f;
      center[id * 3 + 1] = (ymin - k + ymax + k) / 2.0f;
      center[id * 3 + 2] = (zmin - k + zmax + k) / 2.0f;
      radius[id * 3] = x_span / 2.0f;
      radius[id * 3 + 1] = y_span / 2.0f;
      radius[id * 3 + 2] = z_span / 2.0f;
      diagonal[id] = std::sqrt(x_span * x_span + y_span * y_span + z_span * z_span);
    }
  }
}

AABBResult run_evaluation_and_construction_task(EvalAndConstructTask &task, CadStream stream) {
  CAD_PROFILE_SCOPE("run_evaluation_and_construction_task_dist");
  (void)stream;
  int nuv = task.nuv;
  float ustep = (task.ustop - task.ustart) / nuv;
  float vstep = (task.vstop - task.vstart) / nuv;
  int naabb_width = nuv - 1;

  float *results = static_cast<float *>(allocate_from_workspace(ONE_UV_SIZE * nuv * nuv * sizeof(float)));
  float *center = static_cast<float *>(allocate_from_workspace(3 * naabb_width * naabb_width * sizeof(float)));
  float *radius = static_cast<float *>(allocate_from_workspace(3 * naabb_width * naabb_width * sizeof(float)));
  float *diagonal = static_cast<float *>(allocate_from_workspace(naabb_width * naabb_width * sizeof(float)));

  EG_spline2dDeriv_irrational(task.d_ivec, task.d_data, ustep, vstep, task.ustart, task.vstart, nuv, results);
  construct_aabbs_distance_host(nuv, nuv, results, center, radius, diagonal);

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

CheckDistanceResult check_distance_of_aabbs(CheckDistanceTask t) {
  CAD_PROFILE_SCOPE("check_distance_of_aabbs");
  size_t num_elements = t.aabb_cnt1 * t.aabb_cnt2;
  int *mask1 = static_cast<int *>(allocate_from_workspace(t.aabb_cnt1 * sizeof(int)));
  int *mask2 = static_cast<int *>(allocate_from_workspace(t.aabb_cnt2 * sizeof(int)));
  float *distance_ab = static_cast<float *>(allocate_from_workspace(num_elements * sizeof(float)));
  cadMemset(mask1, 0, t.aabb_cnt1 * sizeof(int));
  cadMemset(mask2, 0, t.aabb_cnt2 * sizeof(int));

  float kth_distance = 0.0f;
  // FORMER path: host fills distances then CUB radix sort equivalent (no Ascend shortcut).

  // Pairwise center L2 distance; same expression as FORMER `para_check_aabbs_distance` kernel.
  for (long long idx = 0; idx < static_cast<long long>(num_elements); ++idx) {
    size_t i = static_cast<size_t>(idx) / t.aabb_cnt2;
    size_t j = static_cast<size_t>(idx) % t.aabb_cnt2;
    float dx = std::fabs(t.d_aabb_layer1_center[i * 3] - t.d_aabb_layer2_center[j * 3]);
    float dy = std::fabs(t.d_aabb_layer1_center[i * 3 + 1] - t.d_aabb_layer2_center[j * 3 + 1]);
    float dz = std::fabs(t.d_aabb_layer1_center[i * 3 + 2] - t.d_aabb_layer2_center[j * 3 + 2]);
    distance_ab[idx] = std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  // FORMER: CUB DeviceRadixSort ascending, then k-th key. Same k-th order statistic as sorting scratch[k].
  std::vector<float> scratch(distance_ab, distance_ab + num_elements);
  size_t k = std::min(static_cast<size_t>(k_of_topk), num_elements) - 1;
  const auto k_end = scratch.begin() + static_cast<std::ptrdiff_t>(k) + 1;
  std::partial_sort(scratch.begin(), k_end, scratch.end());
  kth_distance = scratch[k];

  // Same semantics as FORMER `para_mark_mask` (atomicOr): mark side if any pair along the other axis is within kth.
  for (long long i = 0; i < static_cast<long long>(t.aabb_cnt1); ++i) {
    int hit = 0;
    size_t base = static_cast<size_t>(i) * t.aabb_cnt2;
    for (size_t j = 0; j < t.aabb_cnt2; ++j) {
      if (distance_ab[base + j] <= kth_distance) {
        hit = 1;
        break;
      }
    }
    mask1[i] = hit;
  }

  for (long long j = 0; j < static_cast<long long>(t.aabb_cnt2); ++j) {
    int hit = 0;
    for (size_t i = 0; i < t.aabb_cnt1; ++i) {
      if (distance_ab[i * t.aabb_cnt2 + static_cast<size_t>(j)] <= kth_distance) {
        hit = 1;
        break;
      }
    }
    mask2[j] = hit;
  }

  CheckDistanceResult res;
  res.d_mask_1 = mask1;
  res.d_mask_2 = mask2;
  res.naabbs_1 = t.aabb_cnt1;
  res.naabbs_2 = t.aabb_cnt2;
  res.h_kth_distance = kth_distance;
  return res;
}

static std::pair<size_t, size_t> find_aabb(size_t id, std::vector<AABBResult> &aabbs) {
  size_t total = 0;
  for (size_t i = 0; i < aabbs.size(); i++) {
    auto &aabb = aabbs[i];
    total += static_cast<size_t>(aabb.naabb_width) * aabb.naabb_width;
    if (total > id) return {i, id - (total - static_cast<size_t>(aabb.naabb_width) * aabb.naabb_width)};
  }
  return {static_cast<size_t>(-1), static_cast<size_t>(-1)};
}

struct DistPairRef {
  float d2;
  int slot;
  int i;
  int j;
};

struct DistChain {
  std::vector<AABBResult> side1;
  std::vector<AABBResult> side2;
};

static inline float center_distance_sq(size_t i, size_t j, const float *c1, const float *c2) {
  const float dx = c1[i * 3] - c2[j * 3];
  const float dy = c1[i * 3 + 1] - c2[j * 3 + 1];
  const float dz = c1[i * 3 + 2] - c2[j * 3 + 2];
  return dx * dx + dy * dy + dz * dz;
}

#if JSI_CAD_HAS_SVE
static inline svfloat32_t center_distance_sq_sve(size_t i, int j, const float *c1, const float *c2, svbool_t pg) {
  const svint32_t x_idx = svindex_s32(j * 3, 3);
  const svint32_t y_idx = svindex_s32(j * 3 + 1, 3);
  const svint32_t z_idx = svindex_s32(j * 3 + 2, 3);
  const svfloat32_t x1 = svdup_f32(c1[i * 3]);
  const svfloat32_t y1 = svdup_f32(c1[i * 3 + 1]);
  const svfloat32_t z1 = svdup_f32(c1[i * 3 + 2]);
  const svfloat32_t dx = svsub_f32_x(pg, x1, svld1_gather_s32index_f32(pg, c2, x_idx));
  const svfloat32_t dy = svsub_f32_x(pg, y1, svld1_gather_s32index_f32(pg, c2, y_idx));
  const svfloat32_t dz = svsub_f32_x(pg, z1, svld1_gather_s32index_f32(pg, c2, z_idx));
  svfloat32_t d2 = svmul_f32_x(pg, dx, dx);
  d2 = svadd_f32_x(pg, d2, svmul_f32_x(pg, dy, dy));
  return svadd_f32_x(pg, d2, svmul_f32_x(pg, dz, dz));
}
#endif

static bool dist_pair_ref_less(const DistPairRef &a, const DistPairRef &b) {
  if (a.d2 != b.d2) return a.d2 < b.d2;
  if (a.slot != b.slot) return a.slot < b.slot;
  if (a.i != b.i) return a.i < b.i;
  return a.j < b.j;
}

static void dist_pair_top_push(std::vector<DistPairRef> &best, int limit, float d2, int slot, int i, int j) {
  if (limit <= 0) return;
  DistPairRef p{d2, slot, i, j};
  if (static_cast<int>(best.size()) < limit) {
    best.push_back(p);
    std::sort(best.begin(), best.end(), dist_pair_ref_less);
    return;
  }
  if (!dist_pair_ref_less(p, best.back())) return;
  best.back() = p;
  std::sort(best.begin(), best.end(), dist_pair_ref_less);
}

static float dist_collect_top_pairs(const DistChain &ch, std::vector<DistPairRef> &pairs) {
  pairs.clear();
  if (ch.side1.empty() || ch.side1.size() != ch.side2.size()) return 0.0f;
  const int limit = std::max(1, k_of_topk);
#ifdef _OPENMP
  struct DistRowTask {
    int slot;
    int i;
  };
  std::vector<DistRowTask> rows;
  size_t total_pairs = 0;
  for (size_t slot = 0; slot < ch.side1.size(); ++slot) {
    const int n1 = ch.side1[slot].naabb_width * ch.side1[slot].naabb_width;
    const int n2 = ch.side2[slot].naabb_width * ch.side2[slot].naabb_width;
    total_pairs += static_cast<size_t>(n1) * static_cast<size_t>(n2);
    rows.reserve(rows.size() + static_cast<size_t>(n1));
    for (int i = 0; i < n1; ++i) rows.push_back({static_cast<int>(slot), i});
  }
  if (rows.empty()) return 0.0f;
  std::vector<std::vector<DistPairRef>> thread_bests;
#pragma omp parallel if (total_pairs > 4096)
  {
    const int tid = omp_get_thread_num();
#pragma omp single
    { thread_bests.resize(static_cast<size_t>(omp_get_num_threads())); }
    auto &local = thread_bests[static_cast<size_t>(tid)];
    local.reserve(static_cast<size_t>(limit));
#pragma omp for schedule(static) nowait
    for (size_t row = 0; row < rows.size(); ++row) {
      const int slot = rows[row].slot;
      const int i = rows[row].i;
      const auto &a1 = ch.side1[static_cast<size_t>(slot)];
      const auto &a2 = ch.side2[static_cast<size_t>(slot)];
      const int n2 = a2.naabb_width * a2.naabb_width;
#if JSI_CAD_HAS_SVE
      const int vl = static_cast<int>(svcntw());
      alignas(16) float d2_lanes[64];
      alignas(16) int j_lanes[64];
      for (int j = 0; j < n2; j += vl) {
        const svbool_t pg = svwhilelt_b32(j, n2);
        const svfloat32_t d2v =
            center_distance_sq_sve(static_cast<size_t>(i), j, a1.d_aabb_center, a2.d_aabb_center, pg);
        svbool_t hit = pg;
        if (static_cast<int>(local.size()) >= limit) {
          hit = svcmplt_f32(pg, d2v, svdup_f32(local.back().d2));
          if (!svptest_any(pg, hit)) continue;
        }
        const int hits = static_cast<int>(svcntp_b32(pg, hit));
        const svbool_t packed_pg = svwhilelt_b32(0, hits);
        svst1_f32(packed_pg, d2_lanes, svcompact_f32(hit, d2v));
        svst1_s32(packed_pg, j_lanes, svcompact_s32(hit, svindex_s32(j, 1)));
        for (int h = 0; h < hits; ++h) {
          dist_pair_top_push(local, limit, d2_lanes[h], slot, i, j_lanes[h]);
        }
      }
#else
      for (int j = 0; j < n2; ++j) {
        dist_pair_top_push(local, limit,
                           center_distance_sq(static_cast<size_t>(i), static_cast<size_t>(j), a1.d_aabb_center,
                                              a2.d_aabb_center),
                           slot, i, j);
      }
#endif
    }
  }
  std::vector<DistPairRef> best;
  best.reserve(static_cast<size_t>(limit));
  for (auto &local : thread_bests) {
    for (const auto &p : local) {
      dist_pair_top_push(best, limit, p.d2, p.slot, p.i, p.j);
    }
  }
#else
  std::vector<DistPairRef> best;
  best.reserve(static_cast<size_t>(limit));
  for (size_t slot = 0; slot < ch.side1.size(); ++slot) {
    const auto &a1 = ch.side1[slot];
    const auto &a2 = ch.side2[slot];
    const int n1 = a1.naabb_width * a1.naabb_width;
    const int n2 = a2.naabb_width * a2.naabb_width;
    for (int i = 0; i < n1; ++i) {
      for (int j = 0; j < n2; ++j) {
        dist_pair_top_push(best, limit,
                           center_distance_sq(static_cast<size_t>(i), static_cast<size_t>(j), a1.d_aabb_center,
                                              a2.d_aabb_center),
                           static_cast<int>(slot), i, j);
      }
    }
  }
#endif
  pairs = best;
  return pairs.empty() ? 0.0f : std::sqrt(pairs.back().d2);
}

static EvalAndConstructTask make_dist_refine_task(const EvalAndConstructTask &tmpl, const AABBResult &aabb,
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

float minimum_distance(const EvalAndConstructTask &t1, const EvalAndConstructTask &t2) {
  CAD_PROFILE_SCOPE("minimum_distance");
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::vector<DistChain> chains(1);
  EvalAndConstructTask a = t1;
  EvalAndConstructTask b = t2;
  chains[0].side1.push_back(run_evaluation_and_construction_task(a, cuda_streams[current_stream_id]));
  current_stream_id = (current_stream_id + 1) % num_cuda_streams;
  chains[0].side2.push_back(run_evaluation_and_construction_task(b, cuda_streams[current_stream_id]));
  current_stream_id = (current_stream_id + 1) % num_cuda_streams;

  float current_min_dist = 0.0f;
  int round = 0;
  while (true) {
    ++round;
    auto &ch = chains[0];
    printf("Round: %d (pair refine, aligned_slots=%zu)\n", round, ch.side1.size());
    std::vector<DistPairRef> pairs;
    current_min_dist = dist_collect_top_pairs(ch, pairs);
    printf("[current_min_dist] %lf\n", current_min_dist);
    printf("  pairs to refine: %zu\n", pairs.size());

    if (pairs.empty()) break;
    if (round == max_layers) break;
    if (round == 1) global_nuv = next_nuv;

    struct DistRefineKey {
      int slot;
      int bbox;
    };
    struct DistRefineUnit {
      DistRefineKey key;
      EvalAndConstructTask task;
    };
    struct DistRefineJob {
      int side1_unit;
      int side2_unit;
    };
    auto find_or_add_unit = [](std::vector<DistRefineUnit> &units, const DistRefineKey &key,
                               const EvalAndConstructTask &task) {
      for (size_t idx = 0; idx < units.size(); ++idx) {
        if (units[idx].key.slot == key.slot && units[idx].key.bbox == key.bbox) {
          return static_cast<int>(idx);
        }
      }
      units.push_back({key, task});
      return static_cast<int>(units.size() - 1);
    };

    std::vector<DistRefineUnit> unique1;
    std::vector<DistRefineUnit> unique2;
    std::vector<DistRefineJob> jobs;
    unique1.reserve(pairs.size());
    unique2.reserve(pairs.size());
    jobs.reserve(pairs.size());
    for (const auto &p : pairs) {
      if (p.slot < 0 || static_cast<size_t>(p.slot) >= ch.side1.size()) {
        printf("ERROR! bad distance pair slot=%d\n", p.slot);
        std::exit(-1);
      }
      const auto &a1 = ch.side1[static_cast<size_t>(p.slot)];
      const auto &a2 = ch.side2[static_cast<size_t>(p.slot)];
      const int n1 = a1.naabb_width * a1.naabb_width;
      const int n2 = a2.naabb_width * a2.naabb_width;
      if (p.i < 0 || p.i >= n1 || p.j < 0 || p.j >= n2) {
        printf("ERROR! bad distance pair slot=%d i=%d j=%d\n", p.slot, p.i, p.j);
        std::exit(-1);
      }
      const int side1_unit =
          find_or_add_unit(unique1, DistRefineKey{p.slot, p.i}, make_dist_refine_task(t1, a1, p.i, 1, global_nuv));
      const int side2_unit =
          find_or_add_unit(unique2, DistRefineKey{p.slot, p.j}, make_dist_refine_task(t2, a2, p.j, 2, global_nuv));
      jobs.push_back({side1_unit, side2_unit});
    }

    std::vector<AABBResult> unique_res1(unique1.size());
    std::vector<AABBResult> unique_res2(unique2.size());
    for (size_t k = 0; k < unique1.size(); ++k) {
      EvalAndConstructTask task = unique1[k].task;
      unique_res1[k] = run_evaluation_and_construction_task(task, cuda_streams[current_stream_id]);
      current_stream_id = (current_stream_id + 1) % num_cuda_streams;
    }
    for (size_t k = 0; k < unique2.size(); ++k) {
      EvalAndConstructTask task = unique2[k].task;
      unique_res2[k] = run_evaluation_and_construction_task(task, cuda_streams[current_stream_id]);
      current_stream_id = (current_stream_id + 1) % num_cuda_streams;
    }

    DistChain next;
    next.side1.resize(jobs.size());
    next.side2.resize(jobs.size());
    for (size_t k = 0; k < jobs.size(); ++k) {
      next.side1[k] = unique_res1[static_cast<size_t>(jobs[k].side1_unit)];
      next.side2[k] = unique_res2[static_cast<size_t>(jobs[k].side2_unit)];
    }
    chains[0] = std::move(next);
  }

  if (use_time) {
    auto end = std::chrono::steady_clock::now();
    printf("Elapsed: %lfms\n", std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000.0);
  }
  return current_min_dist;
}
