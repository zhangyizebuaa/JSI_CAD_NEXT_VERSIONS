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

#include "evaluation/evaluation.cuh"
#include "jsi_cad.hpp"
#include "common/profiler.hpp"

static int global_nuv = 128;
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

float minimum_distance(const EvalAndConstructTask &t1, const EvalAndConstructTask &t2) {
  CAD_PROFILE_SCOPE("minimum_distance");
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  float current_min_dist = 0.0f;
  std::queue<EvalAndConstructTask> eval_queue;
  std::vector<AABBResult> aabb_1_queue;
  std::vector<AABBResult> aabb_2_queue;
  eval_queue.push(t1);
  eval_queue.push(t2);
  size_t naabbs_1 = static_cast<size_t>(t1.nuv - 1) * (t1.nuv - 1);
  size_t naabbs_2 = static_cast<size_t>(t2.nuv - 1) * (t2.nuv - 1);

  int round = 0;
  while (true) {
    ++round;
    printf("Round: %d\n", round);
    float *d_center_1 = static_cast<float *>(allocate_from_workspace(3 * naabbs_1 * sizeof(float)));
    float *d_radius_1 = static_cast<float *>(allocate_from_workspace(3 * naabbs_1 * sizeof(float)));
    float *d_diagonal_1 = static_cast<float *>(allocate_from_workspace(naabbs_1 * sizeof(float)));
    float *d_center_2 = static_cast<float *>(allocate_from_workspace(3 * naabbs_2 * sizeof(float)));
    float *d_radius_2 = static_cast<float *>(allocate_from_workspace(3 * naabbs_2 * sizeof(float)));
    float *d_diagonal_2 = static_cast<float *>(allocate_from_workspace(naabbs_2 * sizeof(float)));
    int pos_side1 = 0;
    int pos_side2 = 0;

    while (!eval_queue.empty()) {
      auto task = eval_queue.front();
      eval_queue.pop();
      auto result = run_evaluation_and_construction_task(task, cuda_streams[current_stream_id]);
      int elem_cnt = result.naabb_width * result.naabb_width;
      if (result.eval_task_side == 1) {
        aabb_1_queue.push_back(result);
        cadMemcpy(d_center_1 + 3 * pos_side1, result.d_aabb_center, elem_cnt * 3 * sizeof(float), CadMemcpyKind::DeviceToDevice);
        cadMemcpy(d_radius_1 + 3 * pos_side1, result.d_aabb_radius, elem_cnt * 3 * sizeof(float), CadMemcpyKind::DeviceToDevice);
        cadMemcpy(d_diagonal_1 + pos_side1, result.d_aabb_diagonal, elem_cnt * sizeof(float), CadMemcpyKind::DeviceToDevice);
        pos_side1 += elem_cnt;
      } else {
        aabb_2_queue.push_back(result);
        cadMemcpy(d_center_2 + 3 * pos_side2, result.d_aabb_center, elem_cnt * 3 * sizeof(float), CadMemcpyKind::DeviceToDevice);
        cadMemcpy(d_radius_2 + 3 * pos_side2, result.d_aabb_radius, elem_cnt * 3 * sizeof(float), CadMemcpyKind::DeviceToDevice);
        cadMemcpy(d_diagonal_2 + pos_side2, result.d_aabb_diagonal, elem_cnt * sizeof(float), CadMemcpyKind::DeviceToDevice);
        pos_side2 += elem_cnt;
      }
      current_stream_id = (current_stream_id + 1) % num_cuda_streams;
    }

    CheckDistanceTask check_distance;
    check_distance.aabb_cnt1 = naabbs_1;
    check_distance.aabb_cnt2 = naabbs_2;
    check_distance.d_aabb_layer1_center = d_center_1;
    check_distance.d_aabb_layer2_center = d_center_2;
    check_distance.d_aabb_layer1_radius = d_radius_1;
    check_distance.d_aabb_layer2_radius = d_radius_2;
    check_distance.d_aabb_layer1_diagonal = d_diagonal_1;
    check_distance.d_aabb_layer2_diagonal = d_diagonal_2;
    auto checked_result = check_distance_of_aabbs(check_distance);
    if (checked_result.naabbs_1 == 0 || checked_result.naabbs_2 == 0) break;
    current_min_dist = checked_result.h_kth_distance;
    printf("[current_min_dist] %lf\n", checked_result.h_kth_distance);

    int *d_prefix_sum_1 =
        static_cast<int *>(allocate_from_workspace(sizeof(int32_t) * checked_result.naabbs_1));
    int *d_prefix_sum_2 =
        static_cast<int *>(allocate_from_workspace(sizeof(int32_t) * checked_result.naabbs_2));
    int count_possible_dist_1 =
        std::min(PrefixSumAdd(checked_result.d_mask_1, d_prefix_sum_1, checked_result.naabbs_1), k_of_topk);
    int count_possible_dist_2 =
        std::min(PrefixSumAdd(checked_result.d_mask_2, d_prefix_sum_2, checked_result.naabbs_2), k_of_topk);
    int *d_dist_ids_1 = static_cast<int *>(allocate_from_workspace(sizeof(int32_t) * count_possible_dist_1));
    int *d_dist_ids_2 = static_cast<int *>(allocate_from_workspace(sizeof(int32_t) * count_possible_dist_2));
    auto *h_dist_ids_1 = static_cast<int *>(std::malloc(sizeof(int32_t) * count_possible_dist_1));
    auto *h_dist_ids_2 = static_cast<int *>(std::malloc(sizeof(int32_t) * count_possible_dist_2));
    gather_idx_kernel_host(d_dist_ids_1, checked_result.naabbs_1, k_of_topk, d_prefix_sum_1, checked_result.d_mask_1);
    gather_idx_kernel_host(d_dist_ids_2, checked_result.naabbs_2, k_of_topk, d_prefix_sum_2, checked_result.d_mask_2);
    cadMemcpy(h_dist_ids_1, d_dist_ids_1, sizeof(int) * count_possible_dist_1, CadMemcpyKind::DeviceToHost);
    cadMemcpy(h_dist_ids_2, d_dist_ids_2, sizeof(int) * count_possible_dist_2, CadMemcpyKind::DeviceToHost);

    if (round == 1) global_nuv = next_nuv;
    if (round == max_layers) break;

    size_t new_naabbs_1 = 0;
    for (int i = 0; i < count_possible_dist_1; i++) {
      auto res = find_aabb(h_dist_ids_1[i], aabb_1_queue);
      if (res.first == static_cast<size_t>(-1)) {
        printf("ERROR! idx == -1 for %d\n", h_dist_ids_1[i]);
        std::exit(-1);
      }
      auto &aabb = aabb_1_queue[res.first];
      int bbox_id = static_cast<int>(res.second);
      float ustep = (aabb.ustop - aabb.ustart) / aabb.nuv;
      float vstep = (aabb.vstop - aabb.vstart) / aabb.nuv;
      int x = bbox_id % aabb.naabb_width;
      int y = bbox_id / aabb.naabb_width;
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
      new_naabbs_1 += static_cast<size_t>(t.nuv - 1) * (t.nuv - 1);
      eval_queue.push(t);
    }

    size_t new_naabbs_2 = 0;
    for (int i = 0; i < count_possible_dist_2; i++) {
      auto res = find_aabb(h_dist_ids_2[i], aabb_2_queue);
      if (res.first == static_cast<size_t>(-1)) {
        printf("ERROR! idx == -1\n");
        std::exit(-1);
      }
      auto &aabb = aabb_2_queue[res.first];
      int bbox_id = static_cast<int>(res.second);
      float ustep = (aabb.ustop - aabb.ustart) / aabb.nuv;
      float vstep = (aabb.vstop - aabb.vstart) / aabb.nuv;
      int x = bbox_id % aabb.naabb_width;
      int y = bbox_id / aabb.naabb_width;
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
      new_naabbs_2 += static_cast<size_t>(t.nuv - 1) * (t.nuv - 1);
      eval_queue.push(t);
    }

    std::free(h_dist_ids_1);
    std::free(h_dist_ids_2);
    aabb_1_queue.clear();
    aabb_2_queue.clear();
    naabbs_1 = new_naabbs_1;
    naabbs_2 = new_naabbs_2;
    if (eval_queue.empty()) break;
  }

  if (use_time) {
    auto end = std::chrono::steady_clock::now();
    printf("Elapsed: %lfms\n", std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000.0);
  }
  return current_min_dist;
}
