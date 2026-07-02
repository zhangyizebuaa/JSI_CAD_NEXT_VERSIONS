#include <atomic>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef USE_ASCEND_ACL
#include <acl/acl.h>
#endif

#include "common.cuh"
#include "jsi_cad.hpp"

void *workspace = nullptr;
size_t workspace_size = 4LL * 1024 * 1024 * 1024;
std::atomic<size_t> workspace_offset(0);

int current_stream_id = 0;
CadStream cuda_streams[num_cuda_streams] = {0};
static bool g_acl_enabled = false;

void *allocate_from_workspace(size_t size) {
  if (size % 256 != 0) size += 256 - size % 256;
  size_t current_offset = workspace_offset.fetch_add(size);
  if (current_offset + size > workspace_size) {
    printf("Error: Out of workspace memory!\n");
    return nullptr;
  }
  return static_cast<char *>(workspace) + current_offset;
}

void init_cad() {
#ifdef USE_ASCEND_ACL
  if (aclInit(nullptr) == ACL_ERROR_NONE && aclrtSetDevice(0) == ACL_ERROR_NONE) {
    g_acl_enabled = true;
  } else {
    g_acl_enabled = false;
    printf("[warn] ACL init failed, fallback to host runtime\n");
  }
#endif
  if (cadMalloc(&workspace, workspace_size) != 0 || workspace == nullptr) {
    printf("alloc workspace failed\n");
    std::exit(1);
  }
  workspace_offset.store(0);
}

void deinit_cad() {
  cadFree(workspace);
  workspace = nullptr;
  workspace_offset.store(0);
#ifdef USE_ASCEND_ACL
  if (g_acl_enabled) {
    aclrtResetDevice(0);
    aclFinalize();
  }
#endif
}

int cadMemcpy(void *dst, const void *src, size_t size, CadMemcpyKind kind) {
  if (size == 0) return 0;
#ifdef USE_ASCEND_ACL
  if (g_acl_enabled) {
    aclrtMemcpyKind acl_kind = ACL_MEMCPY_HOST_TO_HOST;
    if (kind == CadMemcpyKind::HostToDevice) {
      acl_kind = ACL_MEMCPY_HOST_TO_DEVICE;
    } else if (kind == CadMemcpyKind::DeviceToHost) {
      acl_kind = ACL_MEMCPY_DEVICE_TO_HOST;
    } else if (kind == CadMemcpyKind::DeviceToDevice) {
      acl_kind = ACL_MEMCPY_DEVICE_TO_DEVICE;
    }
    // Current migration keeps host-accessible workspace. If ACL memcpy
    // cannot handle these pointers, gracefully fallback to host memcpy.
    aclError ret = aclrtMemcpy(dst, size, src, size, acl_kind);
    if (ret == ACL_ERROR_NONE) return 0;
  }
#else
  (void)kind;
#endif
  std::memcpy(dst, src, size);
  return 0;
}

int cadMemset(void *ptr, int value, size_t size) {
  if (size == 0) return 0;
#ifdef USE_ASCEND_ACL
  if (g_acl_enabled) {
    // Fallback to host memset if ACL rejects non-device pointers.
    aclError ret = aclrtMemset(ptr, size, value, size);
    if (ret == ACL_ERROR_NONE) return 0;
  }
#endif
  std::memset(ptr, value, size);
  return 0;
}

int cadMalloc(void **ptr, size_t size) {
  if (ptr == nullptr) return -1;
  if (size == 0) {
    *ptr = nullptr;
    return 0;
  }
#ifdef USE_ASCEND_ACL
  if (g_acl_enabled) {
    aclError ret = aclrtMallocHost(ptr, size);
    if (ret == ACL_ERROR_NONE) return 0;
  }
#endif
  *ptr = std::malloc(size);
  return (*ptr == nullptr) ? -1 : 0;
}

int cadFree(void *ptr) {
  if (ptr == nullptr) return 0;
#ifdef USE_ASCEND_ACL
  if (g_acl_enabled) {
    aclError ret = aclrtFreeHost(ptr);
    if (ret == ACL_ERROR_NONE) return 0;
  }
#endif
  std::free(ptr);
  return 0;
}

bool cadIsAclEnabled() {
  return g_acl_enabled;
}

template <typename T>
inline T CadMax(T a, T b) {
  return a > b ? a : b;
}

int PrefixSumAdd(const int32_t *d_input, int32_t *d_output, size_t size) {
  int32_t acc = 0;
  for (size_t i = 0; i < size; ++i) {
    acc += d_input[i];
    d_output[i] = acc;
  }
  return size == 0 ? 0 : d_output[size - 1];
}

// FORMER used `get_SecondDeriv` + multi-stage `reduce_Max` + `reduce_M_K` on device. Max is associative, so
// taking the global max of per-point extrema matches the reduced tree; final scalar matches `reduce_M_K` when
// the reduced grid is 1x1 (see FORMER common.cu).
void compute_K(float *d_evalresults, float *d_deriv, float *d_k, int n) {
  float max_uu = -FLT_MAX;
  float max_uv = -FLT_MAX;
  float max_vv = -FLT_MAX;
  int point_cnt = n * n;

  for (int i = 0; i < point_cnt; ++i) {
    float uu = CadMax(std::fabs(d_evalresults[i * 18 + 9]),
                      CadMax(std::fabs(d_evalresults[i * 18 + 10]), std::fabs(d_evalresults[i * 18 + 11])));
    float uv = CadMax(std::fabs(d_evalresults[i * 18 + 12]),
                      CadMax(std::fabs(d_evalresults[i * 18 + 13]), std::fabs(d_evalresults[i * 18 + 14])));
    float vv = CadMax(std::fabs(d_evalresults[i * 18 + 15]),
                      CadMax(std::fabs(d_evalresults[i * 18 + 16]), std::fabs(d_evalresults[i * 18 + 17])));
    d_deriv[i * 3] = uu;
    d_deriv[i * 3 + 1] = uv;
    d_deriv[i * 3 + 2] = vv;
    max_uu = CadMax(max_uu, uu);
    max_uv = CadMax(max_uv, uv);
    max_vv = CadMax(max_vv, vv);
  }

  int side = n > 1 ? (n - 1) : 1;
  *d_k = (max_uu / (side * side) + 2.0f * max_uv / (side * side) + max_vv / (side * side)) / 8.0f;
}
