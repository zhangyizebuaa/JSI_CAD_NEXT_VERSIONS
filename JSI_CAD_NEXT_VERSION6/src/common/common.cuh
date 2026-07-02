#ifndef _COMMON_UTILS_HPP
#define _COMMON_UTILS_HPP

#include <string>
#include <cstdio>

#include <cstddef>
#include <cstdint>

using CadStream = int;
enum class CadMemcpyKind { HostToDevice, DeviceToHost, DeviceToDevice };

inline constexpr int num_cuda_streams = 4;
extern int current_stream_id;
extern CadStream cuda_streams[num_cuda_streams];

void* allocate_from_workspace(size_t size);

void init_cad();

void deinit_cad();

int cadMalloc(void **ptr, size_t size);
int cadFree(void *ptr);
int cadMemcpy(void *dst, const void *src, size_t size, CadMemcpyKind kind);
int cadMemset(void *ptr, int value, size_t size);
bool cadIsAclEnabled();

int PrefixSumAdd(const int32_t* d_input, int32_t* d_output, size_t size);

void compute_K(float *d_evalresults, //evaluation results
                float *d_deriv,    //max(dxuu,dyuu,duu)
                                    //max(dxuv,dyuv,dzuv)
                                    //max(dxvv, dyvv, dzvv)
                float *d_out, int n); //(M1, M2, M3)
#endif