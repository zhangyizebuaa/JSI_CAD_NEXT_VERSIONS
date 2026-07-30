#ifndef _JSI_CAD_HPP
#define _JSI_CAD_HPP

#include <string>
#include <cstdio>
#include <cstdlib>

using cudaError_t = int;
inline constexpr int cudaSuccess = 0;

struct EvalAndConstructTask {
    int *d_ivec;
    float *d_data;
    int d_ndata;

    float ustart;
    float vstart;
    float ustop;
    float vstop;
    int nuv;

    int eval_task_side;
};

void init_cad();

void deinit_cad();

void* allocate_from_workspace(size_t size);

float minimum_distance(const EvalAndConstructTask &t1, const EvalAndConstructTask &t2);

void intersection(EvalAndConstructTask t1, EvalAndConstructTask t2);

static int default_global_dist_nuv = 512;
static int default_global_intersect_nuv = 16;

inline void checkError(cudaError_t error, std::string msg) {
    if (error != cudaSuccess) {
        printf("%s: %d\n", msg.c_str(), error);
        exit(1);
    }
}

#define ONE_UV_SIZE 18
void evalWithCuda(int *ivec, float *data, int ndata, float ustep, float vstep, float ustart, float vstart, int nuv, float *results);

#endif
