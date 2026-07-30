#include <iostream>
#include <chrono>
#include <cstdlib>
#include <cstdio>

#include "jsi_cad.hpp"
#include "nurbs_data.hpp"

int main() {
    auto& nurbs_obj1 = dist_case1_obj1;

    int nuv = 1024;
    float ustep = (nurbs_obj1.get_uend() - nurbs_obj1.get_ubegin()) / nuv;
    float vstep = (nurbs_obj1.get_vend() - nurbs_obj1.get_vbegin()) / nuv;

    float *pResult = (float *)std::malloc(ONE_UV_SIZE * nuv * nuv * sizeof(float));

    auto start = std::chrono::steady_clock::now();
    evalWithCuda(nurbs_obj1.ivec, nurbs_obj1.data, nurbs_obj1.ndata, ustep, vstep, nurbs_obj1.get_ubegin(), nurbs_obj1.get_vbegin(), nuv, pResult);
    auto end = std::chrono::steady_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
    std::printf("opted time (CPU) : %.6f ms\n", elapsed_ms);

    std::free(pResult);

    printf("RUN Done\n");
    return 0;
}