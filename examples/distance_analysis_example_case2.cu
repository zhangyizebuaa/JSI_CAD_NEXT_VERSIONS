#include <iostream>
#include <chrono>
#include <cstring>

#include "jsi_cad.hpp"
#include "nurbs_data.hpp"
#include "common/common.cuh"

int main() {
    init_cad();

    auto& nurbs_obj1 = dist_case2_obj1;
    auto& nurbs_obj2 = dist_case2_obj2;

    int *d_ivec_1;
    float *d_data_1;
    int *d_ivec_2;
    float *d_data_2;
    
    d_ivec_1 = (int*)allocate_from_workspace(7 * sizeof(int));
    d_ivec_2 = (int*)allocate_from_workspace(7 * sizeof(int));
    d_data_1 = (float*)allocate_from_workspace(nurbs_obj1.ndata * sizeof(float));
    d_data_2 = (float*)allocate_from_workspace(nurbs_obj2.ndata * sizeof(float));

    cadMemcpy(d_ivec_1, nurbs_obj1.ivec, 7 * sizeof(int), CadMemcpyKind::HostToDevice);
    cadMemcpy(d_data_1, nurbs_obj1.data, nurbs_obj1.ndata * sizeof(float), CadMemcpyKind::HostToDevice);
    cadMemcpy(d_ivec_2, nurbs_obj2.ivec, 7 * sizeof(int), CadMemcpyKind::HostToDevice);
    cadMemcpy(d_data_2, nurbs_obj2.data, nurbs_obj2.ndata * sizeof(float), CadMemcpyKind::HostToDevice);
    EvalAndConstructTask t1, t2;
    t1.d_ivec = d_ivec_1;
    t1.d_data = d_data_1;
    t1.d_ndata = nurbs_obj1.ndata;

    // t1.ustart = nurbs_obj1.get_ubegin();
    // t1.vstart = nurbs_obj1.get_vbegin();
    // t1.ustop = nurbs_obj1.get_ubegin() + (nurbs_obj1.get_uend() - nurbs_obj1.get_ubegin()) / 2;
    // t1.vstop = nurbs_obj1.get_vbegin() + (nurbs_obj1.get_vend() - nurbs_obj1.get_vbegin()) / 2;

    // t1.ustart = nurbs_obj1.get_ubegin() + (nurbs_obj1.get_uend() - nurbs_obj1.get_ubegin()) / 2;
    // t1.vstart = nurbs_obj1.get_vbegin();
    // t1.ustop = nurbs_obj1.get_uend();
    // t1.vstop = nurbs_obj1.get_vbegin() + (nurbs_obj1.get_vend() - nurbs_obj1.get_vbegin()) / 2;

    // t1.ustart = nurbs_obj1.get_ubegin();
    // t1.vstart = nurbs_obj1.get_vbegin() + (nurbs_obj1.get_vend() - nurbs_obj1.get_vbegin()) / 2;
    // t1.ustop = nurbs_obj1.get_ubegin() + (nurbs_obj1.get_uend() - nurbs_obj1.get_ubegin()) / 2;
    // t1.vstop = nurbs_obj1.get_vend();

    // t1.ustart = nurbs_obj1.get_ubegin() + (nurbs_obj1.get_uend() - nurbs_obj1.get_ubegin()) / 2;
    // t1.vstart = nurbs_obj1.get_vbegin() + (nurbs_obj1.get_vend() - nurbs_obj1.get_vbegin()) / 2;
    // t1.ustop = nurbs_obj1.get_uend();
    // t1.vstop = nurbs_obj1.get_vend();

    t1.ustart = nurbs_obj1.get_ubegin();
    t1.vstart = nurbs_obj1.get_vbegin();
    t1.ustop = nurbs_obj1.get_uend();
    t1.vstop = nurbs_obj1.get_vend();

    t1.nuv = default_global_dist_nuv;
    t1.eval_task_side = 1;

    std::cout << "face1 u:" << t1.ustart << " - " << t1.ustop << std::endl;
    std::cout << "face1 v:" << t1.vstart << " - " << t1.vstop << std::endl;

    t2.d_ivec = d_ivec_2;
    t2.d_data = d_data_2;
    t2.d_ndata = nurbs_obj2.ndata;

    // t2.ustart = nurbs_obj2.get_ubegin();
    // t2.vstart = nurbs_obj2.get_vbegin();
    // t2.ustop = nurbs_obj2.get_ubegin() + (nurbs_obj2.get_uend() - nurbs_obj2.get_ubegin()) / 2;
    // t2.vstop = nurbs_obj2.get_vbegin() + (nurbs_obj2.get_vend() - nurbs_obj2.get_vbegin()) / 2;
    
    // t2.ustart = nurbs_obj2.get_ubegin() + (nurbs_obj2.get_uend() - nurbs_obj2.get_ubegin()) / 2;
    // t2.vstart = nurbs_obj2.get_vbegin();
    // t2.ustop = nurbs_obj2.get_uend();
    // t2.vstop = nurbs_obj2.get_vbegin() + (nurbs_obj2.get_vend() - nurbs_obj2.get_vbegin()) / 2;
    
    // t2.ustart = nurbs_obj2.get_ubegin();
    // t2.vstart = nurbs_obj2.get_vbegin() + (nurbs_obj2.get_vend() - nurbs_obj2.get_vbegin()) / 2;
    // t2.ustop = nurbs_obj2.get_ubegin() + (nurbs_obj2.get_uend() - nurbs_obj2.get_ubegin()) / 2;
    // t2.vstop = nurbs_obj2.get_vend();

    // t2.ustart = nurbs_obj2.get_ubegin() + (nurbs_obj2.get_uend() - nurbs_obj2.get_ubegin()) / 2;
    // t2.vstart = nurbs_obj2.get_vbegin() + (nurbs_obj2.get_vend() - nurbs_obj2.get_vbegin()) / 2;
    // t2.ustop = nurbs_obj2.get_uend();
    // t2.vstop = nurbs_obj2.get_vend();

    t2.ustart = nurbs_obj2.get_ubegin();
    t2.vstart = nurbs_obj2.get_vbegin();
    t2.ustop = nurbs_obj2.get_uend();
    t2.vstop = nurbs_obj2.get_vend();
    
    t2.nuv = default_global_dist_nuv;
    t2.eval_task_side = 2;

    std::cout << "face2 u:" << t2.ustart << " - " << t2.ustop << std::endl;
    std::cout << "face2 v:" << t2.vstart << " - " << t2.vstop << std::endl;

    // warmup
    // run_evaluation_and_construction_task(t1, 0);
    minimum_distance(t1, t2);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    minimum_distance(t1, t2);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    printf("Check took %lfms\n", std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0);

    printf("RUN Done\n");
    deinit_cad();
    return 0;
}