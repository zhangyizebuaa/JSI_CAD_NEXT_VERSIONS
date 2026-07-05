#include <iostream>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "jsi_cad.hpp"
#include "nurbs_data.hpp"
#include "common/common.cuh"

int main(int argc, char **argv) {
    int case_id = 1;
    if (argc > 1) {
        if (std::strcmp(argv[1], "1") == 0 || std::strcmp(argv[1], "case1") == 0) {
            case_id = 1;
        } else if (std::strcmp(argv[1], "2") == 0 || std::strcmp(argv[1], "case2") == 0) {
            case_id = 2;
        } else if (std::strcmp(argv[1], "3") == 0 || std::strcmp(argv[1], "case3") == 0) {
            case_id = 3;
        } else if (std::strcmp(argv[1], "4") == 0 || std::strcmp(argv[1], "case4") == 0) {
            case_id = 4;
        } else {
            std::cerr << "Usage: " << argv[0] << " [1|2|3|4|case1|case2|case3|case4]" << std::endl;
            return 2;
        }
    }

#ifdef _OPENMP
    if (std::getenv("OMP_NUM_THREADS") == nullptr) {
        omp_set_num_threads(51);
    }
    std::cout << "OpenMP threads: " << omp_get_max_threads() << std::endl;
#endif

    init_cad();

    NurbsFace *cases1[] = {&dist_case1_obj1, &dist_case2_obj1, &dist_hard_wave_obj1, &dist_hard_wave_obj1};
    NurbsFace *cases2[] = {&dist_case1_obj2, &dist_case2_obj2, &dist_hard_wave_obj2, &dist_hard_wave_obj2};
    auto& nurbs_obj1 = *cases1[case_id - 1];
    auto& nurbs_obj2 = *cases2[case_id - 1];
    std::cout << "Running distance case " << case_id << std::endl;

    std::vector<float> scaled_data_1;
    std::vector<float> scaled_data_2;
    const float *host_data_1 = nurbs_obj1.data;
    const float *host_data_2 = nurbs_obj2.data;
    if (case_id == 3) {
        constexpr float geometry_scale = 4.0f;
        scaled_data_1.assign(nurbs_obj1.data, nurbs_obj1.data + nurbs_obj1.ndata);
        scaled_data_2.assign(nurbs_obj2.data, nurbs_obj2.data + nurbs_obj2.ndata);
        const int cp_offset_1 = nurbs_obj1.ivec[3] + nurbs_obj1.ivec[6];
        const int cp_offset_2 = nurbs_obj2.ivec[3] + nurbs_obj2.ivec[6];
        const int xyz_count_1 = 3 * nurbs_obj1.ivec[2] * nurbs_obj1.ivec[5];
        const int xyz_count_2 = 3 * nurbs_obj2.ivec[2] * nurbs_obj2.ivec[5];
        for (int i = 0; i < xyz_count_1; ++i) scaled_data_1[cp_offset_1 + i] *= geometry_scale;
        for (int i = 0; i < xyz_count_2; ++i) scaled_data_2[cp_offset_2 + i] *= geometry_scale;
        host_data_1 = scaled_data_1.data();
        host_data_2 = scaled_data_2.data();
        std::cout << "Scaled hard case: geometry=4x, initial nuv=" << default_global_dist_nuv << std::endl;
    }

    int *d_ivec_1;
    float *d_data_1;
    int *d_ivec_2;
    float *d_data_2;
    
    d_ivec_1 = (int*)allocate_from_workspace(7 * sizeof(int));
    d_ivec_2 = (int*)allocate_from_workspace(7 * sizeof(int));
    d_data_1 = (float*)allocate_from_workspace(nurbs_obj1.ndata * sizeof(float));
    d_data_2 = (float*)allocate_from_workspace(nurbs_obj2.ndata * sizeof(float));

    cadMemcpy(d_ivec_1, nurbs_obj1.ivec, 7 * sizeof(int), CadMemcpyKind::HostToDevice);
    cadMemcpy(d_data_1, host_data_1, nurbs_obj1.ndata * sizeof(float), CadMemcpyKind::HostToDevice);
    cadMemcpy(d_ivec_2, nurbs_obj2.ivec, 7 * sizeof(int), CadMemcpyKind::HostToDevice);
    cadMemcpy(d_data_2, host_data_2, nurbs_obj2.ndata * sizeof(float), CadMemcpyKind::HostToDevice);
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
