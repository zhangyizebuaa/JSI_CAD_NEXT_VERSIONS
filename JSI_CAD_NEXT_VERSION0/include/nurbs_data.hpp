struct NurbsFace {
    int* ivec;
    int ndata;
    float* data;

    float get_ubegin() {return data[0];}
    float get_uend() {return data[ivec[3] - 1];}
    float get_vbegin() {return data[ivec[3]];}
    float get_vend() {return data[ivec[3] + ivec[6] - 1];}
};

extern struct NurbsFace dist_case1_obj1;
extern struct NurbsFace dist_case1_obj2;
extern struct NurbsFace dist_case2_obj1;
extern struct NurbsFace dist_case2_obj2;
extern struct NurbsFace dist_case_hard_parallel_obj1;
extern struct NurbsFace dist_case_hard_parallel_obj2;
extern struct NurbsFace dist_hard_wave_obj1;
extern struct NurbsFace dist_hard_wave_obj2;

extern struct NurbsFace inter_case1_obj1;
extern struct NurbsFace inter_case1_obj2;
extern struct NurbsFace inter_case2_obj1;
extern struct NurbsFace inter_case2_obj2;
extern struct NurbsFace inter_case3_obj1;
extern struct NurbsFace inter_case3_obj2;
extern struct NurbsFace inter_case4_obj1;
extern struct NurbsFace inter_case4_obj2;
