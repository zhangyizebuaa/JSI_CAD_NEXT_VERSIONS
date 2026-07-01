#ifndef _EVALUATION_CUH
#define _EVALUATION_CUH

void
EG_spline2dDeriv_irrational(int *ivec,
                 float *data,
                 float ustep,
                 float vstep,
                 float ustart,
                 float vstart,
                 int nuv,
                 float *results);

void
EG_spline2dDeriv_irrational_opt(int *ivec,
                 float *data,
                 float ustep,
                 float vstep,
                 float ustart,
                 float vstart,
                 int nuv,
                 int ystart,
                 int yend,
                 float *results);

#endif