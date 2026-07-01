# JSI_CAD_NEXT Version 1

CPU implementation with OpenMP parallelization only in the surface sampling stage.

This version is intended to be compared against Version 0, the single-thread CPU baseline.

## Optimization Scope

Version 1 adds OpenMP only to the first compute-heavy step in the pipeline:

```text
NURBS surface input
  -> surface sampling on the u/v grid  [OpenMP added here]
  -> AABB construction                 [serial]
  -> distance/intersection matrix scan [serial]
  -> mask/gather/refine control flow   [serial]
```

The goal of this version is to measure the speedup from parallelizing only the independent `u/v` sampling points.

## OpenMP Build Switches

- `Makefile`
  - `CXXFLAGS` adds `-fopenmp`.
  - `PERF_CXXFLAGS` adds `-fopenmp`.

## OpenMP Parallel Region

- `src/evaluation/evaluation.cu`
  - `EG_spline2dDeriv_irrational_threshold`
    - Parallelizes the `nuv * nuv` surface sampling loop.
    - Each OpenMP iteration evaluates one sampled `(u, v)` point and writes its `xyz` and derivative data.
    - Uses `#pragma omp parallel for schedule(static) collapse(2) if (nuv * nuv > omp_threshold)`.
  - `EG_spline2dDeriv_irrational`
    - Calls `EG_spline2dDeriv_irrational_threshold(..., omp_threshold = 0)`.
    - This means the distance path opens OpenMP by default for sampling.
  - `EG_spline2dDeriv_irrational_opt`
    - Strip-based sampling helper used by the evaluation example.
    - Also parallelizes sampling only.

## Distance Path

- `src/distance_analysis/distance_analysis.cu`
  - `run_evaluation_and_construction_task`
    - Calls `EG_spline2dDeriv_irrational(...)`.
    - Therefore distance sampling opens OpenMP by default.
  - `construct_aabbs_distance_host`
    - Serial in this version.
  - `check_distance_of_aabbs`
    - Distance matrix computation and mask marking are serial in this version.

## Intersection Path

- `src/intersection/intersection.cu`
  - `run_evaluation_and_construction_task`
    - Calls `EG_spline2dDeriv_irrational_threshold(..., omp_threshold = 512)`.
    - Therefore intersection sampling opens OpenMP only when `nuv * nuv > 512`.
  - `construct_aabbs_host`
    - Serial in this version.
  - `para_check_aabbs_host`
    - AABB overlap scanning and mask marking are serial in this version.

## Not Parallelized In This Version

- AABB construction.
- Distance matrix computation.
- Distance mask marking.
- Intersection AABB overlap scanning.
- Intersection mask marking.
- Candidate queue management.
- `find_aabb`.
- `gather_idx_kernel_host`.
- Round/refine control flow.
