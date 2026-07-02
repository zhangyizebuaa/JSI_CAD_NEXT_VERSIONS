# JSI_CAD_NEXT Version 0

Serial CPU baseline used for speedup comparisons.

- OpenMP compile flags and pragmas are disabled.
- SVE is not used.
- Distance and intersection use the original serial refinement flow.

Common correctness update:

- Intersection uses `compute_K` second-derivative bounds to conservatively pad AABBs.
- Padding is enabled by default and can be disabled with `USE_COMPUTE_K=0`.
