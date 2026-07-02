# JSI_CAD_NEXT Version 3

Based on Version 2. This stage applies the change originally introduced in historical Version 4.

- Adds SVE vector distance evaluation during distance top-k collection.
- Adds SVE vector AABB overlap evaluation during intersection candidate collection.
- Retains pair correspondence and all earlier OpenMP changes.
- Batched intersection refinement and duplicate-region reuse are not included.

Common correctness update:

- Intersection uses `compute_K` second-derivative bounds to conservatively pad AABBs.
- Padding is enabled by default and can be disabled with `USE_COMPUTE_K=0`.
