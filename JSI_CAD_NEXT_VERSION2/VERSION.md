# JSI_CAD_NEXT Version 2

Based on Version 1. This stage applies the change originally introduced in historical Version 2.

- Preserves each selected `(side1 AABB, side2 AABB)` pair when constructing the next refinement round.
- Stops recombining independently selected side1 and side2 endpoints as a Cartesian product.
- Retains all OpenMP changes from Version 1.
- SVE, batched intersection refinement, and duplicate-region reuse are not included.

Common correctness update:

- Intersection uses `compute_K` second-derivative bounds to conservatively pad AABBs.
- Padding is enabled by default and can be disabled with `USE_COMPUTE_K=0`.
