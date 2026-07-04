# JSI_CAD_NEXT Version 5

Final paper implementation, based on Version 4.

- Deduplicates refinement tasks that reference the same AABB region.
- Evaluates each unique region once and shares its AABB result among all referencing pairs.
- Applies duplicate-region reuse to both distance and intersection refinement.
- Uses the scaled second-derivative `compute_K` padding for conservative
  intersection AABBs by default.
- Uses `nuv=128` for the first distance refinement round.

Common correctness update:

- Intersection uses `compute_K` second-derivative bounds to conservatively pad AABBs.
- Padding is enabled by default and can be disabled with `USE_COMPUTE_K=0`.
