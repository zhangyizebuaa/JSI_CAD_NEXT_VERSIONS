# JSI_CAD_NEXT Version 3

Based on Version 2.

- Deduplicates refinement tasks that refer to the same AABB region.
- Evaluates each unique region once and shares the generated AABB result among all pairs that reference it.
- Applies the reuse logic to both distance and intersection refinement.
- Retains the OpenMP and SVE optimizations from the previous versions.
