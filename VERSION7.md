# JSI_CAD_NEXT Version 7

Based on Version 6.

Changes in this version:
- Deduplicates refinement tasks that refer to the same AABB region.
- Evaluates each unique region once and shares the generated AABB result among all pairs that reference it.
- Applies the reuse logic to both distance and intersection refinement.

Unchanged from Version 6:
- Distance pair collection uses OpenMP and SVE.
- Intersection pair collection uses slot-level OpenMP parallelism.
- AABB construction remains serial.
- Basis preprocessing is not included.
