# JSI_CAD_NEXT Version 1

Based on Version 0. This version groups all historical OpenMP optimizations into one stage.

OpenMP is used for:

- NURBS surface sampling on the u/v grid.
- Distance top-k pair collection with thread-local candidate lists and deterministic merging.
- Intersection pair collection across aligned slots.
- Batched construction of next-round intersection refinement tasks.

The aligned pair-refinement data structures required by the later OpenMP regions are included in this version.
SVE and duplicate refinement-task reuse are not included.
