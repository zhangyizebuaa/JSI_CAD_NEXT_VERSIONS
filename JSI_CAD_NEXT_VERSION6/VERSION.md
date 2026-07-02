# JSI_CAD_NEXT Version 6

Based on Version 5. This is an experimental second-derivative-aware intersection variant.

- Calls `compute_K` after evaluating each intersection surface patch.
- Scales the bound by the squared relative parameter-domain size so padding shrinks during refinement.
- Uses the scaled second-derivative bound as conservative padding on every AABB radius.
- Keeps the original AABB overlap predicate; only the box extent changes.
- Enable or disable the experiment with `USE_COMPUTE_K=1` or `USE_COMPUTE_K=0`.
