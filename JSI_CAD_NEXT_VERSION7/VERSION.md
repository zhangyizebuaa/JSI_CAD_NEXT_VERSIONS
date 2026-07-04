# JSI_CAD_NEXT Version 7

Final consolidated source snapshot used by the artifact.

- Retains all OpenMP sampling and candidate-collection changes.
- Preserves selected pair correspondence across refinement rounds.
- Uses SVE for distance and intersection matrix screening.
- Batches independent intersection refinement tasks.
- Reuses duplicate distance and intersection subdomain evaluations.
- Applies the scaled second-derivative `compute_K` padding to conservative
  intersection AABBs by default.
- Supports four distance cases and two intersection cases.
- Uses 51 distance threads and 55 intersection threads by default; both can be
  overridden with `THREADS=<n>`.

Disable derivative-aware intersection padding only for a controlled comparison:

```bash
make run_intersect CASE=1 USE_COMPUTE_K=0
```
