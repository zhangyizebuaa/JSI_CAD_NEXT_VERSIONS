# JSI_CAD_NEXT Versions

This repository stores reproducible JSI_CAD_NEXT source-code versions in separate directories.

| Directory | Optimization stage |
| --- | --- |
| `JSI_CAD_NEXT_VERSION0/` | Serial CPU baseline |
| `JSI_CAD_NEXT_VERSION1/` | Historical Versions 1, 3, and 6: grouped OpenMP changes |
| `JSI_CAD_NEXT_VERSION2/` | Historical Version 2: preserve candidate-pair correspondence |
| `JSI_CAD_NEXT_VERSION3/` | Historical Version 4: add SVE |
| `JSI_CAD_NEXT_VERSION4/` | Historical Version 5: batch intersection refinement tasks |
| `JSI_CAD_NEXT_VERSION5/` | Historical Version 7: reuse duplicate refinement regions |
| `JSI_CAD_NEXT_VERSION6/` | Original snapshot of the Compute_K intersection experiment |

Build outputs, IDE files, and performance reports are not tracked.

Run a selected case inside any version directory:

```bash
make run_dist CASE=1
make run_dist CASE=3
make run_intersect CASE=2
make run_all_cases
```

Distance cases 1 and 2 are the original inputs. Case 3 is the hard-parallel input and case 4 is the
hard-wave input. Intersection currently provides cases 1 and 2. `run_all_cases` runs all four
distance cases followed by both intersection cases.

Distance uses 51 OpenMP threads by default and intersection uses 55. Override either command with
`THREADS=<n>`, or set separate defaults with `DIST_THREADS=<n>` and `INTER_THREADS=<n>`.

All versions now include the same second-derivative-aware intersection AABB padding so that
cross-version performance comparisons use the same conservative bound. It is enabled by default
and can be disabled for comparison:

```bash
make run_intersect CASE=1 USE_COMPUTE_K=1
make run_intersect CASE=1 USE_COMPUTE_K=0
```
