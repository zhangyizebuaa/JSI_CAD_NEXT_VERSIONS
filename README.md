# JSI_CAD_NEXT Versions

This repository stores reproducible JSI_CAD_NEXT source-code versions in separate directories.

| Directory | Optimization stage |
| --- | --- |
| `JSI_CAD_NEXT_VERSION0/` | Serial CPU baseline |
| `JSI_CAD_NEXT_VERSION1/` | Historical Versions 1, 3, and 6: grouped OpenMP changes |
| `JSI_CAD_NEXT_VERSION2/` | Historical Version 2: preserve candidate-pair correspondence |
| `JSI_CAD_NEXT_VERSION3/` | Historical Version 4: add SVE |
| `JSI_CAD_NEXT_VERSION4/` | Historical Version 5: batch intersection refinement tasks |
| `JSI_CAD_NEXT_VERSION5/` | Final version: reuse duplicate refinement regions |

Build outputs, IDE files, and performance reports are not tracked.

Run a selected case inside any version directory:

```bash
make run_dist CASE=1
make run_dist CASE=3
make run_intersect CASE=2
make run_all_cases
```

Distance cases 1 and 2 are the original inputs. Case 3 is a 4x geometric scaling of the hard-wave
input, paired with the 4x denser initial grid to retain approximately the same physical sampling
spacing. Case 4 is the unscaled hard-wave input. Intersection currently provides cases 1 and 2.
`run_all_cases` runs all four distance cases followed by both intersection cases.

Distance uses 51 OpenMP threads by default and intersection uses 55. Override either command with
`THREADS=<n>`, or set separate defaults with `DIST_THREADS=<n>` and `INTER_THREADS=<n>`.
Distance uses an initial sampling width of `nuv=512` in every version; later refinement rounds
retain the existing `nuv=16` setting.

All versions now include the same second-derivative-aware intersection AABB padding so that
cross-version performance comparisons use the same conservative bound. It is enabled by default
and can be disabled for comparison:

```bash
make run_intersect CASE=1 USE_COMPUTE_K=1
make run_intersect CASE=1 USE_COMPUTE_K=0
```

For paper reproduction, use `JSI_CAD_NEXT_VERSION5/` as the final implementation.
