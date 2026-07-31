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
spacing. Case 4 is the unscaled hard-wave input. Intersection provides cases 1 through 4.
`run_all_cases` runs all four distance cases followed by all four intersection cases.

Distance uses 51 OpenMP threads for the current paper measurements and intersection uses 55.
Override either command with
`THREADS=<n>`, or set separate defaults with `DIST_THREADS=<n>` and `INTER_THREADS=<n>`.
Distance uses an initial sampling width of `nuv=512` in every version; later refinement rounds
retain the existing `nuv=16` setting.
The historical serial Version 0 materializes the complete first-round distance matrix and therefore
cannot execute the `nuv=512` workload inside the 4 GiB workspace. It remains a source-level
reference; the paper's measured distance ablation starts from Version 1, whose bounded top-k
collection avoids that allocation.

All versions now include the same second-derivative-aware intersection AABB padding so that
cross-version performance comparisons use the same conservative bound. It is enabled by default
and can be disabled for comparison:

```bash
make run_intersect CASE=1 USE_COMPUTE_K=1
make run_intersect CASE=1 USE_COMPUTE_K=0
```

For paper reproduction, use `JSI_CAD_NEXT_VERSION5/` as the final implementation.

## Artifact layout and reproduction

Source-code versions:

- `JSI_CAD_NEXT_VERSION0/`: serial CPU baseline source.
- `JSI_CAD_NEXT_VERSION1/`: OpenMP baseline source.
- `JSI_CAD_NEXT_VERSION2/`: pair-preserving distance/intersection source.
- `JSI_CAD_NEXT_VERSION3/`: SVE-enabled matrix screening source.
- `JSI_CAD_NEXT_VERSION4/`: batched intersection-refinement source.
- `JSI_CAD_NEXT_VERSION5/`: final source with duplicate subregion reuse.

Build a selected version:

```bash
cd JSI_CAD_NEXT_VERSION5
make clean
make all
```

Run one JSI case:

```bash
OMP_PROC_BIND=close OMP_PLACES=cores make run_dist CASE=1 THREADS=51
OMP_PROC_BIND=close OMP_PLACES=cores make run_intersect CASE=1 THREADS=55
```

Run all source versions, four distance cases, four intersection cases, and the
OpenCASCADE comparison. Each measurement is repeated five times and the best
time is reported:

```bash
OMP_PROC_BIND=close OMP_PLACES=cores \
  bash benchmarks/run_all_benchmarks.sh
```

Useful environment overrides:

```bash
REPEATS=5
VERSIONS="1 2 3 4 5"
DISTANCE_CASES="1 2 3"
INTERSECTION_CASES="1 2 3"
DISTANCE_THREADS=264
DISTANCE_CPUSET="0-263"
INTERSECTION_THREADS=56
INTERSECTION_CPUSET="0-55"
RUN_OCCT=1
OCCT_BUILD=/path/to/OCCT-8.0.0-build
```

Test data are stored in every version under `examples/data/`. Distance uses
`dist_case1.cu` through `dist_case3.cu`; the current Case 3 is the scaled
hard-wave workload previously kept as the fifth experimental case. Intersection uses
`inter_case*-1_data.cu` and `inter_case*-2_data.cu`.

Result validation uses the printed query result in each log:

- Distance: compare the final `[current_min_dist]` value across versions for
  the same case.
- Intersection: compare the final `overlap pairs` count across versions for
  the same case.
- OpenCASCADE: compare `distance=` or `intersection_curves=` in the OCCT log.

The script stores full logs plus `raw_runs.csv` and `mean_summary.csv` under
`results/<timestamp>/`.

Set `CAD_TIMING_BREAKDOWN=1` to print both timing scopes used by the paper:

```bash
CAD_TIMING_BREAKDOWN=1 OMP_PROC_BIND=close OMP_PLACES=cores \
  make -C JSI_CAD_NEXT_VERSION5 run_dist CASE=1 THREADS=51
CAD_TIMING_BREAKDOWN=1 OMP_PROC_BIND=close OMP_PLACES=cores \
  make -C JSI_CAD_NEXT_VERSION5 run_intersect CASE=1 THREADS=55
```

`end_to_end_ms` covers the complete query. `kernel_ms` is the accumulated
wall-clock time for surface evaluation/AABB construction and AABB matrix
screening. The remaining `orchestration_ms` covers pair metadata, deduplication,
and result reassembly.
