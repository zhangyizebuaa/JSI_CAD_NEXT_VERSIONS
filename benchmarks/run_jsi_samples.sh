#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
version_dir="${VERSION_DIR:-$repo_dir/JSI_CAD_NEXT_VERSION5}"
repeats="${REPEATS:-5}"
distance_case="${DISTANCE_CASE:-1}"
intersection_case="${INTERSECTION_CASE:-1}"
distance_threads="${DISTANCE_THREADS:-51}"
intersection_threads="${INTERSECTION_THREADS:-55}"

make -C "$version_dir" build/test_dist build/test_intersect

echo "load_before=$(uptime)"
for ((i = 1; i <= repeats; ++i)); do
  echo "jsi_distance_run=$i"
  OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS="$distance_threads" \
    "$version_dir/build/test_dist" "$distance_case"
done

for ((i = 1; i <= repeats; ++i)); do
  echo "jsi_intersection_run=$i"
  CAD_INTERSECT_USE_K=1 OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS="$intersection_threads" \
    "$version_dir/build/test_intersect" "$intersection_case"
done
echo "load_after=$(uptime)"
