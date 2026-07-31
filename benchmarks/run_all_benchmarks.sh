#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
repeats="${REPEATS:-5}"
distance_threads="${DISTANCE_THREADS:-264}"
intersection_threads="${INTERSECTION_THREADS:-56}"
distance_cpuset="${DISTANCE_CPUSET:-0-263}"
intersection_cpuset="${INTERSECTION_CPUSET:-0-55}"
versions="${VERSIONS:-0 1 2 3 4 5}"
distance_cases="${DISTANCE_CASES:-1 2 3}"
intersection_cases="${INTERSECTION_CASES:-1 2 3}"
run_occt="${RUN_OCCT:-1}"
occt_threads="${OCCT_THREADS:-1}"
timestamp="$(date +%Y%m%d_%H%M%S)"
result_dir="${RESULT_DIR:-$repo_dir/results/$timestamp}"
mkdir -p "$result_dir/logs"

raw_csv="$result_dir/raw_runs.csv"
summary_csv="$result_dir/mean_summary.csv"
echo "engine,version,mode,case,repeat,status,elapsed_ms,reported_ms,result,log" > "$raw_csv"

extract_reported_ms() {
  awk '
    /Check took/ {
      for (i = 1; i <= NF; ++i) {
        if ($i ~ /^[0-9.]+ms$/) {
          gsub("ms", "", $i);
          value = $i;
        }
      }
    }
    /best_ms=/ {
      split($0, a, "=");
      value = a[2];
    }
    END {
      if (value == "") value = "NA";
      print value;
    }
  ' "$1"
}

extract_result() {
  awk '
    /\[current_min_dist\]/ { result = $0 }
    /overlap pairs/ { result = $0 }
    /distance=/ { result = $0 }
    /intersection_curves=/ { result = $0 }
    END {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", result);
      gsub(/,/, ";", result);
      if (result == "") result = "NA";
      print result;
    }
  ' "$1"
}

run_one() {
  local engine="$1"
  local version="$2"
  local mode="$3"
  local case_id="$4"
  local repeat_id="$5"
  shift 5
  local log="$result_dir/logs/${engine}_${version}_${mode}_case${case_id}_run${repeat_id}.log"
  local start_ns end_ns elapsed_ms status reported result
  start_ns="$(date +%s%N)"
  set +e
  "$@" > "$log" 2>&1
  status="$?"
  set -e
  end_ns="$(date +%s%N)"
  elapsed_ms="$(awk -v s="$start_ns" -v e="$end_ns" 'BEGIN { printf "%.3f", (e - s) / 1000000.0 }')"
  reported="$(extract_reported_ms "$log")"
  result="$(extract_result "$log")"
  echo "$engine,$version,$mode,$case_id,$repeat_id,$status,$elapsed_ms,$reported,$result,$log" >> "$raw_csv"
  if [[ "$status" != "0" ]]; then
    echo "FAILED engine=$engine version=$version mode=$mode case=$case_id repeat=$repeat_id; see $log"
  fi
}

for version in $versions; do
  version_dir="$repo_dir/JSI_CAD_NEXT_VERSION${version}"
  if [[ ! -f "$version_dir/Makefile" ]]; then
    echo "skip missing version directory: $version_dir"
    continue
  fi

  echo "== build JSI version $version =="
  make -C "$version_dir" build/test_dist build/test_intersect

  for case_id in $distance_cases; do
    for ((run = 1; run <= repeats; ++run)); do
      run_one jsi "v$version" distance "$case_id" "$run" \
        taskset -c "$distance_cpuset" env CAD_TIMING_BREAKDOWN=1 OMP_PROC_BIND=close \
        OMP_PLACES=cores OMP_NUM_THREADS="$distance_threads" \
        "$version_dir/build/test_dist" "$case_id"
    done
  done

  for case_id in $intersection_cases; do
    for ((run = 1; run <= repeats; ++run)); do
      run_one jsi "v$version" intersection "$case_id" "$run" \
        taskset -c "$intersection_cpuset" env CAD_TIMING_BREAKDOWN=1 CAD_INTERSECT_USE_K=1 \
        OMP_PROC_BIND=close OMP_PLACES=cores \
        OMP_NUM_THREADS="$intersection_threads" "$version_dir/build/test_intersect" "$case_id"
    done
  done
done

if [[ "$run_occt" == "1" ]]; then
  echo "== build OpenCASCADE comparison =="
  occt_bin="$("$repo_dir/benchmarks/build_opencascade_compare.sh")"
  for case_id in $distance_cases; do
    run_one occt final distance "$case_id" 1 "$occt_bin" distance "$case_id" "$repeats" "$occt_threads"
  done
  for case_id in $intersection_cases; do
    run_one occt final intersection "$case_id" 1 "$occt_bin" intersection "$case_id" "$repeats" "$occt_threads"
  done
fi

awk -F, '
  NR == 1 { next }
  $6 != 0 { next }
  {
    key = $1 "," $2 "," $3 "," $4;
    t = ($8 == "NA" ? $7 : $8) + 0.0;
    sum[key] += t;
    count[key] += 1;
    if (!(key in minimum) || t < minimum[key]) minimum[key] = t;
    if (!(key in maximum) || t > maximum[key]) maximum[key] = t;
    if (!(key in result)) {
      result[key] = $9;
      stable[key] = "yes";
    } else if (result[key] != $9) {
      stable[key] = "no";
    }
  }
  END {
    print "engine,version,mode,case,runs,mean_ms,min_ms,max_ms,result,result_stable";
    for (key in sum) {
      print key "," count[key] "," sum[key] / count[key] "," minimum[key] "," maximum[key] "," result[key] "," stable[key];
    }
  }
' "$raw_csv" > "$summary_csv"

echo "raw_runs=$raw_csv"
echo "mean_summary=$summary_csv"
