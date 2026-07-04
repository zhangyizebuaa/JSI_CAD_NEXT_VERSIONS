#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
version_dir="${VERSION_DIR:-$repo_dir/JSI_CAD_NEXT_VERSION5}"
occt_build="${OCCT_BUILD:-/home/zyz/src/OCCT-8.0.0-build}"
output_dir="$repo_dir/build/opencascade"

make -C "$version_dir" data
mkdir -p "$output_dir"

g++ -O3 -std=c++17 \
  -I"$version_dir/include" \
  -I"$occt_build/include/opencascade" \
  "$repo_dir/benchmarks/opencascade_compare.cpp" \
  "$version_dir/build/data_i11.o" \
  "$version_dir/build/data_i12.o" \
  "$version_dir/build/data_i21.o" \
  "$version_dir/build/data_i22.o" \
  "$version_dir/build/data_d1.o" \
  "$version_dir/build/data_d2.o" \
  "$version_dir/build/data_d3.o" \
  "$version_dir/build/data_d4.o" \
  -L"$occt_build/lin64/gcc/lib" \
  -Wl,-rpath,"$occt_build/lin64/gcc/lib" \
  -lTKBO -lTKBool -lTKGeomAlgo -lTKTopAlgo -lTKBRep -lTKGeomBase \
  -lTKG3d -lTKG2d -lTKMath -lTKernel \
  -o "$output_dir/opencascade_compare"

echo "$output_dir/opencascade_compare"
