#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}" && pwd)"
BUILD_DIR="${ROOT_DIR}/build"

ASCEND_HOME="${ASCEND_HOME:-/usr/local/Ascend/ascend-toolkit/latest}"
CXX="${CXX:-g++}"

mkdir -p "${BUILD_DIR}"

${CXX} \
  -O3 -fPIC -shared -std=c++17 \
  -march=armv8.2-a \
  -I"${ROOT_DIR}/include" \
  -I"${ASCEND_HOME}/include" \
  "${ROOT_DIR}/src/jsi_ascend_hotspots.cpp" \
  -L"${ASCEND_HOME}/lib64" -lacl \
  -o "${BUILD_DIR}/libjsi_ascend_hotspots.so"

echo "Built: ${BUILD_DIR}/libjsi_ascend_hotspots.so"
