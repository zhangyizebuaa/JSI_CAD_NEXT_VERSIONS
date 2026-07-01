# 910B Hotspots Operator Library

This directory contains a buildable operator-library project that exports:

- `jsi_ascend_distance_topk_mask`
- `jsi_ascend_intersection_mask`

The main project already loads `libjsi_ascend_hotspots.so` dynamically and
calls these symbols when ACL mode is enabled.

## Current status

- ABI is final and compatible with the main project.
- Implementation is CPU/OpenMP fallback to keep behavior stable.
- File structure is ready for replacing internals with AscendC/ACLNN kernels.

## Build on 910B host

```bash
cd src/ascend/hotspots_910b
bash ./build_910b.sh
```

Output:

- `build/libjsi_ascend_hotspots.so`

Deploy this `.so` to runtime library path (or next to executable), e.g.:

```bash
export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH
```

## Next step (device execution)

Replace internals in `src/jsi_ascend_hotspots.cpp` with:

1. AscendC custom kernels (recommended for 910B), or
2. ACLNN operator chains + stream execution.

No changes are needed in the main project loader.
