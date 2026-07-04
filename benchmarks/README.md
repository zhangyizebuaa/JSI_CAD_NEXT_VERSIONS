# OpenCASCADE comparison

This benchmark converts the original JSI NURBS arrays directly into
`Geom_BSplineSurface` objects and measures OpenCASCADE query time separately
from surface construction.

On the Kunpeng server:

```bash
cd ~/JSI_CAD_NEXT_VERSIONS
bash benchmarks/build_opencascade_compare.sh
./build/opencascade/opencascade_compare distance 1 5
./build/opencascade/opencascade_compare intersection 1 5
```

The default paths use `JSI_CAD_NEXT_VERSION5` and the locally built
OpenCASCADE 8.0.0 tree. Override them with `VERSION_DIR` and `OCCT_BUILD`.

The distance benchmark runs OpenCASCADE's exact surface extrema query. The
intersection benchmark constructs exact intersection curves, while JSI's
reported intersection result is an AABB candidate count. Their timings
therefore describe different output contracts and should not be presented as
a strict kernel-for-kernel speedup.
