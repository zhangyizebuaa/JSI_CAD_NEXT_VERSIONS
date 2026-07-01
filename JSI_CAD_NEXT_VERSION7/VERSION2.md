# JSI_CAD_NEXT Version 2

Based on Version 1.

Changes in this version:
- Keep Version 1 sampling OpenMP behavior unchanged.
- Distance refine now uses paired AABB candidates: each selected `(side1 box, side2 box)` is refined as one aligned slot.
- Intersection refine now uses paired AABB candidates: each overlapping `(side1 box, side2 box)` is refined as one aligned slot.

Unchanged from Version 1:
- No SVE path is added.
- No OpenMP is added to AABB construction, distance matrix scan, or intersection matrix scan.
- Distance sampling still uses OpenMP by default.
- Intersection sampling still enables OpenMP only when `nuv * nuv > 512`.
