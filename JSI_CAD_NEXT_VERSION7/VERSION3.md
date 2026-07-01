# JSI_CAD_NEXT Version 3

Based on Version 2.

Changes in this version:
- Adds OpenMP to distance pair collection.
- Each OpenMP thread keeps a local top-k pair list, then local lists are merged into the final top-k pair list.

Unchanged from Version 2:
- Pair refine remains enabled for distance and intersection.
- Sampling OpenMP behavior remains the same as Version 1 and Version 2.
- AABB construction remains serial.
- Intersection pair collection remains serial.
- No SVE path is added.
