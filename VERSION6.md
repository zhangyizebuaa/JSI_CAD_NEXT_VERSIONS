# JSI_CAD_NEXT Version 6

Based on Version 5.

Changes in this version:
- Adds slot-level OpenMP parallelism to intersection pair collection.

Implementation notes:
- `intersect_collect_pairs` parallelizes over aligned slots when there are more than 4 slots.
- Each slot writes overlap pairs into its own local buffer.
- Local buffers are merged back in slot order to keep deterministic pair ordering.
- Next-round intersection refine sampling remains OpenMP-parallel from Version 5.
- Distance logic is unchanged from Version 5.

Unchanged from Version 5:
- Matrix-internal pair collection is not parallelized.
- AABB construction remains serial.
- Basis preprocessing is not included in this version.
