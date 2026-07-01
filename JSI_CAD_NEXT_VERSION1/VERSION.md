# JSI_CAD_NEXT Version 1

Based on Version 0. This stage combines the changes originally introduced in historical Versions 1, 3, and 6.

- OpenMP-parallel NURBS sampling from historical Version 1.
- OpenMP distance candidate collection with thread-local top-k lists from historical Version 3.
- OpenMP intersection candidate collection across slots from historical Version 6.
- Pair/slot support structures are included as dependencies.
- Candidate collection scans the Cartesian product of the independently retained side1 and side2 regions.
- Each retained region is refined and stored once; pair correspondence is not preserved yet.
- SVE, batched intersection refinement, and duplicate-region reuse are not included.
