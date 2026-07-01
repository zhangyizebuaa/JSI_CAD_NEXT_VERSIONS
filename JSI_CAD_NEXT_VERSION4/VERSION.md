# JSI_CAD_NEXT Version 4

Based on Version 3. This stage applies the change originally introduced in historical Version 5.

- Builds next-round intersection refinement tasks first.
- Executes the independent refinement tasks as an OpenMP batch.
- Retains pair correspondence, SVE, and the earlier OpenMP changes.
- Duplicate-region reuse is not included.
