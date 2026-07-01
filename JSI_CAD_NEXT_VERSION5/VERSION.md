# JSI_CAD_NEXT Version 5

Based on Version 4. This stage applies the change originally introduced in historical Version 7.

- Deduplicates refinement tasks that reference the same AABB region.
- Evaluates each unique region once and shares its AABB result among all referencing pairs.
- Applies duplicate-region reuse to both distance and intersection refinement.
