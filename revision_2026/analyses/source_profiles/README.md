# Source-profile PERMANOVA and dispersion diagnostic

The existing R script and two input tables are included here. This module was not rerun during
preparation. Package versions from the saved analysis are in `data/dispersion_software_versions.tsv`.

Run from `revision_2026/` after installing the recorded R/vegan/permute versions:

```bash
Rscript analyses/source_profiles/check_dispersion.R analyses/source_profiles
```

The script writes diagnostic result tables under this module's `data/` directory.
It uses within-species proportions, Bray-Curtis dissimilarity, 999 PERMANOVA permutations,
and 9,999 dispersion permutations. The saved manuscript results remain separately archived in
`source_data/tables/S3/`. A completed numerical reproduction is not claimed for this module here.
