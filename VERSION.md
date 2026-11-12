# Version

Submission reproducibility package (de-MR)

Prepared for the final TP53I3 manuscript and supplementary materials.

Date: 2026-11-13

DOI: pending Zenodo assignment

## Release v4.0 de-MR

- Removes the former exploratory genetic-instrument and regional-overlap
  analyses, including their scripts, result tables, figure sources, software
  entries, and data-source rows.
- Renumbers Supplementary Figures to S1-S13, Tables to S1-S17, and Notes to
  A1-A8 to match the final manuscript.
- Regenerates checksums and the release archive.

## v3.1 statistical consistency revision

- Changes the F0-F4 scale-sensitivity variance from `1/(n - 3)` to the
  manuscript-specified `1.06/(n - 3)`.
- Regenerates the F0-F4 pooled estimate, confidence and prediction intervals,
  cohort weights, and leave-one-out results.
- Routes outputs to `results/` through `TP53I3_RESULT_DIR`.
- Corrects the GSE256398 data-source metadata from bulk to single-nucleus
  transcriptome.
- Updates the analysis map and package checksums.
- Preserves the original v3 package separately for traceability.
