[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23940279.svg)](https://doi.org/10.5281/zenodo.23940279)
# TP53I3 manuscript reproducibility package

This package contains the analysis code, derived result tables, software
versions, random seeds, and frozen plotting inputs needed to verify the reported
results of the TP53I3 manuscript.

## Contents

- `data_sources.csv`: public data resources and accession identifiers.
- `analysis_map.csv`: main-text results mapped to code and result files.
- `supplementary_analysis_map.csv`: supplementary items mapped to supporting
  material.
- `scripts/`: current analysis and figure-generation scripts.
- `results/`: current derived result tables.
- `figure_sources/`: frozen plotting inputs and final assets for main Figures 1-5
  and Supplementary Figures S1-S13.
- `metadata/`: software versions, seeds, and cell-type signature definitions.
- `checksums.csv`: SHA-256 values for all other package files.

## Scope

The package covers bulk cohort processing and meta-analysis, p53 module
scoring, donor-level single-nucleus pseudobulk models, matched raw/DecontX
pseudobulk sensitivity, MASL-referenced sensitivity analyses, composition
adjustment, single-gene AUC, incremental prediction, external-support
summaries, and figure generation.

Raw third-party datasets and large single-cell objects are not redistributed.
They should be obtained from the accessions listed in `data_sources.csv`.
Scripts accept explicit input and output directories; no author workstation path
is required.

## Running the scripts

Obtain the raw inputs first using the accessions in `data_sources.csv`. Most
scripts read an input directory from an environment variable (for example
`TP53_FIG_SOURCE_DIR`, `STATE_SOURCE`, or `PSEUDO_PATH`) and write derived
outputs to a configured output directory or to `results/`; the mapping between
reported results, scripts, and result files is in `analysis_map.csv`.

## Dependencies and installation

Base and recommended R packages (`grid`, `grDevices`, `stats`, `parallel`,
`Matrix`, `utils`) ship with R and are not listed below.

CRAN: `cowplot`, `data.table`, `ggplot2`, `ggrepel`, `gtable` (installed with
ggplot2), `gridExtra`, `lme4`, `metafor`, `patchwork`, `pheatmap`, `png`,
`pROC`, `Seurat`, `SeuratObject`, `svglite`.

Bioconductor: `GSVA`, `SingleCellExperiment`.

GitHub:

- `MuSiC` — `remotes::install_github("xuranw/MuSiC")`

```r
install.packages(c("cowplot", "data.table", "ggplot2", "ggrepel", "gridExtra",
                   "lme4", "metafor", "patchwork", "pheatmap", "png", "pROC",
                   "Seurat", "SeuratObject", "svglite"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GSVA", "SingleCellExperiment"))
```

The manuscript distinguishes association from specificity, diagnostic
performance, and causal inference. Exploratory or sensitivity analyses are
labeled as such where reported.

## Version

Reproducibility package v4.0 was prepared on 2026-11-13. Repository release
v4.0 (2026-11-13) is the de-MR manuscript package; the Zenodo DOI will be assigned
when the new version is deposited. The F0-F4
scale-sensitivity meta-analysis now uses the same `1.06/(n - 3)` Fisher-z
variance approximation as the primary and NAS meta-analyses. See `VERSION.md`.

## License

This repository is distributed under the MIT License. See `LICENSE`.
