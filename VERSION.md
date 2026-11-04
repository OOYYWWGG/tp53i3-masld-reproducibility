# Version

Submission reproducibility package

Prepared for the final TP53I3 manuscript and supplementary materials.

Date: 2026-11-04

## Release v1.0.3: final manuscript figure-number alignment

- Aligns the public figure inventory with the final manuscript: main Figures
  1-5 and Supplementary Figures S1-S15.
- Renames the former main Figure 6 to Supplementary Figure S15 and synchronizes
  its PNG, PDF, SVG, plotting data, and generation scripts.
- Adds Supplementary Figure S14, its interaction and stage-count tables, final
  assets, and plotting script.
- Corrects supplementary-map paths for Figures S7, S8, and S10 and updates
  current Methods and Results section references.
- Changes no estimates, sample sizes, statistical thresholds, or interpretation
  text.

Earlier changelog entries below retain the figure labels used at the time of
those releases. Current mappings are those stated above.

## v3.1 statistical consistency revision

- Changes the F0-F4 scale-sensitivity variance from `1/(n - 3)` to the
  manuscript-specified `1.06/(n - 3)`.
- Regenerates the F0-F4 pooled estimate, confidence and prediction intervals,
  cohort weights, and leave-one-out results.
- Routes outputs to `results/` through `TP53I3_RESULT_DIR`, defaulting to the
  package results directory.
- Corrects the GSE256398 data-source metadata from bulk to single-nucleus
  transcriptome, consistent with the manuscript and supplementary figure.
- Updates the analysis map and package checksums.
- Preserves the original v3 package separately for traceability.

The displayed manuscript values did not change after this correction.

The package contains only current analysis scripts, current derived results,
current figure sources, software versions, random seeds, source accessions,
and integrity checks corresponding to this manuscript version.

## Final figure layout pass

- Reflows the four-panel Figure 3 as a balanced 2 x 2 grid.
- Places Figure 6 on its own document page with its caption on the following page.
- Refines marker, line, tick, panel-spacing, and axis-title readability without changing any plotted value.
- Starts Supplementary section S8 on a new page to prevent a repeated table header from becoming an orphaned final row.
- Synchronizes all final PNG, PDF, and SVG assets with the manuscript-facing package.

## Final legend and panel-alignment pass

- Adds a direct Cohort/Pooled legend to Figure 3D.
- Aligns the visual widths of Figure 6 panels A-C without changing
  the 0-1 axis limits or plotted values.
- Adds explicit Figure S6 legends for the highlighted tissue and cell
  categories.

## Final figure layout pass

- Reflows the four-panel Figure 3 as a balanced 2 x 2 grid.
- Places Figure 6 on its own document page with its caption on the following page.
- Refines marker, line, tick, panel-spacing, and axis-title readability without changing any plotted value.
- Starts Supplementary section S8 on a new page to prevent a repeated table header from becoming an orphaned final row.
- Synchronizes all final PNG, PDF, and SVG assets with the manuscript-facing package.

## Axis-title and facet-label refinement

- Moves the Figure 3 residual definition from a detached rotated
  axis title into the panel subtitle, keeping the numerical axis next
  to the plot.
- Uses one cohort-label scheme across Figure 6 A-C and places the
  model AUC values inside the ROC panels.
- Aligns Figure 6 A/B to the full width used by panel C and keeps all
  main figures displayed at 6.5 inches wide.

## Structural split and table renumbering

- Splits the former Figure 5 into an incremental-effect figure and an out-of-fold performance figure.
- Splits the former Supplementary Figure S12 into S12 and S13.
- Renumbers Supplementary Tables S1-S18 in document order.

## JTM figure compliance pass

- Combines Supplementary Figure S1 into one A/B/C/D image and removes the
  FigureS1_part1/FigureS1_part2 assets.
- Regenerates S2 from the persistent TGYCC-repeat locus table and uses one
  legend.
- Adopts the split S12/S13 panels without an embedded total title and with
  only the used color legend.
- Synchronizes Word embeds, package source assets, submission files, maps,
  checksums, and delivery hashes.
