# Figure editing notes

## Data integrity

- Rebuild figures from the frozen result tables and the R scripts in
  `scripts/figures/`. Do not edit final PNG files by hand.
- Do not rerun statistical models or change estimates, sample sizes,
  endpoints, labels, or model hierarchy when changing presentation.

## Output contract

- Save every changed main figure as a 600-dpi PNG, vector PDF, and SVG.
- Display every main figure at 6.5 inches wide in the manuscript.
- Use single combined image files for multi-panel figures.

## Panel layout

- Use balanced facet widths and an explicit horizontal gap between panels
  in the same row.
- Keep row-level chart titles nearest to the plot, while preserving readable
  space for subtitles and legends.
- For figures with a shared group encoding, use one centered legend at the
  row level and collect the legend from only one panel.

## Figure 3-specific rules

- Keep the short left-axis title `Residual (SD)` beside the numeric tick
  labels, outside the plotting area.
- Keep the full residual definition in the panel subtitle and parenthetical
  x-axis label.
- Keep Figure 3D's Cohort/Pooled legend available.
- Keep `Heterogeneity: I2 = ...` as the C/D subtitle, separate from the
  Cohort/Pooled legend.

## Supplementary Figure S15-specific rules

- Use the same cohort-label scheme in panels A-C.
- Keep AUC values as compact in-panel annotations, not in facet titles.
- Preserve the 0-1 axis range and a small boundary expansion so curves and
  calibration points at 0 or 1 are not clipped.

## Supplementary Figure S14-specific rules

- Panel A shows both conventional and HC3 interaction intervals.
- Panel B shows case/control counts by fibrosis stage and explicitly identifies
  where both source labels are observed.
- Do not relabel the exploratory interaction estimates as evidence of
  equivalence or of stage-independent effects.

## Quality control

- Visually inspect the final Word-exported page at print size.
- Verify embedded PNG hashes against the source files.
- Verify that captions and cross-references remain correct.
- Audit for text overlap, panel clipping, missing legends, and page shifts.

## Supplementary S1 and S2 rules

- S1 is a single A/B/C/D composite image with one caption and no duplicate
  S1 part assets.
- Keep the D title outside the pheatmap plotting matrix and verify vertical
  clearance from the KEGG row.
- S2 reads the persistent `S2_eqtl_loci.csv` plotting table or an explicit
  source path; it does not depend on a temporary query directory.
- S2 uses one region legend and presents the repeat/TSS interval as structural
  context without a regulatory or causal claim.
