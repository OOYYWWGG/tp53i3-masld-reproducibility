suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

result_dir <- Sys.getenv("TP53_FIG_RESULT_DIR")
source_dir <- Sys.getenv("TP53_FIG_SOURCE_DIR")
out_dir <- Sys.getenv("TP53_FIG_OUT_DIR")
stopifnot(
  nzchar(result_dir),
  nzchar(source_dir),
  nzchar(out_dir)
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

palette_group <- c(
  ctrl = "#0072B2",
  MASL = "#E69F00",
  eMASH = "#009E73",
  aMASH = "#6A3D9A"
)

theme_figure <- theme_classic(base_size = 8.2, base_family = "sans") +
  theme(
    axis.text = element_text(size = 7.2, colour = "black"),
    axis.title = element_text(size = 7.8, colour = "black"),
    axis.line = element_line(linewidth = 0.5, colour = "black"),
    plot.title = element_text(
      face = "bold",
      size = 8.8,
      hjust = 0,
      margin = margin(b = 1)
    ),
    plot.subtitle = element_text(
      size = 7.0,
      colour = "#333333",
      margin = margin(b = 1)
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 7.1),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "pt"),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks.length = unit(3.0, "pt"),
    axis.title.y = element_text(margin = margin(r = 3)),
    panel.grid.major = element_line(
      colour = "#E5E5E5",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 7, 5, 5)
  )

rank_resid_on_depth <- function(x, depth) {
  x_rank <- rank(x, ties.method = "average")
  depth_rank <- rank(depth, ties.method = "average")
  as.numeric(scale(residuals(lm(x_rank ~ depth_rank))))
}

donor <- read.csv(
  file.path(result_dir, "state_program_donor_scores.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stopifnot(nrow(donor) == 48)
donor$group <- factor(
  donor$group,
  levels = c("ctrl", "MASL", "eMASH", "aMASH")
)
donor$tp53_resid <- rank_resid_on_depth(
  donor$tp53i3_pct,
  donor$med_nCount
)
donor$p53_resid <- rank_resid_on_depth(
  donor$p53_module,
  donor$med_nCount
)
donor$zonation_resid <- rank_resid_on_depth(
  donor$zonation_module,
  donor$med_nCount
)

x_label <- bquote(
  atop(
    italic("TP53I3") * "-positive hepatocyte fraction",
    "(depth-adjusted rank residual, SD)"
  )
)

scatter_panel <- function(y, y_label, panel, rho, ci) {
  x_span <- diff(range(donor$tp53_resid, na.rm = TRUE))
  x_left <- min(donor$tp53_resid, na.rm = TRUE) - x_span * 0.14
  x_right <- max(donor$tp53_resid, na.rm = TRUE) + x_span * 0.03
  y_mid <- mean(range(donor[[y]], na.rm = TRUE))
  ggplot(
    donor,
    aes(x = tp53_resid, y = .data[[y]], colour = group)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "#C8C8C8") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "#C8C8C8") +
    geom_point(size = 1.4, alpha = 0.84) +
    annotate(
      "text",
      x = x_left - x_span * 0.28,
      y = y_mid,
      label = "Residual (SD)",
      angle = 90,
      size = 2.75,
      colour = "black"
    ) +
    scale_x_continuous(breaks = c(-2, -1, 0, 1)) +
    scale_colour_manual(values = palette_group, drop = FALSE) +
    labs(
      title = sprintf("%s  %s", panel, y_label),
      subtitle = sprintf(
        "Depth-adjusted rank residual (SD)\npartial rho = %s; 95%% CI %s; n = 48",
        rho,
        ci
      ),
      x = x_label,
      y = NULL
    ) +
    coord_cartesian(
      xlim = c(x_left, x_right),
      clip = "off"
    ) +
    theme_figure +
    theme(axis.title.y = element_blank())
}

panel_3a <- scatter_panel(
  "p53_resid",
  "p53 target-gene module",
  "A",
  "0.583",
  "0.334 to 0.768"
)
panel_3b <- scatter_panel(
  "zonation_resid",
  "Perivenous program",
  "B",
  "-0.519",
  "-0.714 to -0.263"
)

module_estimates <- read.csv(
  file.path(source_dir, "Figure3_module_cohort_estimates.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

forest_panel <- function(rows, panel, title, i2) {
  rows$row <- factor(
    rows$row,
    levels = rev(rows$row_order)
  )
  rows$estimate_type <- factor(
    ifelse(rows$pooled, "Pooled", "Cohort"),
    levels = c("Cohort", "Pooled")
  )
  ggplot(
    rows,
    aes(
      x = rho,
      y = row,
      colour = estimate_type
    )
  ) +
    geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
    geom_pointrange(
      aes(xmin = lcl, xmax = ucl),
      linewidth = 0.4,
      size = 0.28
    ) +
    scale_colour_manual(
      values = c(Cohort = "#2166AC", Pooled = "#B2182B"),
      drop = FALSE
    ) +
    scale_x_continuous(
      limits = c(-0.25, 0.9),
      breaks = c(-0.25, 0, 0.25, 0.50, 0.75)
    ) +
    labs(
      title = sprintf("%s  %s", panel, title),
      subtitle = bquote("Heterogeneity: " * I^2 == .(i2)),
      x = "Spearman rho with fibrosis (95% CI)",
      y = NULL
    ) +
    theme_figure +
    theme(
      panel.grid.major.y = element_line(
        colour = "#EEEEEE",
        linewidth = 0.25
      ),
      legend.position = if (panel == "C") "top" else "none",
      legend.direction = "horizontal"
    )
}

rows_3 <- data.frame(
  row = c(
    "GSE135251 (n = 206)",
    "GSE48452 (n = 32)",
    "GSE89632 (n = 39)",
    "Pooled (3 cohorts)"
  ),
  row_order = c(
    "GSE135251 (n = 206)",
    "GSE48452 (n = 32)",
    "GSE89632 (n = 39)",
    "Pooled (3 cohorts)"
  ),
  rho = c(0.483617, 0.473822, 0.487135, 0.483037),
  lcl = c(0.371474, 0.149893, 0.202787, 0.386104),
  ucl = c(0.581851, 0.705890, 0.695721, 0.569418),
  pooled = c(FALSE, FALSE, FALSE, TRUE)
)
rows_4 <- data.frame(
  row = c(
    "GSE281797 (n = 64)",
    "GSE89632 (n = 35)",
    "GSE48452 (n = 30)",
    "GSE135251 (n = 192)",
    "Pooled (4 cohorts)"
  ),
  row_order = c(
    "GSE281797 (n = 64)",
    "GSE89632 (n = 35)",
    "GSE48452 (n = 30)",
    "GSE135251 (n = 192)",
    "Pooled (4 cohorts)"
  ),
  rho = c(0.100173, 0.428907, 0.386283, 0.452950, 0.351000),
  lcl = c(-0.160000, 0.100000, 0.020000, 0.330000, 0.069000),
  ucl = c(0.340000, 0.670000, 0.660000, 0.560000, 0.581000),
  pooled = c(FALSE, FALSE, FALSE, FALSE, TRUE)
)

panel_3c <- forest_panel(
  rows_3,
  "C",
  "F0-F4 module correlations",
  "0"
)
panel_3d <- forest_panel(
  rows_4,
  "D",
  "F0-F3 module correlations",
  "54.3%"
)

top_row_3 <- (
  panel_3a | plot_spacer() | panel_3b
) +
  plot_layout(
    widths = c(1, 0.13, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center"
  )
bottom_row_3 <- (
  panel_3c | plot_spacer() | panel_3d
) +
  plot_layout(
    widths = c(1, 0.13, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center"
  )

figure3 <- top_row_3 / bottom_row_3 +
  plot_layout(heights = c(1, 1))

controls <- read.csv(
  file.path(source_dir, "S10_matched_panel_effects.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
controls$target <- as.logical(controls$target)
controls <- controls[order(controls$OR, decreasing = TRUE), ]
controls$gene <- factor(
  controls$gene,
  levels = rev(controls$gene)
)

panel_4a <- ggplot(
  controls,
  aes(
    x = OR,
    y = gene,
    colour = target
  )
) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#777777") +
  geom_pointrange(
    aes(xmin = lcl, xmax = ucl),
    linewidth = 0.35,
    size = 0.26
  ) +
  scale_colour_manual(
    values = c(`FALSE` = "#5B5B5B", `TRUE` = "#B2182B"),
    labels = c(`FALSE` = "Control genes", `TRUE` = "TP53I3"),
    breaks = c(FALSE, TRUE)
  ) +
  scale_x_log10(
    breaks = c(0.75, 1, 2, 4, 8),
    labels = c("0.75", "1", "2", "4", "8")
  ) +
  labs(
    title = "A  Expression-matched controls",
    subtitle = paste0(
      "TP53I3 rank 20/21; empirical P = 0.095; ",
      "10/20 controls had nominal P < 0.05"
    ),
    x = "Odds ratio, aMASH versus ctrl (95% CI)",
    y = NULL
  ) +
  theme_figure +
  theme(
    legend.position = "top",
    legend.direction = "horizontal"
  )

external <- read.csv(
  file.path(result_dir, "external_validation_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
external <- external[
  external$dataset == "GSE256398",
  ,
  drop = FALSE
]
external$fold <- as.numeric(
  sub(".*; ([0-9.]+)x.*", "\\1", external$effect)
)
external$label <- c(
  MASLD_F0 = "MASLD F0",
  MASH_Fibrosis = "MASH fibrosis",
  MASH_Cirrhosis = "MASH cirrhosis",
  Alcohol_Cirrhosis = "Alcohol cirrhosis",
  Alcohol_Hepatitis = "Alcohol hepatitis"
)[sub(" vs Healthy", "", external$comparison)]
external$display <- sprintf(
  "%s (n = %d)\n%.2fx; BH q = %.3g",
  external$label,
  external$case_n,
  external$fold,
  external$fdr
)
external$significant <- external$fdr < 0.05
external$display <- factor(
  external$display,
  levels = c(
    "MASLD F0 (n = 3)\n2.35x; BH q = 0.245",
    "MASH cirrhosis (n = 4)\n2.15x; BH q = 0.245",
    "Alcohol cirrhosis (n = 4)\n6.63x; BH q = 0.245",
    "MASH fibrosis (n = 4)\n10.45x; BH q = 0.0355",
    "Alcohol hepatitis (n = 5)\n7.80x; BH q = 0.0355"
  )
)

panel_4b <- ggplot(
  external,
  aes(
    x = fold,
    y = display,
    colour = significant
  )
) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#777777") +
  geom_point(size = 1.75) +
  scale_x_log10(
    limits = c(0.75, 13),
    breaks = c(1, 2, 4, 8),
    labels = c("1", "2", "4", "8")
  ) +
  scale_colour_manual(
    values = c(`FALSE` = "#7F7F7F", `TRUE` = "#B2182B"),
    labels = c(`FALSE` = "BH q >= 0.05", `TRUE` = "BH q < 0.05"),
    breaks = c(FALSE, TRUE)
  ) +
  labs(
    title = "B  Disease and alcohol contrasts",
    subtitle = "GSE256398; fold change versus Healthy",
    x = "Fold change versus Healthy (log scale)",
    y = NULL
  ) +
  theme_figure +
  theme(
    axis.text.y = element_text(size = 6.55, lineheight = 0.88),
    legend.position = "top",
    legend.direction = "horizontal"
  )

primary <- read.csv(
  file.path(result_dir, "pseudobulk_results.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
primary <- primary[
  primary$model == "MASL_reference_primary_pseudobulk" &
    primary$contrast %in% c("eMASH", "aMASH"),
  ,
  drop = FALSE
]
decontx <- read.csv(
  file.path(result_dir, "decontx_pseudobulk_raw_vs_corrected.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
decontx <- decontx[
  decontx$reference == "MASL" &
    decontx$contrast %in% c("eMASH", "aMASH"),
  ,
  drop = FALSE
]
masl_rows <- rbind(
  data.frame(
    specification = "Primary model",
    group = primary$contrast,
    ratio = primary$fold_change,
    lcl = primary$fold_low,
    ucl = primary$fold_high,
    q = c(0.00845, 0.00845)
  ),
  data.frame(
    specification = "Matched raw",
    group = decontx$contrast,
    ratio = decontx$ratio_raw,
    lcl = decontx$ratio_low_raw,
    ucl = decontx$ratio_high_raw,
    q = decontx$q_bh_raw
  ),
  data.frame(
    specification = "Matched DecontX",
    group = decontx$contrast,
    ratio = decontx$ratio_decontx,
    lcl = decontx$ratio_low_decontx,
    ucl = decontx$ratio_high_decontx,
    q = decontx$q_bh_decontx
  )
)
masl_rows$group <- factor(masl_rows$group, levels = c("eMASH", "aMASH"))
masl_rows$specification <- factor(
  masl_rows$specification,
  levels = c("Primary model", "Matched raw", "Matched DecontX")
)

panel_4c <- ggplot(
  masl_rows,
  aes(
    x = ratio,
    y = specification,
    colour = group,
    shape = group
  )
) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#777777") +
  geom_pointrange(
    aes(xmin = lcl, xmax = ucl),
    position = position_dodge(width = 0.42),
    linewidth = 0.42,
    size = 0.30
  ) +
  scale_x_log10(
    limits = c(0.8, 5.2),
    breaks = c(1, 2, 3, 4, 5),
    labels = c("1", "2", "3", "4", "5")
  ) +
  scale_colour_manual(values = palette_group[c("eMASH", "aMASH")]) +
  scale_shape_manual(values = c(eMASH = 16, aMASH = 17)) +
  labs(
    title = "C  MASL-referenced sensitivity",
    subtitle = "Ratio of CP10K + 0.01; 95% CI",
    x = "Ratio versus MASL (log scale)",
    y = NULL
  ) +
  theme_figure +
  theme(
    legend.position = "top",
    legend.direction = "horizontal"
  )

lower_row_4 <- (panel_4b | plot_spacer() | panel_4c) +
  plot_layout(widths = c(1, 0.13, 1))

figure4 <- panel_4a / lower_row_4 +
  plot_layout(heights = c(1.72, 1))

save_figure <- function(plot, stem, width, height) {
  ggsave(
    file.path(out_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    device = grDevices::png,
    type = "cairo",
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    file.path(out_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    file.path(out_dir, paste0(stem, ".svg")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = svglite,
    bg = "white",
    limitsize = FALSE
  )
}

save_figure(figure3, "Figure3", 6.5, 6.0)
save_figure(figure4, "Figure4", 6.5, 5.8)
cat("COMBINED_FIG3_SPECIFICITY_FIG4_OK\n")
