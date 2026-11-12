suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: make_FigureS12.R <source-dir> <output-dir>")
}
source_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

estimates <- read.csv(
  file.path(source_dir, "FigureS12_interaction_estimates.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stage_counts <- read.csv(
  file.path(source_dir, "FigureS12_stage_counts.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cohort_order <- c("GSE281797", "GSE135251", "GSE48452", "GSE89632")
estimates$cohort <- factor(estimates$cohort, levels = rev(cohort_order))
estimates$cohort_label <- sprintf(
  "%s (%d/%d)\nP = %.3f; HC3 P = %.3f",
  estimates$cohort,
  estimates$events_label_1,
  estimates$n - estimates$events_label_1,
  estimates$interaction_p,
  estimates$hc3_robust_p
)
estimates$cohort_label <- factor(
  estimates$cohort_label,
  levels = rev(estimates$cohort_label)
)
estimates$ylo <- as.numeric(estimates$cohort)

forest_primary <- data.frame(
  cohort_label = estimates$cohort_label,
  y = estimates$ylo + 0.10,
  estimate = estimates$interaction_beta_sd_per_stage,
  low = estimates$interaction_ci_low,
  high = estimates$interaction_ci_high,
  interval = "Primary 95% CI"
)
forest_hc3 <- data.frame(
  cohort_label = estimates$cohort_label,
  y = estimates$ylo - 0.10,
  estimate = estimates$interaction_beta_sd_per_stage,
  low = estimates$hc3_robust_ci_low,
  high = estimates$hc3_robust_ci_high,
  interval = "HC3 95% CI"
)
forest <- rbind(forest_primary, forest_hc3)
forest$interval <- factor(
  forest$interval,
  levels = c("Primary 95% CI", "HC3 95% CI")
)

palette_interval <- c(
  "Primary 95% CI" = "#2166AC",
  "HC3 95% CI" = "#7F7F7F"
)
shape_interval <- c("Primary 95% CI" = 16, "HC3 95% CI" = 18)

panel_a <- ggplot(
  forest,
  aes(x = estimate, y = cohort_label, colour = interval)
) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_pointrange(
    aes(xmin = low, xmax = high, shape = interval),
    linewidth = 0.45,
    size = 0.36
  ) +
  scale_colour_manual(values = palette_interval, drop = FALSE) +
  scale_shape_manual(values = shape_interval, drop = FALSE) +
  scale_x_continuous(
    limits = c(-1.45, 2.0),
    breaks = c(-1, 0, 1, 2)
  ) +
  labs(
    title = "A  Label-by-fibrosis interaction",
    subtitle = paste0(
      "Change in standardized TP53I3 label difference per fibrosis stage; ",
      "case/control n"
    ),
    x = "Interaction beta (SD per stage; 95% CI)",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 8.2, base_family = "sans") +
  theme(
    axis.text.x = element_text(size = 7.2, colour = "black"),
    axis.text.y = element_text(size = 6.9, colour = "black", lineheight = 0.92),
    axis.title.x = element_text(size = 7.6, colour = "black"),
    axis.line = element_line(linewidth = 0.5, colour = "black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 7.0),
    legend.key.width = unit(12, "pt"),
    panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 9.0,
      hjust = 0,
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(size = 7.1, colour = "#333333"),
    plot.margin = margin(6, 9, 5, 5)
  )

stage_counts$cohort <- factor(stage_counts$cohort, levels = cohort_order)
stage_levels <- c(0, 0.5, 1, 2, 3, 3.5, 4)
stage_labels <- c("F0", "F0.5", "F1", "F2", "F3", "F3.5", "F4")
stage_counts$fibrosis_stage <- factor(
  stage_counts$fibrosis_stage,
  levels = stage_levels,
  labels = stage_labels
)
stage_counts$status <- ifelse(
  stage_counts$case_n > 0 & stage_counts$control_n > 0,
  "Both labels",
  ifelse(stage_counts$case_n > 0, "Case only", "No samples")
)
stage_counts$status <- factor(
  stage_counts$status,
  levels = c("Both labels", "Case only", "No samples")
)
stage_counts$tile_label <- ifelse(
  stage_counts$case_n > 0,
  sprintf("%d/%d", stage_counts$case_n, stage_counts$control_n),
  "—"
)

palette_overlap <- c(
  "Both labels" = "#0072B2",
  "Case only" = "#E69F00",
  "No samples" = "#F0F0F0"
)

panel_b <- ggplot(
  stage_counts,
  aes(x = fibrosis_stage, y = cohort, fill = status)
) +
  geom_tile(colour = "white", linewidth = 0.7, alpha = 0.58) +
  geom_text(
    aes(label = tile_label),
    size = 2.55,
    colour = "#222222"
  ) +
  scale_fill_manual(values = palette_overlap, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "B  Source-label overlap by fibrosis stage",
    subtitle = "MASH/NASH cases / MASL/NAFL controls",
    x = "Fibrosis stage",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 8.2, base_family = "sans") +
  theme(
    axis.text.x = element_text(size = 7.2, colour = "black"),
    axis.text.y = element_text(size = 7.4, colour = "black"),
    axis.title.x = element_text(size = 7.6, colour = "black"),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 7.0),
    legend.key.width = unit(12, "pt"),
    plot.title = element_text(
      face = "bold",
      size = 9.0,
      hjust = 0,
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(size = 7.1, colour = "#333333"),
    plot.margin = margin(5, 9, 5, 5)
  )

figure <- panel_a / panel_b +
  plot_layout(heights = c(1.18, 1.05))

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

save_figure(figure, "FigureS12_label_fibrosis_interaction", 6.5, 6.0)
cat("FIGURES12_INTERACTION_OK\n")
