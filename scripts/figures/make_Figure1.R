suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(svglite)
})

out_dir <- Sys.getenv(
  "TP53_FIG_OUT_DIR",
  unset = file.path(Sys.getenv("TEMP"), "figure_legend_revision_2026-10-23")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fib <- data.frame(
  dataset = c(
    "GSE135251 (F0-F3)", "GSE48452 (F0-F3)", "GSE89632 (F0-F3)",
    "GSE281797 (F0-F3)", "GSE193066 (F0-F3)", "Primary pooled (5 F0-F3)",
    "Sensitivity pooled (3 F0-F4)"
  ),
  rho = c(0.302, 0.443, 0.389, 0.421, 0.240, 0.322, 0.355),
  lcl = c(0.164, 0.087, 0.054, 0.188, 0.046, 0.223, 0.162),
  ucl = c(0.429, 0.698, 0.646, 0.609, 0.416, 0.414, 0.521),
  type = c(
    "Cohort estimate", "Cohort estimate", "Cohort estimate",
    "Cohort estimate", "Cohort estimate",
    "Pooled estimate", "Sensitivity pooled estimate"
  ),
  stringsAsFactors = FALSE
)

nas <- data.frame(
  dataset = c(
    "GSE135251", "GSE48452", "GSE89632", "GSE281797", "GSE193066",
    "Secondary pooled (5 cohorts)"
  ),
  rho = c(0.376, 0.319, 0.567, 0.336, 0.461, 0.405),
  lcl = c(0.25, -0.04, 0.29, 0.09, 0.29, 0.304),
  ucl = c(0.49, 0.61, 0.75, 0.54, 0.60, 0.497),
  type = c(
    "Cohort estimate", "Cohort estimate", "Cohort estimate",
    "Cohort estimate", "Cohort estimate", "Pooled estimate"
  ),
  stringsAsFactors = FALSE
)

legend_levels <- c(
  "Cohort estimate",
  "Pooled estimate",
  "Sensitivity pooled estimate"
)
palette <- c(
  "Cohort estimate" = "#2C7FB8",
  "Pooled estimate" = "#B2182B",
  "Sensitivity pooled estimate" = "#6A3D9A"
)

plot_forest <- function(d, title, levels) {
  d$label <- sprintf("%.2f [%.2f, %.2f]", d$rho, d$lcl, d$ucl)
  d$dataset <- factor(d$dataset, levels = rev(levels))
  d$type <- factor(d$type, levels = legend_levels)
  ggplot(d, aes(x = rho, y = dataset, color = type)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
    geom_pointrange(aes(xmin = lcl, xmax = ucl), size = 0.42) +
    geom_text(
      aes(label = label),
      x = 0.78,
      hjust = 0,
      size = 2.9,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = palette,
      limits = legend_levels,
      drop = FALSE
    ) +
    scale_x_continuous(breaks = c(0, 0.50)) +
    coord_cartesian(xlim = c(-0.15, 1.95)) +
    labs(x = "Spearman rho (95% CI)", y = NULL, title = title) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      legend.key.width = grid::unit(1.2, "lines"),
      legend.margin = margin(0, 0, -2, 0),
      axis.text = element_text(size = 7.1),
      axis.title = element_text(size = 7.4),
      plot.title = element_text(
        face = "bold",
        size = 8.2,
        margin = margin(b = 1)
      )
    )
}

p1 <- plot_forest(
  fib,
  "A  Fibrosis stage",
  c(
    "GSE135251 (F0-F3)", "GSE48452 (F0-F3)", "GSE89632 (F0-F3)",
    "GSE281797 (F0-F3)", "GSE193066 (F0-F3)", "Primary pooled (5 F0-F3)",
    "Sensitivity pooled (3 F0-F4)"
  )
)
p2 <- plot_forest(
  nas,
  "B  NAS",
  c(
    "GSE135251", "GSE48452", "GSE89632", "GSE281797", "GSE193066",
    "Secondary pooled (5 cohorts)"
  )
)

legend <- get_legend(
  p1 + theme(legend.position = "top", legend.justification = "center")
)
panels <- (
  (p1 + theme(legend.position = "none")) |
    plot_spacer() |
    (p2 + theme(legend.position = "none"))
) +
  plot_layout(widths = c(1, 0.12, 1))
legend_row <- plot_grid(
  NULL,
  legend,
  ncol = 2,
  rel_widths = c(0.20, 0.80)
)
figure <- plot_grid(
  legend_row,
  panels,
  ncol = 1,
  rel_heights = c(0.10, 0.90)
)

ggsave(
  file.path(out_dir, "Figure1.png"),
  figure,
  width = 6.5,
  height = 2.80,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(out_dir, "Figure1.pdf"),
  figure,
  width = 6.5,
  height = 2.80,
  device = "pdf",
  bg = "white"
)
ggsave(
  file.path(out_dir, "Figure1.svg"),
  figure,
  width = 6.5,
  height = 2.80,
  units = "in",
  device = svglite,
  bg = "white"
)
cat("FIGURE1_FINAL_OK\n")
