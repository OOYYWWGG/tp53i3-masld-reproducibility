suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)
pseudo_path <- args[1]
sample_level_path <- args[2]
output_dir <- args[3]

pseudo <- read.csv(
  pseudo_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
pseudo <- pseudo[trimws(pseudo[["ct"]]) == "Hepatocyte", ]

sample_level <- read.csv(
  sample_level_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

d <- merge(
  pseudo[, c("sample", "group", "cpm10k")],
  sample_level[, c("sample", "positive_pct", "median_nCount")],
  by = "sample"
)
d$group <- factor(d$group, levels = c("ctrl", "MASL", "eMASH", "aMASH"))
d$log2_cp10k <- log2(d$cpm10k + 0.01)

pal <- c(
  ctrl = "#0072B2",
  MASL = "#E69F00",
  eMASH = "#009E73",
  aMASH = "#6A3D9A"
)

p1 <- ggplot(d, aes(group, cpm10k, color = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = NA) +
  geom_jitter(width = 0.14, size = 1.35, alpha = 0.84) +
  scale_y_log10() +
  scale_color_manual(values = pal) +
  labs(
    x = NULL,
    y = expression(paste("Hepatocyte ", italic("TP53I3"), " CP10K (log10 scale)")),
    title = "A  Hepatocyte pseudobulk\n(uncorrected, primary)"
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 7.1),
    axis.title = element_text(size = 7.4),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks.length = unit(3.0, "pt"),
    plot.title = element_text(
      size = 8.2,
      lineheight = 0.92,
      margin = margin(b = 1)
    )
  )

p2 <- ggplot(d, aes(group, positive_pct, color = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = NA) +
  geom_jitter(width = 0.14, size = 1.35, alpha = 0.84) +
  scale_color_manual(values = pal) +
  labs(
    x = NULL,
    y = expression(paste(italic("TP53I3"), "-positive hepatocyte nuclei (%)")),
    title = "B  Positive fraction (DecontX-corrected)\n(sensitivity)"
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 7.1),
    axis.title = element_text(size = 7.4),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks.length = unit(3.0, "pt"),
    plot.title = element_text(
      size = 8.2,
      lineheight = 0.92,
      margin = margin(b = 1)
    )
  )

rho <- cor.test(
  d$cpm10k,
  d$median_nCount,
  method = "spearman",
  exact = FALSE
)

p3 <- ggplot(d, aes(median_nCount, log2_cp10k, color = group)) +
  geom_point(size = 1.45, alpha = 0.86) +
  geom_smooth(method = "lm", se = FALSE, color = "grey25", linewidth = 0.6) +
  scale_color_manual(values = pal) +
  labs(
    x = "Median nCount per donor",
    y = expression(paste("Hepatocyte ", italic("TP53I3"), " log2(CP10K + 0.01)")),
    title = sprintf(
      "C  Depth check: Spearman rho=%.3f, P=%.3f",
      unname(rho$estimate),
      rho$p.value
    )
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 7.1),
    axis.title = element_text(size = 7.4),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks.length = unit(3.0, "pt"),
    plot.title = element_text(
      size = 8.2,
      lineheight = 0.92,
      margin = margin(b = 1)
    )
  )

legend_source <- ggplot(
  data.frame(
    group = factor(
      c("ctrl", "MASL", "eMASH", "aMASH"),
      levels = c("ctrl", "MASL", "eMASH", "aMASH")
    ),
    x = 1,
    y = 1
  ),
  aes(x, y, color = group)
) +
  geom_point(size = 2) +
  scale_color_manual(values = pal) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 7.1),
    legend.key.size = unit(0.35, "cm")
  )

top_row <- (p1 | plot_spacer() | p2) +
  patchwork::plot_layout(widths = c(1, 0.13, 1))

body <- top_row / p3 +
  patchwork::plot_layout(
    heights = c(1, 1.1)
  )

fig2 <- patchwork::wrap_elements(
  full = cowplot::get_legend(legend_source)
) / body +
  patchwork::plot_layout(heights = c(0.12, 1))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(output_dir, "Figure2.png"),
  plot = fig2,
  width = 6.5,
  height = 5.05,
  units = "in",
  dpi = 600,
  device = grDevices::png,
  type = "cairo",
  bg = "white"
)
ggsave(
  filename = file.path(output_dir, "Figure2.pdf"),
  plot = fig2,
  width = 6.5,
  height = 5.05,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  filename = file.path(output_dir, "Figure2.svg"),
  plot = fig2,
  width = 6.5,
  height = 5.05,
  units = "in",
  device = svglite,
  bg = "white"
)

cat("FIGURE2_600DPI_VECTOR_OK\n")
