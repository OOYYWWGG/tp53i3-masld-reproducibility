suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

SRC <- Sys.getenv("TP53_FIG_SOURCE_DIR")
OUT <- Sys.getenv("TP53_FIG_OUT_DIR")
stopifnot(nzchar(SRC), nzchar(OUT))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

blue <- "#075A8C"
gray <- "#BDBDBD"

theme_pub <- function(base_size = 8.4) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      axis.ticks.length = unit(3.0, "pt"),
      axis.title.y = element_text(margin = margin(r = 3)),
      plot.title = element_text(
        face = "bold",
        color = "black",
        margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        color = "grey30",
        margin = margin(b = 3)
      ),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "left",
      legend.key.width = unit(7, "pt"),
      legend.key.height = unit(8, "pt"),
      legend.text = element_text(size = 6.5),
      panel.spacing = unit(8, "pt"),
      plot.margin = margin(4, 5, 4, 2)
    )
}

tissue <- read.csv(
  file.path(SRC, "S6_hpa_tissue.csv"),
  check.names = FALSE
)
tissue$is_liver <- tissue$tissue == "Liver"
tissue$tissue_label <- factor(
  tissue$tissue_label,
  levels = rev(tissue$tissue_label)
)

panel_a <- ggplot(
  tissue,
  aes(ntpm, tissue_label, fill = is_liver)
) +
  geom_col(width = 0.68) +
  geom_text(
    aes(label = sprintf("%.1f", ntpm)),
    hjust = -0.12,
    size = 2.8
  ) +
  scale_fill_manual(
    values = c(`FALSE` = gray, `TRUE` = blue),
    breaks = c(FALSE, TRUE),
    labels = c("Other tissue", "Liver"),
    name = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(
    title = expression(
      bold("A  ") * bolditalic("TP53I3") * bold(" tissue RNA levels")
    ),
    subtitle = "HPA v25; liver samples highlighted",
    x = "nTPM",
    y = NULL
  ) +
  theme_pub()

cells <- read.csv(
  file.path(SRC, "S4_S6_hpa_single_cell.csv"),
  check.names = FALSE
)
cells$is_highlight <- cells$hpa_cell_type %in%
  c("Hepatocytes", "Cholangiocytes")
cells$hpa_cell_type <- factor(
  cells$hpa_cell_type,
  levels = rev(cells$hpa_cell_type[order(cells$hpa_ncpm)])
)

panel_b <- ggplot(
  cells,
  aes(hpa_ncpm, hpa_cell_type, fill = is_highlight)
) +
  geom_col(width = 0.68) +
  geom_text(
    aes(label = sprintf("%.1f", hpa_ncpm)),
    hjust = -0.12,
    size = 2.8
  ) +
  scale_fill_manual(
    values = c(`FALSE` = gray, `TRUE` = blue),
    breaks = c(FALSE, TRUE),
    labels = c(
      "Other cells",
      "Hepatocytes /\ncholangiocytes"
    ),
    name = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.24))) +
  labs(
    title = "B  Liver-relevant cell types",
    subtitle = "HPA v25 single-cell nCPM",
    x = "nCPM",
    y = NULL
  ) +
  theme_pub()

figure <- panel_a | panel_b

ggsave(
  file.path(OUT, "FigureS6.png"),
  figure,
  width = 6.0,
  height = 3.2,
  units = "in",
  dpi = 600,
  device = grDevices::png,
  type = "cairo",
  bg = "white"
)
ggsave(
  file.path(OUT, "FigureS6.pdf"),
  figure,
  width = 6.0,
  height = 3.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(OUT, "FigureS6.svg"),
  figure,
  width = 6.0,
  height = 3.2,
  units = "in",
  device = svglite,
  bg = "white"
)
