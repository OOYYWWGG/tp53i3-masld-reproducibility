suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

src <- Sys.getenv("TP53_FIG_SOURCE_DIR")
out <- Sys.getenv("TP53_FIG_OUT_DIR")
stopifnot(nzchar(src), nzchar(out))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

rep_start <- 24084335
rep_end <- 24084417
tss <- 24084834

eqtl_source <- Sys.getenv("TP53_S2_EQTL_SOURCE")
if (!nzchar(eqtl_source)) {
  eqtl_source <- file.path(src, "S2_eqtl_loci.csv")
}
d <- read.csv(eqtl_source, check.names = FALSE)
required_columns <- c("snpId", "pos", "pValue", "dist_rep", "dist_tss")
missing_columns <- setdiff(required_columns, names(d))
if (length(missing_columns) > 0) {
  stop(
    "Missing required S2 columns: ",
    paste(missing_columns, collapse = ", ")
  )
}
d$logp <- -log10(d$pValue)
d$region <- ifelse(
  d$pos >= rep_start & d$pos <= rep_end,
  "TGYCC repeat",
  ifelse(d$pos > tss, "Upstream of TSS", "Gene body / 3'")
)
d$region <- factor(
  d$region,
  levels = c("Gene body / 3'", "Upstream of TSS", "TGYCC repeat")
)
region_colours <- c(
  "Gene body / 3'" = "#F8766D",
  "Upstream of TSS" = "#00BFC4",
  "TGYCC repeat" = "#C0392B"
)

theme_s2 <- function(base_size = 8.0) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.30),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black", linewidth = 0.42),
      axis.ticks.length = unit(2.8, "pt"),
      plot.title = element_text(
        face = "bold",
        size = 8.6,
        margin = margin(b = 1)
      ),
      plot.subtitle = element_text(
        size = 6.5,
        colour = "grey30",
        margin = margin(b = 1)
      ),
      legend.text = element_text(size = 6.7),
      legend.title = element_text(size = 6.9),
      plot.margin = margin(4, 5, 4, 4)
    )
}

main <- ggplot(d, aes(pos / 1e6, logp, colour = region)) +
  geom_vline(
    xintercept = c(rep_start, rep_end) / 1e6,
    colour = "#C0392B",
    linewidth = 0.65
  ) +
  geom_vline(
    xintercept = tss / 1e6,
    linetype = "dashed",
    colour = "grey30",
    linewidth = 0.55
  ) +
  geom_point(size = 1.8, alpha = 0.9) +
  annotate(
    "text",
    x = rep_start / 1e6,
    y = max(d$logp) * 0.98,
    label = "Repeat and TSS\nshown in panel B",
    hjust = -0.05,
    vjust = 1,
    size = 2.1,
    colour = "grey25"
  ) +
  scale_colour_manual(
    values = region_colours,
    guide = "none"
  ) +
  labs(
    title = expression(paste(
      "A  GTEx v8 liver ",
      italic(TP53I3),
      " cis-eQTL"
    )),
    subtitle = expression(paste(
      "Red: (TGYCC)16 repeat; dashed: ",
      italic(TP53I3),
      " TSS"
    )),
    x = "Position on chromosome 2 (Mb, hg38)",
    y = expression(-log[10](italic(P)))
  ) +
  theme_s2() +
  theme(legend.position = "none")

zoom_data <- d[d$pos >= rep_start - 800 & d$pos <= tss + 150, ]
zoom <- ggplot(zoom_data, aes(pos - rep_start, logp, colour = region)) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = rep_end - rep_start,
      ymin = -Inf,
      ymax = Inf
    ),
    inherit.aes = FALSE,
    fill = "#C0392B",
    alpha = 0.16,
    colour = NA
  ) +
  geom_vline(
    xintercept = tss - rep_start,
    linetype = "dashed",
    colour = "grey25",
    linewidth = 0.58
  ) +
  geom_point(size = 1.9, alpha = 0.9) +
  annotate(
    "text",
    x = (rep_end - rep_start) / 2,
    y = Inf,
    label = "Repeat",
    vjust = 1.4,
    size = 2.2,
    colour = "#9E2A1B"
  ) +
  annotate(
    "text",
    x = tss - rep_start,
    y = Inf,
    label = "TSS",
    hjust = -0.08,
    vjust = 1.4,
    size = 2.2,
    colour = "grey20"
  ) +
  coord_cartesian(
    xlim = c(-760, 650),
    clip = "off"
  ) +
  scale_colour_manual(
    values = region_colours,
    drop = TRUE,
    guide = guide_legend(
      title = "eQTL region",
      nrow = 1,
      override.aes = list(size = 2.1, alpha = 0.9)
    )
  ) +
  labs(
    title = "B  Repeat-to-TSS local zoom",
    subtitle = paste(
      "Nearest reported eQTL: 0.62 kb upstream of the repeat;",
      "repeat-to-TSS interval: 0.50 kb"
    ),
    x = "Distance from repeat start (bp)",
    y = expression(-log[10](italic(P)))
  ) +
  theme_s2() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center"
  )

figure <- main / zoom +
  plot_layout(heights = c(1.45, 1), guides = "keep")

ggsave(
  file.path(out, "FigureS2.png"),
  figure,
  width = 6.2,
  height = 4.8,
  units = "in",
  dpi = 600,
  device = grDevices::png,
  type = "cairo",
  bg = "white"
)
ggsave(
  file.path(out, "FigureS2.pdf"),
  figure,
  width = 6.2,
  height = 4.8,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(out, "FigureS2.svg"),
  figure,
  width = 6.2,
  height = 4.8,
  units = "in",
  device = svglite,
  bg = "white"
)
