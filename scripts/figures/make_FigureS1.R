suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(svglite)
})

src <- Sys.getenv("TP53_FIG_SOURCE_DIR")
out <- Sys.getenv("TP53_FIG_OUT_DIR")
stopifnot(nzchar(src), nzchar(out))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

theme_s1 <- function(base_size = 8.0) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.28),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(2.4, "pt"),
      plot.title = element_text(face = "bold", size = 8.4, margin = margin(b = 1)),
      plot.subtitle = element_text(size = 6.2, margin = margin(b = 1)),
      legend.title = element_text(size = 6.4),
      legend.text = element_text(size = 6.2),
      legend.key.width = unit(0.65, "lines"),
      panel.spacing = unit(6, "pt"),
      plot.margin = margin(4, 5, 4, 4)
    )
}

gs <- read.csv(file.path(src, "S1_gsva_counts.csv"), check.names = FALSE)
comparison_order <- c(
  "GSE126848: NASH vs healthy",
  "GSE126848: NAFL vs obese",
  "GSE126848: NAFL vs healthy",
  "GSE126848: NASH vs obese",
  "GSE135251: NAFL vs control",
  "GSE135251: NASH vs control",
  "GSE135251: NASH vs NAFL",
  "GSE126848: NASH vs NAFL",
  "GSE126848: obese vs healthy"
)
gs$label <- paste0(gs$dataset, ": ", gs$comparison)
gs$label <- factor(gs$label, levels = rev(comparison_order))
gs$geneset <- factor(
  gs$geneset,
  levels = c("GO gene sets", "Hallmark gene sets")
)
pA <- ggplot(gs, aes(count, label, fill = geneset)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.72),
    hjust = -0.12,
    size = 2.25,
    colour = "grey20"
  ) +
  scale_fill_manual(
    values = c("GO gene sets" = "#4C78A8",
               "Hallmark gene sets" = "#E45756"),
    labels = c("GO", "Hallmark"),
    name = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = "A  Significant GSVA pathways",
    subtitle = "FDR < 0.05",
    x = "Pathways",
    y = NULL
  ) +
  theme_s1() +
  theme(
    axis.text.y = element_text(size = 5.8),
    legend.position = "top"
  )

kegg <- read.csv(file.path(src, "S1_kegg_top20.csv"), check.names = FALSE)
kegg$Description <- factor(kegg$Description, levels = rev(kegg$Description))
pB <- ggplot(kegg, aes(FoldEnrichment, Description)) +
  geom_segment(
    aes(x = 0, xend = FoldEnrichment, y = Description, yend = Description),
    colour = "grey72",
    linewidth = 0.32
  ) +
  geom_point(aes(size = Count, colour = p.adjust)) +
  scale_colour_gradientn(
    colours = c("#2166AC", "#F7F7F7", "#B2182B"),
    name = "FDR",
    guide = guide_colorbar(
      barwidth = unit(0.55, "lines"),
      barheight = unit(2.2, "lines")
    )
  ) +
  scale_size_continuous(range = c(0.9, 3.2), name = "Genes") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "B  Candidate-pool KEGG",
    subtitle = "Top 20 terms",
    x = "Fold enrichment",
    y = NULL
  ) +
  theme_s1() +
  theme(
    axis.text.y = element_text(size = 5.8),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 5.8),
    legend.key.height = unit(0.42, "lines"),
    legend.spacing.y = unit(1, "pt"),
    legend.margin = margin(0, 0, 0, 0)
  )

gsea <- read.csv(file.path(src, "S1_gsea_counts.csv"), check.names = FALSE)
gsea$comparison <- gsub("_", " ", gsea$comparison)
gsea$geneset <- factor(
  gsea$geneset,
  levels = c("GO gene sets", "Hallmark gene sets")
)
gsea_order <- c(
  "GSE126848: healthy vs NAFL",
  "GSE126848: healthy vs NASH",
  "GSE126848: healthy vs obese",
  "GSE126848: NAFL vs NASH",
  "GSE126848: obese vs NAFL",
  "GSE126848: obese vs NASH",
  "GSE135251: control vs NAFL",
  "GSE135251: control vs NASH",
  "GSE135251: NAFL vs NASH"
)
gsea$label <- paste0(gsea$dataset, ": ", gsea$comparison)
gsea <- gsea[gsea$label %in% gsea_order, ]
gsea$comparison <- factor(
  gsea$comparison,
  levels = rev(unique(sub("^[^:]+: ", "", gsea_order)))
)
pC <- ggplot(gsea, aes(count, comparison, fill = geneset)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.72),
    hjust = -0.10,
    size = 2.1
  ) +
  scale_fill_manual(
    values = c("GO gene sets" = "#4C78A8",
               "Hallmark gene sets" = "#E45756"),
    labels = c("GO", "Hallmark"),
    name = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  facet_grid(dataset ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "C  Significant GSEA pathways",
    subtitle = "FDR <= 0.25",
    x = "Pathways",
    y = NULL
  ) +
  theme_s1(8.0) +
  theme(
    axis.text.y = element_text(size = 6.1),
    strip.text = element_text(size = 6.5, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 5.9)
  )

hall <- read.csv(file.path(src, "S1_hallmark_nes.csv"), check.names = FALSE)
mat <- as.matrix(hall[, -1, drop = FALSE])
rownames(mat) <- hall$pathway
colnames(mat) <- c(
  "126848 NASH/NAFL",
  "126848 healthy/NAFL",
  "126848 healthy/NASH",
  "126848 healthy/obese",
  "126848 obese/NAFL",
  "126848 obese/NASH",
  "135251 NAFL/NASH",
  "135251 control/NAFL",
  "135251 control/NASH"
)
ph <- pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks = seq(-3, 3, length.out = 102),
  na_col = "grey85",
  border_color = "grey90",
  fontsize = 6.2,
  fontsize_row = 5.7,
  fontsize_col = 5.7,
  cellwidth = 9.5,
  cellheight = 5.2,
  treeheight_row = 12,
  treeheight_col = 0,
  angle_col = 90,
  silent = TRUE
)

top_row <- pA | pB
pD_table <- ph$gtable
pD_table <- gtable::gtable_add_rows(
  pD_table,
  heights = unit(10, "pt"),
  pos = 0
)
pD_table <- gtable::gtable_add_grob(
  pD_table,
  grid::textGrob(
    "D  Hallmark NES heatmap",
    x = 0,
    hjust = 0,
    gp = grid::gpar(
      fontsize = 8.4,
      fontface = "bold",
      col = "black"
    )
  ),
  t = 1,
  l = 1,
  r = ncol(pD_table),
  clip = "off",
  name = "panel-D-title"
)
pD <- wrap_elements(full = pD_table, clip = FALSE)
bottom_row <- (pC | pD) +
  plot_layout(widths = c(0.8, 1.2))
figure <- top_row / plot_spacer() / bottom_row +
  plot_layout(heights = c(2.70, 0.18, 4.20)) &
  theme(plot.title = element_text(margin = margin(b = 1)))

ggsave(
  file.path(out, "FigureS1.png"),
  figure,
  width = 6.2,
  height = 8.20,
  units = "in",
  dpi = 600,
  device = grDevices::png,
  type = "cairo",
  bg = "white"
)
ggsave(
  file.path(out, "FigureS1.pdf"),
  figure,
  width = 6.2,
  height = 8.20,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(out, "FigureS1.svg"),
  figure,
  width = 6.2,
  height = 8.20,
  units = "in",
  device = svglite,
  bg = "white"
)
