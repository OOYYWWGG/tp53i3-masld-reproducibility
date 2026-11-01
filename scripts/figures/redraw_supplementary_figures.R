suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
  library(svglite)
  library(grid)
})

SRC <- Sys.getenv("TP53_FIG_SOURCE_DIR")
OUT <- Sys.getenv("TP53_FIG_OUT_DIR")
stopifnot(nzchar(SRC), nzchar(OUT))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

blue <- "#075A8C"
orange <- "#E69F00"
green <- "#009E73"
red <- "#B2182B"
gray <- "#BDBDBD"

theme_pub <- function(base_size = 9) {
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
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      strip.text = element_text(face = "bold", color = "black"),
      panel.spacing = unit(8, "pt"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

save_plot <- function(plot, stem, width, height) {
  ggsave(
    file.path(OUT, paste0(stem, ".png")),
    plot = plot, width = width, height = height, units = "in",
    dpi = 600, device = grDevices::png, type = "cairo", bg = "white"
  )
  ggsave(
    file.path(OUT, paste0(stem, ".pdf")),
    plot = plot, width = width, height = height, units = "in",
    device = cairo_pdf, bg = "white"
  )
  ggsave(
    file.path(OUT, paste0(stem, ".svg")),
    plot = plot, width = width, height = height, units = "in",
    device = svglite, bg = "white"
  )
}

save_grid <- function(stem, width, height, draw_fun) {
  png_file <- file.path(OUT, paste0(stem, ".png"))
  png(
    png_file,
    width = width, height = height, units = "in",
    res = 600, type = "cairo-png", bg = "white"
  )
  draw_fun()
  dev.off()
  png::writePNG(
    png::readPNG(png_file),
    target = png_file,
    dpi = 600
  )

  cairo_pdf(
    file.path(OUT, paste0(stem, ".pdf")),
    width = width, height = height, bg = "white"
  )
  draw_fun()
  dev.off()

  svglite(
    file.path(OUT, paste0(stem, ".svg")),
    width = width, height = height, bg = "white"
  )
  draw_fun()
  dev.off()
}

## S1A and S1B.
gs <- read.csv(file.path(SRC, "S1_gsva_counts.csv"), check.names = FALSE)
comparison_order_s1 <- c(
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
gs_full <- expand.grid(
  dataset = c("GSE126848", "GSE135251"),
  comparison = unique(gs$comparison),
  geneset = c("GO gene sets", "Hallmark gene sets"),
  stringsAsFactors = FALSE
)
gs_full <- merge(
  gs_full,
  gs[, c("dataset", "comparison", "geneset", "count")],
  by = c("dataset", "comparison", "geneset"),
  all.x = TRUE
)
gs_full$count[is.na(gs_full$count)] <- 0
gs_full$label <- paste0(gs_full$dataset, ": ", gs_full$comparison)
gs_full <- gs_full[gs_full$label %in% comparison_order_s1, ]
gs_full$label <- factor(gs_full$label, levels = rev(comparison_order_s1))
gs_full$geneset <- factor(
  gs_full$geneset,
  levels = c("GO gene sets", "Hallmark gene sets")
)

pA <- ggplot(gs_full, aes(count, label, fill = geneset)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.72),
    hjust = -0.12, size = 2.8, color = "grey20"
  ) +
  scale_fill_manual(values = c("GO gene sets" = "#4C78A8",
                               "Hallmark gene sets" = "#E45756")) +
  scale_x_continuous(limits = c(0, 2500), expand = expansion(mult = c(0, 0.01))) +
  labs(
    title = "A  Significant GSVA pathways per comparison",
    subtitle = "FDR < 0.05; intentionally zero combinations are omitted",
    x = "Number of significant pathways",
    y = NULL, fill = NULL
  ) +
  theme_pub()

kegg <- read.csv(file.path(SRC, "S1_kegg_top20.csv"), check.names = FALSE)
kegg$Description <- factor(kegg$Description, levels = rev(kegg$Description))
pB <- ggplot(kegg, aes(FoldEnrichment, Description)) +
  geom_segment(aes(x = 0, xend = FoldEnrichment, y = Description,
                   yend = Description), color = "grey72", linewidth = 0.35) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradientn(
    colours = c("#2166AC", "#F7F7F7", "#B2182B"),
    name = "FDR",
    guide = guide_colorbar(display = "rectangles", nbin = 101)
  ) +
  scale_size_continuous(range = c(1.2, 4.2), name = "Gene count") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "B  Candidate-pool KEGG enrichment",
    subtitle = "Top 20 terms ranked by Benjamini-Hochberg FDR",
    x = "Fold enrichment", y = NULL
  ) +
  theme_pub() +
  theme(legend.position = "right")

s1_part1 <- pA / pB +
  plot_layout(heights = c(1.08, 1.15))
save_plot(s1_part1, "FigureS1_part1", 5.984, 8.698)

## S1C and S1D.
gsea_counts <- read.csv(file.path(SRC, "S1_gsea_counts.csv"), check.names = FALSE)
gsea_counts$comparison <- gsub("_", " ", gsea_counts$comparison)
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
gsea_counts$label <- paste0(gsea_counts$dataset, ": ", gsea_counts$comparison)
gsea_counts <- gsea_counts[gsea_counts$label %in% gsea_order, ]
gsea_counts$comparison_short <- factor(
  gsea_counts$comparison,
  levels = unique(sub("^[^:]+: ", "", gsea_order))
)
gsea_counts$geneset <- factor(
  gsea_counts$geneset,
  levels = c("GO gene sets", "Hallmark gene sets")
)
pC <- ggplot(gsea_counts, aes(comparison_short, count, fill = geneset)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.76),
    vjust = -0.22, size = 2.4
  ) +
  scale_fill_manual(values = c("GO gene sets" = "#4C78A8",
                               "Hallmark gene sets" = "#E45756")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  facet_grid(. ~ dataset, scales = "free_x", space = "free_x") +
  labs(
    title = "C  Significant GSEA pathways per comparison",
    subtitle = "FDR <= 0.25",
    x = NULL, y = "Significant pathways", fill = NULL
  ) +
  theme_pub(8.4) +
  theme(
    axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1, size = 6.5),
    legend.position = "top"
  )

hall <- read.csv(file.path(SRC, "S1_hallmark_nes.csv"), check.names = FALSE)
pathways <- hall$pathway
mat <- as.matrix(hall[, -1, drop = FALSE])
rownames(mat) <- pathways
colnames(mat) <- gsub("_", " ", colnames(mat))
ph <- pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks = seq(-3, 3, length.out = 102),
  na_col = "grey85",
  border_color = "grey92",
  fontsize_row = 3.6,
  fontsize_col = 4.8,
  cellwidth = 25,
  cellheight = 5.0,
  treeheight_row = 22,
  main = "D  Hallmark NES heatmap",
  angle_col = 90,
  silent = TRUE
)

draw_part2 <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    2, 1, heights = unit(c(0.30, 0.70), "npc")
  )))
  pushViewport(viewport(layout.pos.row = 1))
  grid.draw(ggplotGrob(pC))
  popViewport()
  pushViewport(viewport(layout.pos.row = 2))
  grid.draw(ph$gtable)
  popViewport(2)
}
save_grid("FigureS1_part2", 5.984, 8.698, draw_part2)

## S4.
hpa_cells <- read.csv(file.path(SRC, "S4_S6_hpa_single_cell.csv"),
                       check.names = FALSE)
s4 <- hpa_cells[
  !is.na(hpa_cells$study_cell_type) & !is.na(hpa_cells$study_positive_fraction),
]
s4$study_percent <- as.numeric(s4$study_positive_fraction) * 100
s4$study_label <- gsub("_", " ", s4$study_cell_type)
ct4 <- cor.test(
  s4$hpa_ncpm, s4$study_percent,
  method = "spearman", exact = FALSE
)
pS4 <- ggplot(s4, aes(hpa_ncpm, study_percent)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              color = "grey45", fill = "grey85", linewidth = 0.55) +
  geom_point(size = 1.9, color = blue) +
  geom_text_repel(aes(label = study_label), size = 2.4,
                  min.segment.length = 0.15, seed = 20261026) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = expression(bold("Cell-type ranking of ") * bolditalic("TP53I3") *
                         bold(" across HPA and GSE289173")),
    subtitle = sprintf(
      "Spearman rho = %.3f, P = %.2g, n = %d shared cell types",
      unname(ct4$estimate), ct4$p.value, nrow(s4)
    ),
    x = "HPA v25 single-cell expression (nCPM)",
    y = expression("GSE289173 " * italic("TP53I3") * "-positive fraction (%)")
  ) +
  theme_pub()
save_plot(pS4, "FigureS4", 5.984, 4.445)

## S5.
gnomad <- read.csv(file.path(SRC, "S5_gnomad_metrics.csv"),
                   check.names = FALSE)
metric_value <- function(key) {
  as.numeric(gnomad$value[match(key, gnomad$metric)])
}
s5 <- data.frame(
  class = factor(c("Expected LoF", "Observed LoF"),
                 levels = c("Expected LoF", "Observed LoF")),
  count = c(metric_value("exp_lof"), metric_value("obs_lof"))
)
pS5 <- ggplot(s5, aes(class, count, fill = class)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.1f", count)), vjust = -0.45, size = 4) +
  scale_fill_manual(values = c("Expected LoF" = "#BDBDBD",
                               "Observed LoF" = blue), guide = "none") +
  scale_y_continuous(limits = c(0, 45), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = expression(bolditalic("TP53I3") *
                         bold(" loss-of-function constraint")),
    subtitle = sprintf(
      "gnomAD v4 | pLI = %.1e (~0) | LOEUF = %.2f | mis_z = %.2f",
      metric_value("pLI"), metric_value("LOEUF (oe_lof_upper)"),
      metric_value("mis_z")
    ),
    x = NULL, y = "Number of loss-of-function variants"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold", color = "black")
  )
save_plot(pS5, "FigureS5", 5.984, 4.352)

## S6.
tis <- read.csv(file.path(SRC, "S6_hpa_tissue.csv"), check.names = FALSE)
tis$is_liver <- tis$tissue == "Liver"
tis$tissue_label <- factor(tis$tissue_label, levels = rev(tis$tissue_label))
pS6a <- ggplot(tis, aes(ntpm, tissue_label, fill = is_liver)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.1f", ntpm)), hjust = -0.12, size = 2.8) +
  scale_fill_manual(values = c(`FALSE` = gray, `TRUE` = blue),
                    breaks = c(FALSE, TRUE),
                    labels = c("Other tissue", "Liver"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(
    title = expression(bold("A  ") * bolditalic("TP53I3") *
                         bold(" tissue RNA levels")),
    subtitle = "HPA v25; liver samples highlighted",
    x = "nTPM", y = NULL
  ) +
  theme_pub(8.4) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.key.width = unit(8, "pt"),
    legend.key.height = unit(8, "pt"),
    legend.text = element_text(size = 6.8),
    plot.margin = margin(4, 5, 4, 2)
  )

cell6 <- hpa_cells
cell6$is_highlight <- cell6$hpa_cell_type %in% c("Hepatocytes", "Cholangiocytes")
cell6$hpa_cell_type <- factor(
  cell6$hpa_cell_type,
  levels = rev(cell6$hpa_cell_type[order(cell6$hpa_ncpm)])
)
pS6b <- ggplot(cell6, aes(hpa_ncpm, hpa_cell_type, fill = is_highlight)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.1f", hpa_ncpm)), hjust = -0.12, size = 2.8) +
  scale_fill_manual(values = c(`FALSE` = gray, `TRUE` = blue),
                    breaks = c(FALSE, TRUE),
                    labels = c(
                      "Other cells",
                      "Hepatocytes /\ncholangiocytes"
                    ),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.24))) +
  labs(
    title = "B  Liver-relevant cell types",
    subtitle = "HPA v25 single-cell nCPM",
    x = "nCPM", y = NULL
  ) +
  theme_pub(8.4) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.key.width = unit(7, "pt"),
    legend.key.height = unit(8, "pt"),
    legend.text = element_text(size = 6.5),
    plot.margin = margin(4, 5, 4, 2)
  )
save_plot(pS6a | pS6b, "FigureS6", 6.0, 3.2)

## S7.
s7 <- read.csv(file.path(SRC, "S7_donor_pseudobulk.csv"),
               check.names = FALSE)
s7$group <- factor(
  s7$group,
  levels = c("Healthy", "MASLD_F0", "MASH_Fibrosis", "MASH_Cirrhosis")
)
levels(s7$group) <- c("Healthy", "MASLD F0", "MASH fibrosis", "MASH cirrhosis")
pal7 <- c(
  "Healthy" = "#4C93C6",
  "MASLD F0" = "#F2C14E",
  "MASH fibrosis" = "#2CB7A5",
  "MASH cirrhosis" = "#D65A5A"
)
pS7a <- ggplot(s7, aes(group, cpm10k, fill = group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.68,
               color = "grey25", linewidth = 0.35) +
  geom_jitter(position = position_jitter(width = 0.12, height = 0,
                                         seed = 202610261),
              size = 1.3, shape = 21, stroke = 0.22,
              fill = "white", color = "grey15") +
  scale_fill_manual(values = pal7, guide = "none") +
  labs(
    title = expression(bold("(A) Hepatocyte ") * bolditalic("TP53I3") *
                         bold(" by group")),
    subtitle = "GSE256398; one point per donor (n = 17)",
    x = NULL, y = expression(italic("TP53I3") * " expression (CP10K)")
  ) +
  theme_pub(7.5) +
  theme(
    axis.text.x = element_text(angle = 22, hjust = 1),
    axis.title.y = element_text(margin = margin(r = -5)),
    plot.margin = margin(4, 4, 4, 1)
  )

med7 <- aggregate(cpm10k ~ group, data = s7, FUN = median)
pS7b <- ggplot(s7, aes(group, cpm10k)) +
  geom_line(data = med7, aes(group = 1), color = "grey35", linewidth = 0.55) +
  geom_point(data = med7, size = 1.9, shape = 18, color = red) +
  geom_jitter(aes(fill = group),
              position = position_jitter(width = 0.12, height = 0,
                                         seed = 202610262),
              size = 1.15,
              shape = 21, stroke = 0.22, color = "grey15") +
  scale_fill_manual(values = pal7, guide = "none") +
  labs(
    title = "(B) Fibrosis-strata trend",
    subtitle = "Jonckheere-Terpstra P = 0.0389",
    x = NULL, y = expression("Median " * italic("TP53I3") * " expression (CP10K)")
  ) +
  theme_pub(7.5) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1))

ct7 <- read.csv(file.path(SRC, "S7_celltype.csv"), check.names = FALSE)
ct7$is_hep <- ct7$cell_type == "Hepatocyte"
ct7$cell_type <- factor(
  ct7$cell_type,
  levels = rev(ct7$cell_type[order(ct7$positive_fraction)])
)
pS7c <- ggplot(ct7, aes(positive_fraction * 100, cell_type, fill = is_hep)) +
  geom_col(width = 0.64, color = "grey30", linewidth = 0.22) +
  geom_text(
    aes(label = sprintf("%.2f%% (n=%s)", positive_fraction * 100,
                        format(n_cells, big.mark = ","))),
    hjust = -0.06, size = 2.2
  ) +
  scale_fill_manual(values = c(`FALSE` = gray, `TRUE` = "#6A3D9A"),
                    breaks = c("FALSE", "TRUE"),
                    labels = c("Other cell types", "Hepatocyte"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.52))) +
  labs(
    title = "(C) Positive fraction by cell type",
    subtitle = "GSE256398; 111,435 nuclei",
    x = expression("% " * italic("TP53I3") * "-positive cells"), y = NULL
  ) +
  theme_pub(7.5) +
  theme(
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal"
  )

cmp7 <- read.csv(file.path(SRC, "S7_comparisons.csv"), check.names = FALSE)
cmp7$label <- sub(" vs Healthy$", "", cmp7$comparison)
cmp7$label <- gsub("_", " ", cmp7$label)
cmp7$label <- factor(
  cmp7$label,
  levels = cmp7$label[order(cmp7$fold_change)]
)
cmp7$lab <- sprintf("%.2fx  BH P = %.3f", cmp7$fold_change, cmp7$bh_padj)
pS7d <- ggplot(cmp7, aes(fold_change, label)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = 1, xend = fold_change, y = label, yend = label),
               color = "grey55", linewidth = 0.55) +
  geom_point(aes(color = significance), size = 2.05) +
  geom_text(aes(label = lab), hjust = -0.08, size = 2.3) +
  scale_color_manual(
    values = c("BH padj < 0.05" = red, "ns" = "grey60"),
    breaks = c("BH padj < 0.05", "ns"),
    labels = c("BH-adjusted P < 0.05", "Not significant"),
    name = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.88))) +
  labs(
    title = "(D) Fold change versus Healthy",
    subtitle = "BH-adjusted three-comparison family",
    x = "Fold change (median pseudobulk)", y = NULL
  ) +
  theme_pub(7.5) +
  theme(
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal"
  )

s7_fig <- (pS7a | pS7b) / (pS7c | pS7d) +
  plot_layout(heights = c(1, 1))
save_plot(s7_fig, "FigureS7", 6.2, 5.3)

## S8.
roi <- read.csv(file.path(SRC, "S8_roi_plotdata.csv"), check.names = FALSE)
roi$group <- factor(roi$group, levels = c("Healthy", "MASH"))
rho8a <- cor.test(roi$TP53I3_log2, roi$p53_13gene_zmean,
                  method = "spearman", exact = FALSE)
rho8b <- cor.test(roi$TP53I3_log2, roi$fibrosis_program_zmean,
                  method = "spearman", exact = FALSE)
pS8a <- ggplot(roi, aes(p53_13gene_zmean, TP53I3_log2, color = group)) +
  geom_point(size = 1.75, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x,
              se = TRUE, color = "grey45", fill = "grey86",
              linewidth = 0.55) +
  scale_color_manual(values = c(Healthy = "#2C7FB8", MASH = "#D95F02")) +
  labs(
    title = expression(bolditalic("TP53I3") *
                         bold(" versus 13-gene p53 module")),
    subtitle = sprintf("Spearman rho = %.2f, P = %.3f",
                       unname(rho8a$estimate), rho8a$p.value),
    x = expression("13-gene p53 module (z-mean; " * italic("TP53I3") *
                     " excluded)"),
    y = expression(italic("TP53I3") * " (log2)"), color = NULL
  ) +
  theme_pub(8.2)
pS8b <- ggplot(roi, aes(fibrosis_program_zmean, TP53I3_log2, color = group)) +
  geom_point(size = 1.75, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x,
              se = TRUE, color = "grey45", fill = "grey86",
              linewidth = 0.55) +
  scale_color_manual(values = c(Healthy = "#2C7FB8", MASH = "#D95F02")) +
  labs(
    title = expression(bolditalic("TP53I3") *
                         bold(" versus fibrosis program")),
    subtitle = sprintf("Spearman rho = %.2f, P = %.3f",
                       unname(rho8b$estimate), rho8b$p.value),
    x = "Fibrosis program (z-mean)",
    y = expression(italic("TP53I3") * " (log2)"), color = NULL
  ) +
  theme_pub(8.2)
s8_fig <- (pS8a | pS8b) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
save_plot(s8_fig, "FigureS8", 6.0, 3.2)

## S10.
s10 <- read.csv(file.path(SRC, "S10_matched_panel_effects.csv"),
                check.names = FALSE)
rank_test <- read.csv(file.path(SRC, "S10_rank_test.csv"),
                      check.names = FALSE)
s10$target <- s10$target %in% c(TRUE, "True", "TRUE", "true")
gene_levels <- s10$gene[order(s10$OR)]
s10$gene <- factor(s10$gene, levels = gene_levels)
one_sided_p <- as.numeric(
  rank_test$empirical_p_one_sided[1]
)
target_rank <- as.integer(rank_test$tp53i3_rank_high_to_low[1])
pS10 <- ggplot(s10, aes(OR, gene, color = target)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0,
                 linewidth = 0.55) +
  geom_point(size = 1.95) +
  scale_x_log10(breaks = c(1, 2, 4, 8), labels = c("1", "2", "4", "8")) +
  scale_color_manual(
    values = c(`FALSE` = "grey30", `TRUE` = red),
    labels = c("Control genes", expression(italic("TP53I3") * " (target)")),
    name = NULL
  ) +
  labs(
    title = "Expression-matched non-module control panel",
    subtitle = bquote(
      italic("TP53I3") * " ascending OR rank " * .(target_rank) * "/" *
        .(nrow(s10)) * "; one-sided empirical P = " *
        .(sprintf("%.3f", one_sided_p))
    ),
    x = "OR in aMASH versus ctrl", y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold", color = "black"),
    plot.subtitle = element_text(color = "black", size = 8),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "right"
  )
save_plot(pS10, "FigureS10", 6.0, 4.7)

cat("ALL_EIGHT_VECTOR_FIGURES_OK\n")
