suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source_path <- Sys.getenv("STATE_SOURCE")
out_dir <- Sys.getenv("STATE_OUT")
stopifnot(nzchar(source_path), nzchar(out_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(source_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c(
  "group", "tp53i3_pct", "med_nCount",
  "p53_module", "zonation_module"
)
stopifnot(all(required %in% names(d)), nrow(d) == 48)

d$group <- factor(d$group, levels = c("ctrl", "MASL", "eMASH", "aMASH"))
donor_scores_columns <- c(
  "sample", "group", "n_hep", "tp53i3_pct", "med_nCount",
  "p53_targets_mean", "p53_targets_pos", "Oxidative_mean", "Oxidative_pos",
  "Perivenous_mean", "Perivenous_pos", "Mature_hep_mean", "Mature_hep_pos",
  "p53_module", "zonation_module", "p53_module_z", "zonation_module_z",
  "state_module_score"
)
write.csv(
  d[, intersect(donor_scores_columns, names(d)), drop = FALSE],
  file.path(out_dir, "state_program_donor_scores.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

rank_resid_on_depth <- function(x, depth) {
  x_rank <- rank(x, ties.method = "average")
  depth_rank <- rank(depth, ties.method = "average")
  as.numeric(scale(residuals(lm(x_rank ~ depth_rank))))
}

d$tp53_resid <- rank_resid_on_depth(d$tp53i3_pct, d$med_nCount)
d$p53_resid <- rank_resid_on_depth(d$p53_module, d$med_nCount)
d$zonation_resid <- rank_resid_on_depth(d$zonation_module, d$med_nCount)

rho_p53 <- cor(
  d$tp53_resid,
  d$p53_resid,
  method = "pearson",
  use = "complete.obs"
)
rho_zonation <- cor(
  d$tp53_resid,
  d$zonation_resid,
  method = "pearson",
  use = "complete.obs"
)
stopifnot(
  abs(rho_p53 - 0.583) < 0.001,
  abs(rho_zonation + 0.519) < 0.001
)
## Confidence intervals are frozen values from the 2,000-draw bootstrap
## (see the manuscript verification bundle); the script recomputes the point
## estimates and persists the donor-level table and correlation summary.
correlations <- data.frame(
  program = c("p53 target-gene module", "perivenous program"),
  partial_rho = c(round(rho_p53, 3), round(rho_zonation, 3)),
  ci_lower = c(0.334, -0.714),
  ci_upper = c(0.768, -0.263),
  n = c(48, 48),
  adjustment = "rank(median nCount)",
  result_role = "descriptive donor-level association",
  stringsAsFactors = FALSE
)
write.csv(
  correlations,
  file.path(out_dir, "state_program_correlations.csv"),
  row.names = FALSE,
  quote = FALSE,
  fileEncoding = "UTF-8"
)

palette <- c(
  ctrl = "#0072B2",
  MASL = "#E69F00",
  eMASH = "#009E73",
  aMASH = "#6A3D9A"
)

x_label <- bquote(
  atop(
    italic("TP53I3") * "-positive hepatocyte fraction",
    "depth-adjusted rank residual, SD"
  )
)

base_theme <- theme_classic(base_size = 8, base_family = "sans") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 7.3),
    legend.key.height = grid::unit(0.8, "lines"),
    axis.text = element_text(size = 7.2, colour = "black"),
    axis.title = element_text(size = 7.8, colour = "black"),
    axis.line = element_line(linewidth = 0.55, colour = "black"),
    axis.ticks = element_line(linewidth = 0.5, colour = "black"),
    plot.title = element_text(face = "bold", size = 8.8, hjust = 0),
    plot.subtitle = element_text(size = 7.1, colour = "#333333", hjust = 0),
    plot.margin = margin(2, 4, 2, 2)
  )

scatter_panel <- function(y, y_label, panel, rho, ci) {
  subtitle <- sprintf(
    "partial rho = %s; 95%% CI %s; n = 48",
    rho,
    ci
  )
  ggplot(d, aes(x = tp53_resid, y = .data[[y]], colour = group)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "#C8C8C8") +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#C8C8C8") +
    geom_point(size = 2.0, alpha = 0.82, stroke = 0.35) +
    scale_colour_manual(values = palette, drop = FALSE) +
    labs(
      title = sprintf("%s  %s", panel, y_label),
      subtitle = subtitle,
      x = x_label,
      y = bquote(
        atop(
          .(y_label),
          "depth-adjusted rank residual, SD"
        )
      )
    ) +
    base_theme
}

p53_plot <- scatter_panel(
  "p53_resid",
  "p53 target-gene module",
  "A",
  "0.583",
  "0.334 to 0.768"
)
zonation_plot <- scatter_panel(
  "zonation_resid",
  "Perivenous program",
  "B",
  "-0.519",
  "-0.714 to -0.263"
)

figure <- (p53_plot | zonation_plot) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center"
  )

png_path <- file.path(out_dir, "Figure3_state_programs.png")
pdf_path <- file.path(out_dir, "Figure3_state_programs.pdf")
pdf_tmp <- tempfile(pattern = "Figure3_state_programs_", fileext = ".pdf")
ggsave(
  png_path,
  figure,
  width = 6.5,
  height = 2.9,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  pdf_tmp,
  figure,
  width = 6.5,
  height = 2.9,
  units = "in",
  device = "pdf",
  bg = "white"
)
stopifnot(file.copy(pdf_tmp, pdf_path, overwrite = TRUE))
unlink(pdf_tmp)
cat(png_path, "\n")
cat(pdf_path, "\n")
