suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(patchwork)
})

result_dir <- Sys.getenv("FIG_RESULT_DIR")
out_dir <- Sys.getenv("FIG_OUT_DIR")
package_figure_dir <- Sys.getenv("FIG_PACKAGE_DIR")
stopifnot(
  nzchar(result_dir),
  nzchar(out_dir),
  nzchar(package_figure_dir)
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
render_dir <- tempfile("tp53i3_split_figures_")
dir.create(render_dir, recursive = TRUE, showWarnings = FALSE)

delta <- read.csv(
  file.path(result_dir, "incremental_delta_auc.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
predictions <- read.csv(
  file.path(result_dir, "incremental_predictions.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
auc_results <- read.csv(
  file.path(result_dir, "incremental_auc.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
dca <- read.csv(
  file.path(result_dir, "incremental_dca.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

plot_delta <- delta[
  delta$comparison %in% c("M12-M1", "M2-M1", "M12c-M1c"),
]
plot_delta$cohort <- factor(
  plot_delta$dataset,
  levels = c("GSE135251", "GSE48452", "GSE89632", "GSE281797")
)
plot_delta$comparison_label <- factor(
  plot_delta$comparison,
  levels = c("M12-M1", "M2-M1", "M12c-M1c"),
  labels = c("M12 - M1", "M2 - M1", "M12c - M1c")
)

figure5 <- ggplot(
  plot_delta,
  aes(x = delta_auc, y = cohort, color = comparison_label)
) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_pointrange(
    aes(xmin = lcl_95, xmax = ucl_95),
    position = position_dodge(width = 0.62),
    linewidth = 0.5,
    size = 0.36
  ) +
  scale_color_manual(
    values = c("#B2182B", "#2166AC", "#7F7F7F"),
    name = NULL
  ) +
  labs(
    title = "Incremental value (subject-level bootstrap)",
    x = "Delta AUC (95% CI)",
    y = NULL
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 8.5),
    axis.text = element_text(size = 8.5),
    axis.title = element_text(size = 9),
    plot.title = element_text(
      face = "bold",
      size = 10,
      margin = margin(b = 1)
    ),
    legend.margin = margin(0, 0, 2, 0)
  )

analyses <- c("GSE281797_5x20", "GSE48452_5x20", "GSE89632_LOO")

roc_rows <- list()
calibration_rows <- list()
for (analysis_name in analyses) {
  x <- predictions[predictions$analysis == analysis_name, , drop = FALSE]
  for (model_name in c("M1", "M12")) {
    valid <- is.finite(x[[model_name]])
    if (sum(valid) < 8 || length(unique(x$y[valid])) < 2) next
    curve <- roc(
      x$y[valid],
      x[[model_name]][valid],
      quiet = TRUE,
      direction = "<"
    )
    roc_rows[[length(roc_rows) + 1]] <- data.frame(
      analysis = analysis_name,
      model = model_name,
      specificity = rev(curve$specificities),
      sensitivity = rev(curve$sensitivities),
      stringsAsFactors = FALSE
    )
  }

  for (model_name in c("M1", "M12")) {
    prediction <- pmin(pmax(x[[model_name]], 1e-6), 1 - 1e-6)
    valid <- is.finite(prediction)
    prediction_valid <- prediction[valid]
    outcome <- x$y[valid]
    breaks <- unique(quantile(
      prediction_valid,
      probs = seq(0, 1, 0.2)
    ))
    if (length(breaks) < 3) next
    bins <- cut(
      prediction_valid,
      breaks = breaks,
      include.lowest = TRUE,
      labels = FALSE
    )
    for (bin in sort(unique(bins))) {
      rows <- which(bins == bin)
      calibration_rows[[length(calibration_rows) + 1]] <- data.frame(
        analysis = analysis_name,
        model = model_name,
        predicted = mean(prediction_valid[rows]),
        observed = mean(outcome[rows]),
        stringsAsFactors = FALSE
      )
    }
  }
}

roc_plot <- do.call(rbind, roc_rows)
calibration_plot <- do.call(rbind, calibration_rows)

auc_values <- function(analysis_name) {
  values <- auc_results[
    auc_results$analysis == analysis_name &
      auc_results$model %in% c("M1", "M12"),
    ,
    drop = FALSE
  ]
  rownames(values) <- values$model
  stopifnot(all(c("M1", "M12") %in% rownames(values)))
  values
}
auc_281797 <- auc_values("GSE281797_5x20")
auc_48452 <- auc_values("GSE48452_5x20")
auc_89632 <- auc_values("GSE89632_LOO")
facet_labels <- c(
  GSE281797_5x20 = "GSE281797\n(5 events)",
  GSE48452_5x20 = "GSE48452",
  GSE89632_LOO = "GSE89632"
)

auc_labels <- data.frame(
  analysis = factor(analyses, levels = analyses),
  label = c(
    sprintf(
      "AUC\nM1 = %.3f\nM12 = %.3f",
      auc_281797["M1", "auc"],
      auc_281797["M12", "auc"]
    ),
    sprintf(
      "AUC\nM1 = %.3f\nM12 = %.3f",
      auc_48452["M1", "auc"],
      auc_48452["M12", "auc"]
    ),
    sprintf(
      "AUC\nM1 = %.3f\nM12 = %.3f",
      auc_89632["M1", "auc"],
      auc_89632["M12", "auc"]
    )
  )
)
model_colors <- c(M1 = "#2166AC", M12 = "#B2182B")

panel_a <- ggplot(
  roc_plot,
  aes(x = 1 - specificity, y = sensitivity, color = model)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey75") +
  geom_path(linewidth = 0.55) +
  facet_wrap(~analysis, labeller = as_labeller(facet_labels)) +
  geom_text(
    data = auc_labels,
    aes(x = 0.03, y = 0.97, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 1.63,
    lineheight = 0.9,
    colour = "grey20"
  ) +
  scale_color_manual(values = model_colors, name = NULL) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = TRUE
  ) +
  labs(
    title = "A  Out-of-fold ROC",
    x = "1 - specificity",
    y = "sensitivity"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.66, "cm"),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 8.5),
    axis.text = element_text(size = 8.5),
    axis.title = element_text(size = 9),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks.length = unit(3.0, "pt"),
    strip.text = element_text(size = 8.3, lineheight = 0.86),
    plot.title = element_text(
      face = "bold",
      size = 10,
      margin = margin(b = 0)
    ),
    legend.margin = margin(0, 0, 2, 0),
    plot.margin = margin(5, 7, 5, 3)
  )

panel_b <- ggplot(
  calibration_plot,
  aes(x = predicted, y = observed, color = model)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey75") +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.4) +
  facet_wrap(
    ~analysis,
    labeller = as_labeller(facet_labels)
  ) +
  scale_color_manual(values = model_colors, name = NULL) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = TRUE
  ) +
  labs(
    title = "B  Calibration (quintiles)",
    x = "Predicted probability",
    y = "Observed proportion"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.66, "cm"),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 8.5),
    axis.text = element_text(size = 8.5),
    axis.title = element_text(size = 9),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks.length = unit(3.0, "pt"),
    strip.text = element_text(size = 8.8),
    plot.title = element_text(
      face = "bold",
      size = 10,
      margin = margin(b = 0)
    ),
    legend.margin = margin(0, 0, 2, 0),
    plot.margin = margin(5, 7, 5, 3)
  )

dca_plot <- dca[dca$analysis %in% analyses, , drop = FALSE]
panel_c <- ggplot(
  dca_plot,
  aes(x = threshold, y = net_benefit, color = model)
) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey55") +
  geom_line(linewidth = 0.5) +
  facet_wrap(
    ~analysis,
    labeller = as_labeller(facet_labels)
  ) +
  scale_color_manual(
    values = c(
      M1 = "#2166AC",
      M12 = "#B2182B",
      treat_all = "#777777"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(0.2, 0.4, 0.6, 0.8),
    labels = c("0.2", "0.4", "0.6", "0.8")
  ) +
  labs(
    title = "C  Decision curve",
    x = "Threshold probability",
    y = "Net benefit (treat_none = 0)"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0.66, "cm"),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 8.5),
    axis.text = element_text(size = 8.5),
    axis.title = element_text(size = 9),
    axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks.length = unit(3.0, "pt"),
    strip.text = element_text(size = 8.8),
    plot.title = element_text(
      face = "bold",
      size = 10,
      margin = margin(b = 0)
    ),
    legend.margin = margin(0, 0, 2, 0)
  )

figure6 <- panel_a / panel_b / panel_c +
  plot_layout(heights = c(1, 1, 1))

save_figure <- function(plot, stem, width, height) {
  rendered_png <- file.path(render_dir, paste0(stem, ".png"))
  rendered_pdf <- file.path(render_dir, paste0(stem, ".pdf"))
  rendered_svg <- file.path(render_dir, paste0(stem, ".svg"))
  ggsave(
    rendered_png,
    plot,
    width = width,
    height = height,
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    rendered_pdf,
    plot,
    width = width,
    height = height,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    rendered_svg,
    plot,
    width = width,
    height = height,
    bg = "white",
    limitsize = FALSE
  )
  stopifnot(
    file.copy(rendered_png, file.path(out_dir, paste0(stem, ".png")), overwrite = TRUE),
    file.copy(rendered_pdf, file.path(out_dir, paste0(stem, ".pdf")), overwrite = TRUE),
    file.copy(rendered_svg, file.path(out_dir, paste0(stem, ".svg")), overwrite = TRUE)
  )
}

save_figure(figure5, "Figure5", 6.5, 3.0)
save_figure(figure6, "Figure6", 6.5, 8.1)

file.copy(
  file.path(render_dir, "Figure5.png"),
  file.path(package_figure_dir, "Figure5_incremental.png"),
  overwrite = TRUE
)
file.copy(
  file.path(render_dir, "Figure5.pdf"),
  file.path(package_figure_dir, "Figure5_incremental.pdf"),
  overwrite = TRUE
)
file.copy(
  file.path(render_dir, "Figure5.svg"),
  file.path(package_figure_dir, "Figure5_incremental.svg"),
  overwrite = TRUE
)
file.copy(
  file.path(render_dir, "Figure6.png"),
  file.path(package_figure_dir, "Figure6_incremental_performance.png"),
  overwrite = TRUE
)
file.copy(
  file.path(render_dir, "Figure6.pdf"),
  file.path(package_figure_dir, "Figure6_incremental_performance.pdf"),
  overwrite = TRUE
)
file.copy(
  file.path(render_dir, "Figure6.svg"),
  file.path(package_figure_dir, "Figure6_incremental_performance.svg"),
  overwrite = TRUE
)

cat("Generated split Figure 5 and Figure 6\n")
