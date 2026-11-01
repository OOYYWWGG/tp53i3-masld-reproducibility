suppressPackageStartupMessages(library(pROC))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: run_gse135251_auc.R <analysis-output-directory>")
}
base_dir <- args[1]

data_path <- file.path(base_dir, "TP53I3_analysis_dataset.csv")
d <- read.csv(data_path, fileEncoding = "UTF-8-BOM", check.names = FALSE,
              stringsAsFactors = FALSE)
x <- d[d$dataset == "GSE135251", , drop = FALSE]
x$tp53i3_z <- (x$tp53i3 - mean(x$tp53i3, na.rm = TRUE)) / sd(x$tp53i3, na.rm = TRUE)

make_result <- function(comparison_id, case_definition, control_definition,
                        population, n_case, n_control, n_total_cohort,
                        auc = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                        direction = "", threshold = NA_real_,
                        sensitivity = NA_real_, specificity = NA_real_,
                        p_value = NA_real_, notes = "") {
  data.frame(
    comparison_id = comparison_id,
    case_definition = case_definition,
    control_definition = control_definition,
    population = population,
    n_case = n_case,
    n_control = n_control,
    n_included = n_case + n_control,
    n_total_cohort = n_total_cohort,
    n_excluded = n_total_cohort - n_case - n_control,
    auc = auc,
    auc_ci_low = ci_low,
    auc_ci_high = ci_high,
    direction = direction,
    threshold_log2cpm = threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    youden = ifelse(is.na(sensitivity) || is.na(specificity),
                    NA_real_, sensitivity + specificity - 1),
    p_value = p_value,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

run_auc <- function(comparison_id, case_definition, control_definition,
                    population, case_index, control_index,
                    exclusion_reason, notes = "") {
  included <- case_index | control_index
  dd <- x[included, , drop = FALSE]
  dd$case_control <- ifelse(case_index[included], 1, 0)
  n_case <- sum(dd$case_control == 1)
  n_control <- sum(dd$case_control == 0)
  n_total <- nrow(x)

  if (n_case < 2 || n_control < 2 ||
      sum(!is.na(dd$tp53i3[dd$case_control == 1])) < 2 ||
      sum(!is.na(dd$tp53i3[dd$case_control == 0])) < 2) {
    return(list(
      result = make_result(comparison_id, case_definition, control_definition,
                           population, n_case, n_control, n_total,
                           notes = paste(notes, exclusion_reason, "insufficient_n")),
      roc = NULL,
      data = dd
    ))
  }

  roc_obj <- pROC::roc(
    response = dd$case_control,
    predictor = dd$tp53i3,
    levels = c(0, 1),
    direction = "auto",
    quiet = TRUE
  )
  ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  best <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  best <- best[1, , drop = FALSE]
  direction <- ifelse(roc_obj$direction == "<", "higher_in_case", "lower_in_case")
  p_value <- suppressWarnings(
    wilcox.test(dd$tp53i3 ~ dd$case_control, exact = FALSE)$p.value
  )

  result <- make_result(
    comparison_id, case_definition, control_definition, population,
    n_case, n_control, n_total,
    auc = as.numeric(pROC::auc(roc_obj)),
    ci_low = ci[1], ci_high = ci[3],
    direction = direction,
    threshold = best$threshold,
    sensitivity = best$sensitivity,
    specificity = best$specificity,
    p_value = p_value,
    notes = paste(notes, exclusion_reason)
  )
  list(result = result, roc = roc_obj, data = dd)
}

nafl <- x$group_unified == "NAFL"
nash <- x$group_unified == "NASH"
fib <- x$fibrosis_stage_numeric

comparisons <- list()
comparisons[[1]] <- run_auc(
  "NAFL_vs_NASH_F0_F1",
  "NASH with fibrosis F0-F1",
  "NAFL",
  "NAFLD only; controls excluded",
  nash & fib %in% c(0, 1),
  nafl,
  "excluded NASH F2-F4 and controls"
)
comparisons[[2]] <- run_auc(
  "NAFL_vs_NASH_F2_F4",
  "NASH with fibrosis F2-F4",
  "NAFL",
  "NAFLD only; controls excluded",
  nash & fib %in% c(2, 3, 4),
  nafl,
  "excluded NASH F0-F1 and controls"
)
comparisons[[3]] <- run_auc(
  "F0_F1_vs_F2_F4",
  "Fibrosis F2-F4",
  "Fibrosis F0-F1",
  "NAFLD only; controls excluded",
  fib %in% c(2, 3, 4) & x$group_unified %in% c("NAFL", "NASH"),
  fib %in% c(0, 1) & x$group_unified %in% c("NAFL", "NASH"),
  "excluded controls; F0-F1 includes NAFL and NASH F0-F1"
)
comparisons[[4]] <- run_auc(
  "NAS_lt5_vs_ge5_all",
  "NAS >= 5",
  "NAS < 5",
  "all GSE135251 samples with NAS",
  !is.na(x$nas) & x$nas >= 5,
  !is.na(x$nas) & x$nas < 5,
  "no additional exclusion"
)

primary <- do.call(rbind, lapply(comparisons, `[[`, "result"))
primary$p_fdr <- p.adjust(primary$p_value, method = "BH")
primary$fdr_scope <- "GSE135251_primary_4_tests"

# Sensitivity: NAS threshold among NAFLD only, excluding controls
comparisons[[5]] <- run_auc(
  "NAS_lt5_vs_ge5_NAFLD_only",
  "NAS >= 5",
  "NAS < 5",
  "NAFLD only; controls excluded",
  x$group_unified %in% c("NAFL", "NASH") & !is.na(x$nas) & x$nas >= 5,
  x$group_unified %in% c("NAFL", "NASH") & !is.na(x$nas) & x$nas < 5,
  "excluded controls to avoid NAS=0 inflation"
)
sensitivity <- comparisons[[5]]$result
sensitivity$p_fdr <- p.adjust(sensitivity$p_value, method = "BH")
sensitivity$fdr_scope <- "GSE135251_NAS_sensitivity"

primary <- primary[, c("comparison_id", "case_definition", "control_definition",
                       "population", "n_case", "n_control", "n_included",
                       "n_total_cohort", "n_excluded", "auc", "auc_ci_low",
                       "auc_ci_high", "direction", "threshold_log2cpm",
                       "sensitivity", "specificity", "youden", "p_value",
                       "p_fdr", "fdr_scope", "notes")]
sensitivity <- sensitivity[, names(primary)[names(primary) != "fdr_scope"]]
sensitivity$fdr_scope <- "GSE135251_NAS_sensitivity"
sensitivity <- sensitivity[, names(primary)]

write.csv(primary, file.path(base_dir, "GSE135251_AUC_primary_results.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sensitivity, file.path(base_dir, "GSE135251_AUC_sensitivity_results.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

sample_flags <- x[, c("gsm", "group_unified", "disease_stage",
                      "fibrosis_stage_numeric", "nas", "tp53i3", "tp53i3_z")]
sample_flags$NAFL_vs_NASH_F0_F1 <- ifelse(
  comparisons[[1]]$data$case_control[match(sample_flags$gsm, comparisons[[1]]$data$gsm)] == 1,
  "case",
  ifelse(comparisons[[1]]$data$case_control[match(sample_flags$gsm, comparisons[[1]]$data$gsm)] == 0,
         "control", NA)
)
sample_flags$NAFL_vs_NASH_F2_F4 <- ifelse(
  comparisons[[2]]$data$case_control[match(sample_flags$gsm, comparisons[[2]]$data$gsm)] == 1,
  "case",
  ifelse(comparisons[[2]]$data$case_control[match(sample_flags$gsm, comparisons[[2]]$data$gsm)] == 0,
         "control", NA)
)
sample_flags$F0_F1_vs_F2_F4 <- ifelse(
  comparisons[[3]]$data$case_control[match(sample_flags$gsm, comparisons[[3]]$data$gsm)] == 1,
  "case",
  ifelse(comparisons[[3]]$data$case_control[match(sample_flags$gsm, comparisons[[3]]$data$gsm)] == 0,
         "control", NA)
)
sample_flags$NAS_lt5_vs_ge5_all <- ifelse(
  comparisons[[4]]$data$case_control[match(sample_flags$gsm, comparisons[[4]]$data$gsm)] == 1,
  "case",
  ifelse(comparisons[[4]]$data$case_control[match(sample_flags$gsm, comparisons[[4]]$data$gsm)] == 0,
         "control", NA)
)
sample_flags$NAS_lt5_vs_ge5_NAFLD_only <- ifelse(
  comparisons[[5]]$data$case_control[match(sample_flags$gsm, comparisons[[5]]$data$gsm)] == 1,
  "case",
  ifelse(comparisons[[5]]$data$case_control[match(sample_flags$gsm, comparisons[[5]]$data$gsm)] == 0,
         "control", NA)
)
write.csv(sample_flags, file.path(base_dir, "GSE135251_AUC_analysis_dataset.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Four-panel ROC figure for the four requested primary comparisons
plot_labels <- c(
  "NAFL vs NASH F0-F1",
  "NAFL vs NASH F2-F4",
  "F0-F1 vs F2-F4",
  "NAS <5 vs >=5"
)
plot_colors <- c("#2E86AB", "#D1495B", "#6A994E", "#7B2CBF")

draw_rocs <- function() {
  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.0, 1.5), mgp = c(2.4, 0.8, 0))
  for (i in 1:4) {
    roc_obj <- comparisons[[i]]$roc
    if (is.null(roc_obj)) {
      plot.new()
      title(main = plot_labels[i])
      next
    }
    pROC::plot.roc(
      roc_obj,
      col = plot_colors[i],
      lwd = 2.2,
      legacy.axes = TRUE,
      print.auc = TRUE,
      print.auc.cex = 0.9,
      print.auc.y = 0.12,
      xlab = "1 - Specificity",
      ylab = "Sensitivity",
      main = plot_labels[i],
      cex.main = 1.0
    )
    best <- pROC::coords(
      roc_obj, x = "best", best.method = "youden",
      ret = c("specificity", "sensitivity"), transpose = FALSE
    )[1, , drop = FALSE]
    points(1 - best$specificity, best$sensitivity, pch = 19, col = plot_colors[i], cex = 1.1)
  }
}

png(file.path(base_dir, "GSE135251_ROC_curves.png"),
    width = 2400, height = 2000, res = 300)
draw_rocs()
dev.off()

pdf(file.path(base_dir, "GSE135251_ROC_curves.pdf"),
    width = 8, height = 6.7)
draw_rocs()
dev.off()
