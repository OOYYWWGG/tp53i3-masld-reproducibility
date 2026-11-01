args <- commandArgs(trailingOnly = TRUE)
donor_path <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

donor_data <- read.csv(
  donor_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

fit_family <- function(data, reference, matrix_type) {
  data$group <- factor(
    data$group,
    levels = c(
      reference,
      setdiff(c("ctrl", "MASL", "eMASH", "aMASH"), reference)
    )
  )
  fit <- glm(
    log2cpm ~ group + log(median_nCount),
    data = data
  )
  coefficients <- summary(fit)$coefficients
  intervals <- confint.default(fit)
  included_contrasts <- if (identical(reference, "MASL")) {
    c("eMASH", "aMASH")
  } else {
    c("MASL", "eMASH", "aMASH")
  }
  terms <- paste0("group", included_contrasts)
  result <- do.call(rbind, lapply(terms, function(term) {
    data.frame(
      matrix_type = matrix_type,
      reference = reference,
      contrast = sub("^group", "", term),
      n = nrow(data),
      beta_log2 = coefficients[term, "Estimate"],
      se = coefficients[term, "Std. Error"],
      ci_low = intervals[term, 1],
      ci_high = intervals[term, 2],
      p_raw = coefficients[term, ncol(coefficients)],
      ratio = 2^coefficients[term, "Estimate"],
      ratio_low = 2^intervals[term, 1],
      ratio_high = 2^intervals[term, 2],
      stringsAsFactors = FALSE
    )
  }))
  result$q_bh <- p.adjust(result$p_raw, method = "BH")
  result
}

model_rows <- list()
for (matrix_type in unique(donor_data$matrix_type)) {
  x <- donor_data[donor_data$matrix_type == matrix_type, ]
  model_rows[[paste(matrix_type, "MASL")]] <- fit_family(
    x,
    reference = "MASL",
    matrix_type = matrix_type
  )
  model_rows[[paste(matrix_type, "ctrl")]] <- fit_family(
    x,
    reference = "ctrl",
    matrix_type = matrix_type
  )
}
models <- do.call(rbind, model_rows)
rownames(models) <- NULL

comparison <- merge(
  models[
    models$matrix_type == "matched_subset_raw",
    c("reference", "contrast", "n", "ratio", "ratio_low", "ratio_high",
      "p_raw", "q_bh"),
  ],
  models[
    models$matrix_type == "matched_subset_decontx",
    c("reference", "contrast", "n", "ratio", "ratio_low", "ratio_high",
      "p_raw", "q_bh"),
  ],
  by = c("reference", "contrast", "n"),
  suffixes = c("_raw", "_decontx")
)
comparison$ratio_change_decontx_vs_raw <-
  comparison$ratio_decontx / comparison$ratio_raw
comparison <- comparison[
  order(comparison$reference, comparison$contrast),
]
rownames(comparison) <- NULL

write.csv(
  models,
  file.path(out_dir, "decontx_pseudobulk_models.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  comparison,
  file.path(out_dir, "decontx_pseudobulk_raw_vs_corrected.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
