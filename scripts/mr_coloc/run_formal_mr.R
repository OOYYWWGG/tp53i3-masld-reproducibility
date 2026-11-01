suppressPackageStartupMessages(library(TwoSampleMR))
suppressPackageStartupMessages(library(MRPRESSO))

data_dir <- "data/mr"
harm_dir <- "results/mr"
out_dir <- "results/mr"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

outcome_specs <- data.frame(
  label = c(
    "FinnGen_NAFLD",
    "FinnGen_CIRRHOSIS_BROAD",
    "FinnGen_K11_FIBROCHIRLIV",
    "GCST90054782_NAFLD_4761cases",
    "GCST90091033_NAFLD_8434cases",
    "GCST90319877_CIRRHOSIS_15225cases",
    "GCST90728570_MASLD_MTAG_N122644"
  ),
  harmonised_file = c(
    "harmonised_FinnGen_NAFLD.csv",
    "harmonised_FinnGen_CIRRHOSIS_BROAD.csv",
    "harmonised_FinnGen_K11_FIBROCHIRLIV.csv",
    "harmonised_GCST90054782_NAFLD_4761cases.csv",
    "harmonised_GCST90091033_NAFLD_8434cases.csv",
    "harmonised_GCST90319877_CIRRHOSIS_15225cases.csv",
    "harmonised_GCST90728570_MASLD_MTAG_N122644.csv"
  ),
  samplesize_outcome = c(
    3504 + 496844,
    5545 + 494803,
    2653 + 485213,
    4761 + 373227,
    8434 + 770180,
    15225 + 1564786,
    122644
  ),
  outcome_scale = c(
    "binary_log_odds",
    "binary_log_odds",
    "binary_log_odds",
    "binary_log_odds",
    "binary_log_odds",
    "binary_log_odds",
    "effective_N_log_odds"
  ),
  stringsAsFactors = FALSE
)

numeric_columns <- c(
  "beta.exposure", "beta.outcome", "se.exposure", "se.outcome",
  "pval.exposure", "pval.outcome", "eaf.exposure", "eaf.outcome",
  "samplesize.exposure", "samplesize.outcome"
)

method_rows <- list()
heterogeneity_rows <- list()
pleiotropy_rows <- list()
leave_one_out_rows <- list()
single_snp_rows <- list()
presso_rows <- list()
directionality_rows <- list()
analysis_input_rows <- list()
run_notes <- character(0)

for (i in seq_len(nrow(outcome_specs))) {
  label <- outcome_specs$label[i]
  input_path <- file.path(harm_dir, outcome_specs$harmonised_file[i])
  if (!file.exists(input_path)) stop("Missing harmonized input: ", input_path)

  dat <- read.csv(
    input_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM",
    colClasses = "character"
  )
  for (column in intersect(numeric_columns, names(dat))) {
    dat[[column]] <- suppressWarnings(as.numeric(dat[[column]]))
  }
  keep <- as.logical(dat$mr_keep)
  keep[is.na(keep)] <- FALSE
  dat <- dat[keep, , drop = FALSE]
  dat$mr_keep <- TRUE
  if (!nrow(dat)) stop(label, " has no harmonized instruments with mr_keep=TRUE")
  if (anyDuplicated(dat$SNP)) stop(label, " has duplicate retained SNPs")

  dat$samplesize.outcome <- outcome_specs$samplesize_outcome[i]
  dat$units.exposure <- "SD"
  dat$units.outcome <- "log odds"
  analysis_input_path <- file.path(
    out_dir, paste0("mr_analysis_input_", label, ".csv")
  )
  write.csv(dat, analysis_input_path, row.names = FALSE, fileEncoding = "UTF-8")
  analysis_input_rows[[i]] <- data.frame(
    outcome = label,
    harmonised_file = input_path,
    analysis_input_file = analysis_input_path,
    n_snp = nrow(dat),
    samplesize_outcome = outcome_specs$samplesize_outcome[i],
    outcome_scale = outcome_specs$outcome_scale[i],
    min_F = min((dat$beta.exposure / dat$se.exposure)^2),
    max_F = max((dat$beta.exposure / dat$se.exposure)^2),
    total_R2 = sum(
      2 * dat$eaf.exposure * (1 - dat$eaf.exposure) * dat$beta.exposure^2
    ),
    stringsAsFactors = FALSE
  )

  method_list <- c(
    "mr_ivw",
    "mr_weighted_median",
    "mr_egger_regression",
    "mr_raps"
  )
  res <- withCallingHandlers(
    TwoSampleMR::mr(dat, method_list = method_list),
    warning = function(w) {
      run_notes <<- c(run_notes, paste(label, "method warning:", conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  res <- TwoSampleMR::generate_odds_ratios(res)
  res$outcome_label <- label
  res$n_snp <- if ("nsnp" %in% names(res)) res$nsnp else nrow(dat)
  method_rows[[i]] <- res

  het <- TwoSampleMR::mr_heterogeneity(
    dat,
    method_list = c("mr_ivw", "mr_egger_regression")
  )
  het$outcome_label <- label
  heterogeneity_rows[[i]] <- het

  pleio <- TwoSampleMR::mr_pleiotropy_test(dat)
  pleio$outcome_label <- label
  pleiotropy_rows[[i]] <- pleio

  loo <- TwoSampleMR::mr_leaveoneout(dat)
  loo$outcome_label <- label
  leave_one_out_rows[[i]] <- loo

  singlesnp <- TwoSampleMR::mr_singlesnp(dat)
  singlesnp$outcome_label <- label
  single_snp_rows[[i]] <- singlesnp

  presso <- tryCatch(
    suppressWarnings(
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        data = dat,
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        NbDistribution = 1000,
        seed = 20261025
      )
    ),
    error = function(e) e
  )

  if (inherits(presso, "error")) {
    presso_rows[[i]] <- data.frame(
      outcome_label = label,
      n_snp = nrow(dat),
      presso_raw_beta = NA_real_,
      presso_raw_se = NA_real_,
      presso_raw_p = NA_real_,
      presso_outlier_corrected_beta = NA_real_,
      presso_outlier_corrected_p = NA_real_,
      presso_global_rss = NA_real_,
      presso_global_p = NA_real_,
      presso_outliers = NA_character_,
      status = "ERROR",
      message = conditionMessage(presso),
      stringsAsFactors = FALSE
    )
  } else {
    main <- presso[["Main MR results"]]
    global <- presso[["MR-PRESSO results"]][["Global Test"]]
    outlier_test <- presso[["MR-PRESSO results"]][["Outlier Test"]]
    raw <- main[tolower(trimws(main[["MR Analysis"]])) == "raw", , drop = FALSE]
    corrected <- main[
      tolower(trimws(main[["MR Analysis"]])) == "outlier-corrected",
      ,
      drop = FALSE
    ]
    outliers <- character(0)
    if (!is.null(outlier_test) && nrow(outlier_test)) {
      p_col <- grep("P.value|P-value|P", names(outlier_test), value = TRUE)[1]
      if (!is.na(p_col)) {
        outliers <- rownames(outlier_test)[
          !is.na(outlier_test[[p_col]]) & outlier_test[[p_col]] < 0.05
        ]
      }
    }
    presso_rows[[i]] <- data.frame(
      outcome_label = label,
      n_snp = nrow(dat),
      presso_raw_beta = raw[["Causal Estimate"]][1],
      presso_raw_se = raw[["Sd"]][1],
      presso_raw_p = raw[["P-value"]][1],
      presso_outlier_corrected_beta = corrected[["Causal Estimate"]][1],
      presso_outlier_corrected_p = corrected[["P-value"]][1],
      presso_global_rss = global[["RSSobs"]],
      presso_global_p = global[["Pvalue"]],
      presso_outliers = paste(outliers, collapse = ";"),
      status = "OK",
      message = "",
      stringsAsFactors = FALSE
    )
  }

  direction <- tryCatch(
    TwoSampleMR::directionality_test(dat),
    error = function(e) e
  )
  if (inherits(direction, "error")) {
    directionality_rows[[i]] <- data.frame(
      outcome_label = label,
      n_snp = nrow(dat),
      snp_r2.exposure = NA_real_,
      snp_r2.outcome = NA_real_,
      correct_causal_direction = NA,
      steiger_pval = NA_real_,
      status = "ERROR",
      message = conditionMessage(direction),
      stringsAsFactors = FALSE
    )
  } else {
    direction$outcome_label <- label
    direction$n_snp <- nrow(dat)
    direction$status <- "OK"
    direction$message <- ""
    directionality_rows[[i]] <- direction
  }
}

method_table <- do.call(rbind, method_rows)
heterogeneity_table <- do.call(rbind, heterogeneity_rows)
pleiotropy_table <- do.call(rbind, pleiotropy_rows)
leave_one_out_table <- do.call(rbind, leave_one_out_rows)
single_snp_table <- do.call(rbind, single_snp_rows)
presso_table <- do.call(rbind, presso_rows)
directionality_table <- do.call(rbind, directionality_rows)
directionality_table$method_note <- paste(
  "Approximate SNP-exposure and SNP-outcome correlations;",
  "binary-outcome r values were not supplied because prevalence was unavailable."
)
analysis_input_table <- do.call(rbind, analysis_input_rows)

write.csv(
  analysis_input_table,
  file.path(out_dir, "mr_analysis_inputs.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  method_table,
  file.path(out_dir, "mr_method_results_all.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  method_table[method_table$method == "Inverse variance weighted", ],
  file.path(out_dir, "mr_ivw_results.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  heterogeneity_table,
  file.path(out_dir, "mr_heterogeneity.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  pleiotropy_table,
  file.path(out_dir, "mr_egger_pleiotropy.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  leave_one_out_table,
  file.path(out_dir, "mr_leave_one_out_all.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  single_snp_table,
  file.path(out_dir, "mr_single_snp_all.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  presso_table,
  file.path(out_dir, "mr_presso.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  directionality_table,
  file.path(out_dir, "mr_directionality.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

report_lines <- c(
  "TP53I3 formal MR using summary statistics",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("R version: ", R.version.string),
  paste0("TwoSampleMR: ", as.character(packageVersion("TwoSampleMR"))),
  paste0("MRPRESSO: ", as.character(packageVersion("MRPRESSO"))),
  paste0("mr.raps: ", as.character(packageVersion("mr.raps"))),
  "Instrument source: formal PLINK LD clumping",
  "Harmonization action: 2",
  "MR-RAPS was run through TwoSampleMR::mr(method_list='mr_raps')",
  "MR-PRESSO seed: 20261025; NbDistribution: 1000",
  "",
  "MR method results:",
  capture.output(print(method_table, row.names = FALSE)),
  "",
  "Heterogeneity:",
  capture.output(print(heterogeneity_table, row.names = FALSE)),
  "",
  "Egger pleiotropy:",
  capture.output(print(pleiotropy_table, row.names = FALSE)),
  "",
  "MR-PRESSO:",
  capture.output(print(presso_table, row.names = FALSE)),
  "",
  "Directionality checks:",
  capture.output(print(directionality_table, row.names = FALSE)),
  "",
  "Warnings and notes:",
  if (length(run_notes)) run_notes else "none"
)
writeLines(
  report_lines,
  file.path(out_dir, "formal_mr_report.txt"),
  useBytes = TRUE
)

cat("Formal MR complete.\n")
