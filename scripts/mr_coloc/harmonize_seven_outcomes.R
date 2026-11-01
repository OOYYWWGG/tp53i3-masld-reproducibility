suppressPackageStartupMessages(library(TwoSampleMR))

data_dir <- "data/mr"
out_dir <- "results/mr"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

instrument_file <- file.path(data_dir, "retained_instruments.csv")

instruments <- read.csv(
  instrument_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM",
  colClasses = "character"
)
required_instrument_columns <- c(
  "SNP", "effect_allele", "other_allele", "eaf",
  "beta", "se", "pval", "samplesize"
)
missing_instrument_columns <- setdiff(required_instrument_columns, names(instruments))
if (length(missing_instrument_columns)) {
  stop("Missing instrument columns: ", paste(missing_instrument_columns, collapse = ", "))
}
if (nrow(instruments) != 9 || anyDuplicated(instruments$SNP)) {
  stop("Expected exactly nine unique retained instruments")
}

instruments$beta <- as.numeric(instruments$beta)
instruments$se <- as.numeric(instruments$se)
instruments$pval <- as.numeric(instruments$pval)
instruments$eaf <- as.numeric(instruments$eaf)
instruments$samplesize <- as.numeric(instruments$samplesize)

exposure <- data.frame(
  SNP = as.character(instruments$SNP),
  beta.exposure = as.numeric(instruments$beta),
  se.exposure = as.numeric(instruments$se),
  pval.exposure = as.numeric(instruments$pval),
  effect_allele.exposure = as.character(instruments$effect_allele),
  other_allele.exposure = as.character(instruments$other_allele),
  eaf.exposure = as.numeric(instruments$eaf),
  samplesize.exposure = as.numeric(instruments$samplesize),
  exposure = rep("TP53I3_eQTLGen_blood", nrow(instruments)),
  id.exposure = rep("TP53I3_eQTLGen_blood", nrow(instruments)),
  stringsAsFactors = FALSE
)

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
  file = c(
    "outcome_FinnGen_NAFLD.csv",
    "outcome_FinnGen_CIRRHOSIS_BROAD.csv",
    "outcome_FinnGen_K11_FIBROCHIRLIV.csv",
    "outcome_GCST90054782_NAFLD_4761cases.csv",
    "outcome_GCST90091033_NAFLD_8434cases.csv",
    "outcome_GCST90319877_CIRRHOSIS_15225cases.csv",
    "outcome_GCST90728570_MASLD_MTAG_N122644.csv"
  ),
  stringsAsFactors = FALSE
)

summary_rows <- vector("list", nrow(outcome_specs))
row_audit <- vector("list", nrow(outcome_specs))

for (i in seq_len(nrow(outcome_specs))) {
  label <- outcome_specs$label[i]
  outcome_file <- file.path(data_dir, outcome_specs$file[i])
  if (!file.exists(outcome_file)) stop("Missing outcome file: ", outcome_file)

  outcome_raw <- read.csv(
    outcome_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM",
    colClasses = "character"
  )
  required_outcome_columns <- c(
    "SNP", "beta.outcome", "se.outcome", "pval.outcome",
    "effect_allele.outcome", "other_allele.outcome", "eaf.outcome"
  )
  missing_outcome_columns <- setdiff(required_outcome_columns, names(outcome_raw))
  if (length(missing_outcome_columns)) {
    stop(
      label,
      " is missing columns: ",
      paste(missing_outcome_columns, collapse = ", ")
    )
  }

  numeric_outcome_columns <- c(
    "beta.outcome", "se.outcome", "pval.outcome", "eaf.outcome"
  )
  for (column in numeric_outcome_columns) {
    outcome_raw[[column]] <- suppressWarnings(as.numeric(outcome_raw[[column]]))
  }

  overlap <- outcome_raw[outcome_raw$SNP %in% exposure$SNP, , drop = FALSE]
  duplicate_instrument_rows <- overlap$SNP[duplicated(overlap$SNP)]
  if (length(duplicate_instrument_rows)) {
    stop(
      label,
      " has duplicate retained-instrument SNPs: ",
      paste(unique(duplicate_instrument_rows), collapse = ", ")
    )
  }
  invalid_numeric_rows <- sum(
    !is.finite(overlap$beta.outcome) |
      !is.finite(overlap$se.outcome) |
      !is.finite(overlap$pval.outcome) |
      !is.finite(overlap$eaf.outcome)
  )
  if (invalid_numeric_rows) {
    stop(label, " has invalid numeric values in retained-instrument rows")
  }

  outcome <- overlap
  outcome$outcome <- label
  outcome$id.outcome <- label

  harmonised <- TwoSampleMR::harmonise_data(
    exposure_dat = exposure,
    outcome_dat = outcome,
    action = 2
  )

  output_path <- file.path(out_dir, paste0("harmonised_", label, ".csv"))
  write.csv(harmonised, output_path, row.names = FALSE, fileEncoding = "UTF-8")

  kept_logical <- as.logical(harmonised$mr_keep)
  kept_logical[is.na(kept_logical)] <- FALSE

  palindromic <- if ("palindromic" %in% names(harmonised)) {
    as.logical(harmonised$palindromic)
  } else {
    rep(FALSE, nrow(harmonised))
  }
  ambiguous <- if ("ambiguous" %in% names(harmonised)) {
    as.logical(harmonised$ambiguous)
  } else {
    rep(FALSE, nrow(harmonised))
  }
  palindromic[is.na(palindromic)] <- FALSE
  ambiguous[is.na(ambiguous)] <- FALSE

  kept_harmonised <- harmonised[kept_logical, , drop = FALSE]
  missing_snp_ids <- setdiff(exposure$SNP, overlap$SNP)
  excluded_snp_ids <- if (nrow(harmonised)) {
    harmonised$SNP[!kept_logical]
  } else {
    character(0)
  }
  post_f <- (as.numeric(kept_harmonised$beta.exposure) /
               as.numeric(kept_harmonised$se.exposure))^2
  post_r2 <- sum(
    2 * as.numeric(kept_harmonised$eaf.exposure) *
      (1 - as.numeric(kept_harmonised$eaf.exposure)) *
      as.numeric(kept_harmonised$beta.exposure)^2,
    na.rm = TRUE
  )

  summary_rows[[i]] <- data.frame(
    label = label,
    outcome_file = basename(outcome_file),
    outcome_total_snps = nrow(outcome_raw),
    expected_instruments = nrow(exposure),
    overlap_snps = nrow(overlap),
    missing_from_outcome = nrow(exposure) - nrow(overlap),
    duplicate_instrument_snps = length(unique(duplicate_instrument_rows)),
    missing_snp_ids = paste(missing_snp_ids, collapse = ";"),
    excluded_snp_ids = paste(excluded_snp_ids, collapse = ";"),
    harmonised_rows = nrow(harmonised),
    kept_snps = sum(kept_logical),
    excluded_after_harmonisation = nrow(harmonised) - sum(kept_logical),
    palindromic_snps = sum(palindromic),
    ambiguous_snps = sum(ambiguous),
    min_F_kept = if (sum(kept_logical)) min(post_f, na.rm = TRUE) else NA_real_,
    max_F_kept = if (sum(kept_logical)) max(post_f, na.rm = TRUE) else NA_real_,
    total_R2_kept = if (sum(kept_logical)) post_r2 else 0,
    status = if (length(duplicate_instrument_rows)) {
      "REVIEW"
    } else if (sum(kept_logical) == nrow(exposure)) {
      "PASS"
    } else if (sum(kept_logical) == nrow(overlap)) {
      "PASS_WITH_OUTCOME_MISSING"
    } else {
      "REVIEW"
    },
    output_file = basename(output_path),
    stringsAsFactors = FALSE
  )

  audit <- data.frame(
    label = label,
    SNP = exposure$SNP,
    in_outcome = exposure$SNP %in% outcome_raw$SNP,
    in_harmonised = exposure$SNP %in% harmonised$SNP,
    palindromic = FALSE,
    ambiguous = FALSE,
    mr_keep = FALSE,
    status = "missing_from_outcome",
    stringsAsFactors = FALSE
  )
  match_index <- match(audit$SNP, harmonised$SNP)
  present <- !is.na(match_index)
  audit$palindromic[present] <- palindromic[match_index[present]]
  audit$ambiguous[present] <- ambiguous[match_index[present]]
  audit$mr_keep[present] <- kept_logical[match_index[present]]
  audit$status[present & audit$mr_keep] <- "kept"
  audit$status[present & !audit$mr_keep & audit$ambiguous] <- "excluded_ambiguous"
  audit$status[present & !audit$mr_keep & !audit$ambiguous] <- "excluded_other"
  row_audit[[i]] <- audit
}

summary_table <- do.call(rbind, summary_rows)
row_audit_table <- do.call(rbind, row_audit)
write.csv(
  summary_table,
  file.path(out_dir, "harmonization_outcome_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  row_audit_table,
  file.path(out_dir, "harmonization_snp_status.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

report <- c(
  "TP53I3 formal harmonization of nine LD-clumped instruments",
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("R version: ", R.version.string),
  paste0("TwoSampleMR version: ", as.character(packageVersion("TwoSampleMR"))),
  "Harmonization action: 2",
  "Instrument file: retained_instruments.csv",
  sprintf("Instruments: %d", nrow(exposure)),
  "",
  capture.output(print(summary_table, row.names = FALSE))
)
writeLines(
  report,
  file.path(out_dir, "harmonization_report.txt"),
  useBytes = TRUE
)

cat("Harmonization complete.\n")
