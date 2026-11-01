suppressPackageStartupMessages({
  library(stats)
})

OUT <- Sys.getenv("TP53I3_MASL_MASH_OUTPUT_DIR", unset = "results/MASL_MASH_sensitivity/reproduction")
INPUT <- Sys.getenv("TP53I3_MASL_MASH_INPUT_DIR", unset = "results/MASL_MASH_sensitivity/inputs")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PSEUDO_PATH <- file.path(
  INPUT,
  "celltype_pseudobulk.csv"
)
SAMPLE_PATH <- file.path(
  INPUT,
  "sample_level_hepatocyte.csv"
)
PUBLISHED_PATH <- file.path(INPUT, "published_primary_results.csv")

B_PERM <- 20000L
PERM_SEED <- 202610261L
N_MDE <- 10000L
MDE_SEED <- 202610262L
MDE_GRID <- seq(0, 3, by = 0.025)
CONTRASTS <- c("eMASH", "aMASH")

pseudo <- read.csv(PSEUDO_PATH, check.names = FALSE, stringsAsFactors = FALSE)
pseudo <- pseudo[trimws(pseudo$ct) == "Hepatocyte", c("sample", "group", "cpm10k")]
sl <- read.csv(
  SAMPLE_PATH,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

d <- merge(
  pseudo,
  sl[, c("sample", "median_nCount")],
  by = "sample",
  all = FALSE,
  sort = FALSE
)
d$group <- factor(d$group, levels = c("MASL", "ctrl", "eMASH", "aMASH"))
d$log2cpm <- log2(d$cpm10k + 0.01)
d$log_ncount <- log(d$median_nCount)

if (nrow(d) != 48L) stop("Expected 48 donor-level hepatocyte rows")
if (!identical(as.integer(table(d$group)), c(6L, 23L, 9L, 10L))) {
  stop("Unexpected group counts")
}
if (!("FL511" %in% d$sample)) stop("FL511 is missing")

fit_primary <- function(dat) {
  dat$group <- factor(dat$group, levels = c("MASL", "ctrl", "eMASH", "aMASH"))
  lm(log2cpm ~ group + log_ncount, data = dat)
}

extract_contrasts <- function(model, label) {
  co <- summary(model)$coefficients
  ci <- confint.default(model)
  do.call(
    rbind,
    lapply(CONTRASTS, function(target) {
      term <- paste0("group", target)
      data.frame(
        analysis = label,
        contrast = paste(target, "vs MASL"),
        target_group = target,
        n_total = nobs(model),
        n_target = sum(model$model$group == target),
        n_masl = sum(model$model$group == "MASL"),
        beta_log2 = unname(co[term, "Estimate"]),
        se = unname(co[term, "Std. Error"]),
        ci_low_log2 = unname(ci[term, 1]),
        ci_high_log2 = unname(ci[term, 2]),
        fold_change = 2^unname(co[term, "Estimate"]),
        fold_low = 2^unname(ci[term, 1]),
        fold_high = 2^unname(ci[term, 2]),
        raw_p = unname(co[term, 4]),
        df_residual = df.residual(model),
        residual_sigma = sigma(model),
        stringsAsFactors = FALSE
      )
    })
  )
}

primary_model <- fit_primary(d)
primary <- extract_contrasts(primary_model, "primary_full_data")
primary$bh_q <- p.adjust(primary$raw_p, method = "BH")

published <- read.csv(
  PUBLISHED_PATH,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
published_masl <- published[
  published$model == "MASL_reference_primary_pseudobulk" &
    published$contrast %in% CONTRASTS,
]
published_masl <- published_masl[match(CONTRASTS, published_masl$contrast), ]
primary$published_beta_log2 <- published_masl$beta_log2
primary$published_raw_p <- published_masl$p
primary$reproduction_max_abs_diff <- pmax(
  abs(primary$beta_log2 - primary$published_beta_log2),
  abs(primary$raw_p - primary$published_raw_p)
)
if (any(primary$reproduction_max_abs_diff > 1e-10)) {
  stop("Primary model reproduction failed")
}

exact_wilcoxon <- function(dat, target) {
  x <- dat$cpm10k[dat$group == target]
  y <- dat$cpm10k[dat$group == "MASL"]
  tied <- any(duplicated(c(x, y)))
  if (!tied) {
    test <- wilcox.test(x, y, alternative = "two.sided", exact = TRUE)
    return(data.frame(
      method = "exact Wilcoxon rank-sum",
      statistic = unname(test$statistic),
      p = test$p.value,
      actual_alpha = 0.05,
      stringsAsFactors = FALSE
    ))
  }
  set.seed(PERM_SEED + match(target, CONTRASTS))
  observed <- median(x) - median(y)
  combined <- c(x, y)
  n_x <- length(x)
  perm <- replicate(B_PERM, {
    idx <- sample.int(length(combined), n_x)
    median(combined[idx]) - median(combined[-idx])
  })
  data.frame(
    method = "Monte Carlo median-difference permutation (ties present)",
    statistic = observed,
    p = (1 + sum(abs(perm) >= abs(observed))) / (B_PERM + 1),
    actual_alpha = 0.05,
    stringsAsFactors = FALSE
  )
}

full_wilcoxon <- do.call(
  rbind,
  lapply(CONTRASTS, function(target) {
    cbind(
      data.frame(contrast = paste(target, "vs MASL"), target_group = target,
                 analysis = "full_data"),
      exact_wilcoxon(d, target)
    )
  })
)
full_wilcoxon$direct_pair_bh_q <- p.adjust(full_wilcoxon$p, method = "BH")
full_pairwise <- pairwise.wilcox.test(
  d$cpm10k,
  d$group,
  p.adjust.method = "BH",
  exact = TRUE
)
full_wilcoxon$six_pair_bh_q <- vapply(
  CONTRASTS,
  function(target) full_pairwise$p.value[target, "MASL"],
  numeric(1)
)

d_loo <- d[d$sample != "FL511", , drop = FALSE]
loo_model <- fit_primary(d_loo)
loo <- extract_contrasts(loo_model, "loo_without_FL511")
loo$bh_q <- p.adjust(loo$raw_p, method = "BH")
loo_wilcoxon <- do.call(
  rbind,
  lapply(CONTRASTS, function(target) {
    cbind(
      data.frame(contrast = paste(target, "vs MASL"), target_group = target,
                 analysis = "loo_without_FL511"),
      exact_wilcoxon(d_loo, target)
    )
  })
)
loo_wilcoxon$direct_pair_bh_q <- p.adjust(loo_wilcoxon$p, method = "BH")
loo_pairwise <- pairwise.wilcox.test(
  d_loo$cpm10k,
  d_loo$group,
  p.adjust.method = "BH",
  exact = TRUE
)
loo_wilcoxon$six_pair_bh_q <- vapply(
  CONTRASTS,
  function(target) loo_pairwise$p.value[target, "MASL"],
  numeric(1)
)

X_full <- model.matrix(~ group + log_ncount, data = d)
y <- d$log2cpm
qr_full <- qr(X_full)
coef_full <- qr.coef(qr_full, y)
Q <- qr.Q(qr_full)
R <- qr.R(qr_full)

target_specific_null <- function(target) {
  target_col <- match(paste0("group", target), colnames(X_full))
  reduced_cols <- setdiff(seq_len(ncol(X_full)), target_col)
  X_reduced <- X_full[, reduced_cols, drop = FALSE]
  coef_reduced <- qr.coef(qr(X_reduced), y)
  fitted_reduced <- as.numeric(X_reduced %*% coef_reduced)
  residual_reduced <- y - fitted_reduced
  set.seed(PERM_SEED + match(target, CONTRASTS))
  perm <- numeric(B_PERM)
  for (b in seq_len(B_PERM)) {
    y_perm <- fitted_reduced + sample(residual_reduced, replace = FALSE)
    perm[b] <- unname(qr.coef(qr_full, y_perm)[target_col])
  }
  list(
    target = target,
    target_col = target_col,
    fitted_reduced = fitted_reduced,
    residual_reduced = residual_reduced,
    observed = unname(coef_full[target_col]),
    perm = perm,
    null_center = mean(perm),
    permutation_p = (
      1 + sum(abs(perm - mean(perm)) >= abs(coef_full[target_col]))
    ) / (B_PERM + 1),
    null_abs_q95 = as.numeric(quantile(
      abs(perm - mean(perm)),
      probs = 0.95,
      type = 1,
      names = FALSE
    ))
  )
}

nulls <- setNames(lapply(CONTRASTS, target_specific_null), CONTRASTS)
permutation <- do.call(
  rbind,
  lapply(CONTRASTS, function(target) {
    z <- nulls[[target]]
    data.frame(
      analysis = "target_specific_Freedman_Lane_permutation",
      contrast = paste(target, "vs MASL"),
      target_group = target,
      observed_beta_log2 = z$observed,
      permutation_null_center = z$null_center,
      permutation_p = z$permutation_p,
      null_abs_q95_log2 = z$null_abs_q95,
      permutations = B_PERM,
      seed = PERM_SEED + match(target, CONTRASTS),
      stringsAsFactors = FALSE
    )
  })
)
permutation$bh_q <- p.adjust(permutation$permutation_p, method = "BH")

simulate_permutation_power <- function(target) {
  z <- nulls[[target]]
  set.seed(MDE_SEED + match(target, CONTRASTS))
  n <- length(y)
  residual_centered <- z$residual_reduced - mean(z$residual_reduced)
  residual_matrix <- matrix(
    sample(residual_centered, n * N_MDE, replace = TRUE),
    nrow = n,
    ncol = N_MDE
  )
  base_y <- z$fitted_reduced + residual_matrix
  x_target <- X_full[, z$target_col]
  rows <- vector("list", length(MDE_GRID))
  for (i in seq_along(MDE_GRID)) {
    beta <- MDE_GRID[i]
    y_sim <- base_y + x_target * beta
    coef_sim <- backsolve(R, crossprod(Q, y_sim))
    beta_sim <- coef_sim[z$target_col, ]
    reject <- abs(beta_sim - z$null_center) >= z$null_abs_q95
    rows[[i]] <- data.frame(
      target_group = target,
      contrast = paste(target, "vs MASL"),
      method = "permutation_calibrated_model_design",
      beta_log2 = beta,
      fold_change = 2^beta,
      power = mean(reject),
      simulations = N_MDE,
      seed = MDE_SEED + match(target, CONTRASTS),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

simulate_wilcoxon_power <- function(target) {
  set.seed(MDE_SEED + 100 + match(target, CONTRASTS))
  n_target <- sum(d$group == target)
  n_masl <- sum(d$group == "MASL")
  lower_critical <- qwilcox(0.025, n_target, n_masl)
  upper_critical <- qwilcox(0.975, n_target, n_masl)
  alpha_actual <- pwilcox(lower_critical, n_target, n_masl) +
    (1 - pwilcox(upper_critical - 1, n_target, n_masl))
  target_idx <- seq_len(n_target)
  pooled_log2 <- c(
    d$log2cpm[d$group == target],
    d$log2cpm[d$group == "MASL"]
  )
  pooled_centered <- pooled_log2 - median(pooled_log2)
  n_pair <- n_target + n_masl
  simulated_rank_values <- matrix(
    sample(pooled_centered, n_pair * N_MDE, replace = TRUE),
    nrow = n_pair,
    ncol = N_MDE
  )
  rows <- vector("list", length(MDE_GRID))
  for (i in seq_along(MDE_GRID)) {
    beta <- MDE_GRID[i]
    y_sim <- simulated_rank_values
    y_sim[seq_len(n_target), ] <- y_sim[seq_len(n_target), , drop = FALSE] + beta
    reject <- vapply(
      seq_len(N_MDE),
      function(j) {
        values <- y_sim[, j]
        ranked <- rank(values, ties.method = "average")
        u_target <- sum(ranked[target_idx]) - n_target * (n_target + 1) / 2
        u_target <= lower_critical || u_target >= upper_critical
      },
      logical(1)
    )
    rows[[i]] <- data.frame(
      target_group = target,
      contrast = paste(target, "vs MASL"),
      method = "exact_Wilcoxon_null_distribution",
      beta_log2 = beta,
      fold_change = 2^beta,
      power = mean(reject),
      simulations = N_MDE,
      seed = MDE_SEED + 100 + match(target, CONTRASTS),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  attr(out, "lower_critical") <- lower_critical
  attr(out, "upper_critical") <- upper_critical
  attr(out, "alpha_actual") <- alpha_actual
  out
}

model_power_curves <- do.call(
  rbind,
  lapply(CONTRASTS, simulate_permutation_power)
)
wilcoxon_curves <- do.call(
  rbind,
  lapply(CONTRASTS, simulate_wilcoxon_power)
)

mde_from_curve <- function(curve, method) {
  do.call(
    rbind,
    lapply(CONTRASTS, function(target) {
      z <- curve[curve$target_group == target, , drop = FALSE]
      hit <- z[is.finite(z$power) & z$power >= 0.80, , drop = FALSE]
      if (!nrow(hit)) {
        return(data.frame(
          contrast = paste(target, "vs MASL"),
          target_group = target,
          method = method,
          mde80_log2 = NA_real_,
          mde80_fold_change = NA_real_,
          simulations = N_MDE,
          stringsAsFactors = FALSE
        ))
      }
      minimum_beta <- min(hit$beta_log2, na.rm = TRUE)
      data.frame(
        contrast = paste(target, "vs MASL"),
        target_group = target,
        method = method,
        mde80_log2 = minimum_beta,
        mde80_fold_change = 2^minimum_beta,
        simulations = N_MDE,
        stringsAsFactors = FALSE
      )
    })
  )
}

mde <- rbind(
  mde_from_curve(model_power_curves, "permutation_calibrated_model_design"),
  mde_from_curve(wilcoxon_curves, "exact_Wilcoxon_null_distribution")
)

summary_rows <- do.call(
  rbind,
  lapply(CONTRASTS, function(target) {
    p <- primary[primary$target_group == target, ]
    pm <- permutation[permutation$target_group == target, ]
    lo <- loo[loo$target_group == target, ]
    wf <- full_wilcoxon[full_wilcoxon$target_group == target, ]
    wl <- loo_wilcoxon[loo_wilcoxon$target_group == target, ]
    mm <- mde[
      mde$target_group == target &
        mde$method == "permutation_calibrated_model_design",
    ]
    mw <- mde[
      mde$target_group == target &
        mde$method == "exact_Wilcoxon_null_distribution",
    ]
    model_concordant <- (
      p$bh_q < 0.05 &&
        pm$bh_q < 0.05 &&
        lo$bh_q < 0.05 &&
        all(c(p$beta_log2, pm$observed_beta_log2, lo$beta_log2) > 0)
    )
    rank_concordant <- (
      wf$six_pair_bh_q < 0.05 &&
        wl$six_pair_bh_q < 0.05 &&
        wf$statistic > qwilcox(
          0.5,
          sum(d$group == target),
          sum(d$group == "MASL")
        ) &&
        wl$statistic > qwilcox(
          0.5,
          sum(d_loo$group == target),
          sum(d_loo$group == "MASL")
        )
    )
    classification <- if (model_concordant && rank_concordant) {
      "concordant_across_checks"
    } else if (model_concordant) {
      "model_supported_rank_inconclusive"
    } else {
      "model_dependent_or_unresolved_precision"
    }
    data.frame(
      contrast = paste(target, "vs MASL"),
      target_group = target,
      masl_n = sum(d$group == "MASL"),
      target_n = sum(d$group == target),
      primary_beta_log2 = p$beta_log2,
      primary_fold_change = p$fold_change,
      primary_p = p$raw_p,
      primary_bh_q = p$bh_q,
      permutation_p = pm$permutation_p,
      permutation_bh_q = pm$bh_q,
      loo_model_beta_log2 = lo$beta_log2,
      loo_model_fold_change = lo$fold_change,
      loo_model_p = lo$raw_p,
      loo_model_bh_q = lo$bh_q,
      full_exact_wilcoxon_p = wf$p,
      full_exact_wilcoxon_direct_pair_bh_q = wf$direct_pair_bh_q,
      full_exact_wilcoxon_six_pair_bh_q = wf$six_pair_bh_q,
      loo_exact_wilcoxon_p = wl$p,
      loo_exact_wilcoxon_direct_pair_bh_q = wl$direct_pair_bh_q,
      loo_exact_wilcoxon_six_pair_bh_q = wl$six_pair_bh_q,
      model_mde80_log2 = mm$mde80_log2,
      model_mde80_fold_change = mm$mde80_fold_change,
      wilcoxon_mde80_log2 = mw$mde80_log2,
      wilcoxon_mde80_fold_change = mw$mde80_fold_change,
      model_concordant = model_concordant,
      rank_concordant = rank_concordant,
      classification = classification,
      stringsAsFactors = FALSE
    )
  })
)

input_checks <- data.frame(
  item = c(
    "pseudobulk_file",
    "sample_depth_file",
    "n_total",
    "n_ctrl",
    "n_MASL",
    "n_eMASH",
    "n_aMASH",
    "FL511_cpm10k",
    "FL511_rank_within_MASL",
    "MASL_median_cpm10k",
    "eMASH_median_cpm10k",
    "aMASH_median_cpm10k"
  ),
  value = c(
    PSEUDO_PATH,
    SAMPLE_PATH,
    nrow(d),
    sum(d$group == "ctrl"),
    sum(d$group == "MASL"),
    sum(d$group == "eMASH"),
    sum(d$group == "aMASH"),
    d$cpm10k[d$sample == "FL511"],
    rank(-d$cpm10k[d$sample == "FL511"], ties.method = "min"),
    median(d$cpm10k[d$group == "MASL"]),
    median(d$cpm10k[d$group == "eMASH"]),
    median(d$cpm10k[d$group == "aMASH"])
  ),
  stringsAsFactors = FALSE
)

permutation_detail <- do.call(
  rbind,
  lapply(CONTRASTS, function(target) {
    z <- nulls[[target]]
    data.frame(
      target_group = target,
      permutations = B_PERM,
      seed = PERM_SEED + match(target, CONTRASTS),
      observed_beta_log2 = z$observed,
      null_center = z$null_center,
      null_q025 = as.numeric(quantile(z$perm, 0.025)),
      null_q975 = as.numeric(quantile(z$perm, 0.975)),
      null_abs_q95 = z$null_abs_q95,
      p_two_sided = z$permutation_p,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(input_checks, file.path(OUT, "input_checks.csv"), row.names = FALSE)
write.csv(primary, file.path(OUT, "primary_reproduction.csv"), row.names = FALSE)
write.csv(full_wilcoxon, file.path(OUT, "exact_wilcoxon_full_data.csv"), row.names = FALSE)
write.csv(loo, file.path(OUT, "FL511_leave_one_donor_out_model.csv"), row.names = FALSE)
write.csv(
  loo_wilcoxon,
  file.path(OUT, "FL511_leave_one_donor_out_wilcoxon.csv"),
  row.names = FALSE
)
write.csv(
  permutation,
  file.path(OUT, "coefficient_permutation_results.csv"),
  row.names = FALSE
)
write.csv(
  permutation_detail,
  file.path(OUT, "coefficient_permutation_null.csv"),
  row.names = FALSE
)
write.csv(model_power_curves, file.path(OUT, "model_MDE_power_curve.csv"), row.names = FALSE)
write.csv(
  wilcoxon_curves,
  file.path(OUT, "exact_wilcoxon_MDE_power_curve.csv"),
  row.names = FALSE
)
write.csv(mde, file.path(OUT, "design_MDE_summary.csv"), row.names = FALSE)
write.csv(summary_rows, file.path(OUT, "MASL_MASH_sensitivity_summary.csv"), row.names = FALSE)

report_lines <- c(
  "# MASL-MASH Core Comparison Sensitivity Analysis",
  "",
  "## Analysis specification",
  "",
  "- Primary model: donor-level `log2(CP10K + 0.01) ~ group + log(median nCount)`.",
  "- Reference group: MASL. Direct contrasts: eMASH vs MASL and aMASH vs MASL.",
  "- Multiplicity: BH across the two direct MASL-reference contrasts.",
  "- Permutation: target-specific Freedman-Lane residual permutation, 20,000 permutations.",
  "- FL511 rule: one leave-one-donor-out analysis with all other donors retained.",
  "- Design MDE: 80% power at unadjusted two-sided alpha 0.05; reported in log2 expression units and fold-change units.",
  "",
  "## Input",
  "",
  sprintf(
    "- %d donors: ctrl %d, MASL %d, eMASH %d, aMASH %d.",
    nrow(d),
    sum(d$group == "ctrl"),
    sum(d$group == "MASL"),
    sum(d$group == "eMASH"),
    sum(d$group == "aMASH")
  ),
  sprintf(
    "- FL511 CP10K = %.6f, the highest of the six MASL donors.",
    d$cpm10k[d$sample == "FL511"]
  ),
  "",
  "## Results",
  "",
  paste0(
    "- ",
    summary_rows$contrast,
    ": primary fold change ",
    sprintf("%.2f", summary_rows$primary_fold_change),
    " (BH q = ",
    format.pval(summary_rows$primary_bh_q, digits = 3, eps = 1e-4),
    "); permutation P = ",
    format.pval(summary_rows$permutation_p, digits = 3, eps = 1e-4),
    " (BH q = ",
    format.pval(summary_rows$permutation_bh_q, digits = 3, eps = 1e-4),
    "); after removing FL511, model fold change ",
    sprintf("%.2f", summary_rows$loo_model_fold_change),
    " (BH q = ",
    format.pval(summary_rows$loo_model_bh_q, digits = 3, eps = 1e-4),
    "); exact Wilcoxon P = ",
    format.pval(summary_rows$full_exact_wilcoxon_p, digits = 3, eps = 1e-4),
    " (BH q = ",
    format.pval(
      summary_rows$full_exact_wilcoxon_direct_pair_bh_q,
      digits = 3,
      eps = 1e-4
    ),
    "; six-comparison rank-family q = ",
    format.pval(
      summary_rows$full_exact_wilcoxon_six_pair_bh_q,
      digits = 3,
      eps = 1e-4
    ),
    "); after removing FL511, exact Wilcoxon P = ",
    format.pval(summary_rows$loo_exact_wilcoxon_p, digits = 3, eps = 1e-4),
    " (BH q = ",
    format.pval(
      summary_rows$loo_exact_wilcoxon_direct_pair_bh_q,
      digits = 3,
      eps = 1e-4
    ),
    "; six-comparison rank-family q = ",
    format.pval(
      summary_rows$loo_exact_wilcoxon_six_pair_bh_q,
      digits = 3,
      eps = 1e-4
    ),
    "); model MDE80 = ",
    sprintf("%.2f", summary_rows$model_mde80_fold_change),
    "-fold; Wilcoxon MDE80 = ",
    sprintf("%.2f", summary_rows$wilcoxon_mde80_fold_change),
    "-fold."
  ),
  "",
  "## Interpretation",
  "",
  paste0(
    "- ",
    summary_rows$contrast,
    ": ",
    summary_rows$classification,
    "."
  ),
  "",
  "The primary model is unchanged. The Wilcoxon result is a rank-based",
  "sensitivity check, not a replacement analysis. P values above 0.05 are",
  "interpreted as failure to establish a difference, not evidence of",
  "equivalence. No equivalence boundary is introduced.",
  ""
)
writeLines(report_lines, file.path(OUT, "MASL_MASH_sensitivity_report.md"), useBytes = TRUE)

cat("MASL_MASH_SENSITIVITY_COMPLETE\n")
cat("OUTPUT_DIR=", OUT, "\n", sep = "")
