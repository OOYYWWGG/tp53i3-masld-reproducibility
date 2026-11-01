suppressPackageStartupMessages({
  library(metafor)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else Sys.getenv("TEMP")
out_dir <- file.path(base_dir, "module_score_output")

d <- read.csv(file.path(out_dir, "bulk_module_scores.csv"),
              check.names = FALSE, stringsAsFactors = FALSE)

for (ds in unique(d$dataset)) {
  i <- d$dataset == ds
  d$p53_mean_std[i] <- as.numeric(scale(d$p53_mean_z[i]))
  d$zonation_mean_std[i] <- as.numeric(scale(d$zonation_mean_z[i]))
  d$state_mean_std[i] <- as.numeric(scale(d$state_mean_z[i]))
  d$ductular_z[i] <- as.numeric(scale(d$ductular_score[i]))
}

sp_ci <- function(x, y) {
  ok <- complete.cases(x, y)
  n <- sum(ok)
  if (n < 8 || length(unique(y[ok])) < 2) {
    return(c(n = n, rho = NA, lcl = NA, ucl = NA, p = NA))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman",
                                  exact = FALSE))
  z <- atanh(unname(ct$estimate))
  se <- 1 / sqrt(n - 3)
  c(n = n, rho = unname(ct$estimate),
    lcl = tanh(z - qnorm(0.975) * se),
    ucl = tanh(z + qnorm(0.975) * se),
    p = ct$p.value)
}

partial_sp <- function(x, y, z) {
  ok <- complete.cases(x, y, z)
  n <- sum(ok)
  if (n < 8 || length(unique(y[ok])) < 2) {
    return(c(n = n, rho = NA, p = NA))
  }
  rx <- residuals(lm(rank(x[ok]) ~ rank(z[ok])))
  ry <- residuals(lm(rank(y[ok]) ~ rank(z[ok])))
  ct <- suppressWarnings(cor.test(rx, ry, method = "pearson"))
  c(n = n, rho = unname(ct$estimate), p = ct$p.value)
}

score_map <- data.frame(
  score = c("p53_ssgsea_z", "zonation_ssgsea_z", "state_ssgsea_z",
            "exploratory_ssgsea_z", "p53_mean_std", "zonation_mean_std",
            "state_mean_std"),
  score_type = c(rep("ssGSEA", 4), rep("mean_z", 3)),
  component = c("p53", "zonation", "state", "exploratory",
                "p53", "zonation", "state"),
  stringsAsFactors = FALSE
)

cor_rows <- list()
for (ds in unique(d$dataset)) {
  x <- d[d$dataset == ds, ]
  endpoints <- c("tp53i3", "nas")
  if (ds == "GSE49541") {
    endpoints <- c(endpoints, "fib_bin")
  } else {
    endpoints <- c(endpoints, "fib")
  }
  for (k in seq_len(nrow(score_map))) {
    s <- score_map$score[k]
    for (ep in endpoints) {
      a <- sp_ci(x[[s]], x[[ep]])
      b <- partial_sp(x[[s]], x[[ep]], x$hepatocyte_fraction)
      cor_rows[[length(cor_rows) + 1]] <- data.frame(
        dataset = ds,
        score = s,
        score_type = score_map$score_type[k],
        component = score_map$component[k],
        endpoint = ep,
        n = a["n"],
        rho = a["rho"],
        rho_lcl = a["lcl"],
        rho_ucl = a["ucl"],
        p = a["p"],
        partial_n = b["n"],
        partial_rho_hep = b["rho"],
        partial_p_hep = b["p"],
        stringsAsFactors = FALSE
      )
    }
  }
}
cor_tab <- do.call(rbind, cor_rows)
write.csv(cor_tab, file.path(out_dir, "module_correlations.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Composition attenuation using standardized linear models.
att_rows <- list()
for (ds in unique(d$dataset)) {
  x <- d[d$dataset == ds, ]
  endpoints <- if (ds == "GSE49541") "fib_bin" else c("fib", "nas")
  for (k in seq_len(nrow(score_map))) {
    s <- score_map$score[k]
    for (ep in endpoints) {
      x$endpoint_z <- as.numeric(scale(x[[ep]]))
      f0 <- as.formula(paste(s, "~ endpoint_z"))
      f1 <- as.formula(paste(s, "~ endpoint_z + hepatocyte_fraction"))
      f2 <- as.formula(paste(s, "~ endpoint_z + hepatocyte_fraction + ductular_z"))
      m0 <- try(lm(f0, data = x), silent = TRUE)
      m1 <- try(lm(f1, data = x), silent = TRUE)
      m2 <- try(lm(f2, data = x), silent = TRUE)
      if (any(vapply(list(m0, m1, m2), inherits, logical(1), "try-error"))) next
      get_term <- function(m, term) {
        co <- summary(m)$coefficients
        if (!(term %in% rownames(co))) return(c(beta = NA, se = NA, p = NA))
        c(beta = co[term, "Estimate"], se = co[term, "Std. Error"],
          p = co[term, "Pr(>|t|)"])
      }
      b0 <- get_term(m0, "endpoint_z")
      b1 <- get_term(m1, "endpoint_z")
      b2 <- get_term(m2, "endpoint_z")
      h1 <- get_term(m1, "hepatocyte_fraction")
      h2 <- get_term(m2, "hepatocyte_fraction")
      ci <- function(b, se) c(b - qnorm(0.975) * se, b + qnorm(0.975) * se)
      att_rows[[length(att_rows) + 1]] <- data.frame(
        dataset = ds,
        score = s,
        score_type = score_map$score_type[k],
        component = score_map$component[k],
        endpoint = ep,
        n = nobs(m0),
        beta_base = b0["beta"],
        beta_base_lcl = ci(b0["beta"], b0["se"])[1],
        beta_base_ucl = ci(b0["beta"], b0["se"])[2],
        p_base = b0["p"],
        beta_hep = b1["beta"],
        beta_hep_lcl = ci(b1["beta"], b1["se"])[1],
        beta_hep_ucl = ci(b1["beta"], b1["se"])[2],
        p_hep = b1["p"],
        beta_hep_duct = b2["beta"],
        beta_hep_duct_lcl = ci(b2["beta"], b2["se"])[1],
        beta_hep_duct_ucl = ci(b2["beta"], b2["se"])[2],
        p_hep_duct = b2["p"],
        hepatocyte_beta = h2["beta"],
        hepatocyte_p = h2["p"],
        attenuation_hep_pct = 100 * (b0["beta"] - b1["beta"]) / b0["beta"],
        attenuation_full_pct = 100 * (b0["beta"] - b2["beta"]) / b0["beta"],
        stringsAsFactors = FALSE
      )
    }
  }
}
att_tab <- do.call(rbind, att_rows)
write.csv(att_tab, file.path(out_dir, "composition_attenuation.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Random-effects meta of Spearman rho (Fisher z).
meta_rho <- function(rho, n) {
  yi <- atanh(rho)
  vi <- 1 / (n - 3)
  fit <- metafor::rma(yi, vi, method = "REML")
  c(k = length(rho), pooled_rho = tanh(as.numeric(coef(fit))),
    lcl = tanh(fit$ci.lb), ucl = tanh(fit$ci.ub),
    p = fit$pval, I2 = fit$I2)
}

meta_rows <- list()
for (ep in c("fib", "nas")) {
  for (k in seq_len(nrow(score_map))) {
    s <- score_map$score[k]
    sub <- cor_tab[cor_tab$endpoint == ep & cor_tab$score == s &
                     cor_tab$dataset %in% c("GSE135251", "GSE48452", "GSE89632"), ]
    sub <- sub[is.finite(sub$rho) & sub$n > 3, ]
    if (nrow(sub) < 2) next
    m <- meta_rho(sub$rho, sub$n)
    meta_rows[[length(meta_rows) + 1]] <- data.frame(
      endpoint = ep,
      score = s,
      score_type = score_map$score_type[k],
      component = score_map$component[k],
      k = m["k"], pooled_rho = m["pooled_rho"],
      lcl = m["lcl"], ucl = m["ucl"], p = m["p"], I2 = m["I2"],
      cohorts = paste(sub$dataset, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
}
meta_tab <- do.call(rbind, meta_rows)
write.csv(meta_tab, file.path(out_dir, "meta_correlations.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
