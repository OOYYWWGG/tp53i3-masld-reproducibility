suppressPackageStartupMessages(library(metafor))

dat <- data.frame(
  dataset = c("GSE135251", "GSE48452", "GSE89632"),
  n = c(206, 32, 39),
  rho = c(0.337309283912838, 0.5114015850822, 0.314520575569248)
)
dat$yi <- atanh(dat$rho)
dat$vi <- 1.06 / (dat$n - 3)

fit <- rma(yi, vi, data = dat, method = "REML", test = "knha")
pred <- predict(fit)
loo <- leave1out(fit)

res <- data.frame(
  analysis = "TP53I3_fibrosis_F0_F4_3cohorts",
  k = fit$k,
  pooled_rho = tanh(as.numeric(coef(fit))),
  ci_low = tanh(fit$ci.lb),
  ci_high = tanh(fit$ci.ub),
  p = fit$pval,
  tau2 = fit$tau2,
  I2 = fit$I2,
  pi_low = tanh(pred$pi.lb),
  pi_high = tanh(pred$pi.ub),
  method = paste(
    "Fisher z; variance = 1.06/(n - 3);",
    "REML; Hartung-Knapp; prediction interval"
  ),
  stringsAsFactors = FALSE
)

weights <- data.frame(
  dataset = dat$dataset,
  n = dat$n,
  rho = dat$rho,
  weight_percent = as.numeric(weights(fit)) / sum(weights(fit)) * 100
)

loo_tab <- data.frame(
  omitted = dat$dataset,
  pooled_rho = tanh(loo$estimate),
  ci_low = tanh(loo$ci.lb),
  ci_high = tanh(loo$ci.ub),
  p = loo$pval,
  tau2 = loo$tau2,
  I2 = loo$I2
)

result_dir <- Sys.getenv("TP53I3_RESULT_DIR", "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  res,
  file.path(result_dir, "tp53i3_meta_F0F4_3cohorts.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  weights,
  file.path(result_dir, "tp53i3_meta_F0F4_3cohort_weights.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  loo_tab,
  file.path(result_dir, "tp53i3_meta_F0F4_3cohort_LOO.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
