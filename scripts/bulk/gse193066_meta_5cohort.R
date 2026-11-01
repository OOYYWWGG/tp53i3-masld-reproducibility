suppressPackageStartupMessages(library(metafor))

fib <- data.frame(
  dataset = c("GSE135251", "GSE48452", "GSE89632", "GSE281797", "GSE193066"),
  rho = c(0.3022, 0.4426, 0.3893, 0.4206, 0.2396),
  n = c(192, 30, 35, 64, 106)
)
nas <- data.frame(
  dataset = c("GSE135251", "GSE48452", "GSE89632", "GSE281797", "GSE193066"),
  rho = c(0.376, 0.319, 0.567, 0.336, 0.4614),
  n = c(206, 32, 38, 64, 106)
)

fit <- function(d, label) {
  m <- rma(
    yi = atanh(d$rho),
    vi = 1.06 / (d$n - 3),
    method = "REML",
    test = "knha"
  )
  data.frame(
    endpoint = label,
    k = m$k,
    pooled_rho = tanh(as.numeric(m$b)),
    ci_low = tanh(m$ci.lb),
    ci_high = tanh(m$ci.ub),
    p = m$pval,
    tau2 = m$tau2,
    I2 = m$I2,
    prediction_low = tanh(predict(m)$pi.lb),
    prediction_high = tanh(predict(m)$pi.ub),
    stringsAsFactors = FALSE
  )
}

res <- rbind(fit(fib, "fibrosis_F0_F3_5cohort"), fit(nas, "NAS_5cohort"))

m <- rma(yi = atanh(fib$rho), vi = 1.06 / (fib$n - 3), method = "REML", test = "knha")
weights_df <- data.frame(dataset = fib$dataset, weight_pct = as.numeric(weights(m)))
loo <- leave1out(m)
loo_df <- data.frame(
  dataset = fib$dataset,
  rho = tanh(loo$estimate),
  lcl = tanh(loo$ci.lb),
  ucl = tanh(loo$ci.ub)
)

write.csv(fib, "INPUT_DIR/gse193066_fib_inputs.csv", row.names = FALSE)
write.csv(res, "INPUT_DIR/gse193066_meta_results.csv", row.names = FALSE)
write.csv(weights_df, "INPUT_DIR/gse193066_meta_weights.csv", row.names = FALSE)
write.csv(loo_df, "INPUT_DIR/gse193066_meta_loo.csv", row.names = FALSE)
