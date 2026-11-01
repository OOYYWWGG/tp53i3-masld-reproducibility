sl <- read.csv(
  "INPUT_DIR/sc_nagq/sample_level_hepatocyte.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

sl$group <- factor(sl$group, levels = c("MASL", "ctrl", "eMASH", "aMASH"))
sl$logit_raw <- qlogis(sl$n_positive / sl$n_cells)
sl$logit_cc <- qlogis((sl$n_positive + 0.5) / (sl$n_cells + 1))
sl$log_ncount <- log(sl$median_nCount)

fit_contrast <- function(outcome, label) {
  model <- glm(
    reformulate(c("group", "log_ncount"), response = outcome),
    data = sl
  )
  coefs <- summary(model)$coefficients
  ci <- confint.default(model)
  rows <- lapply(grep("^group", rownames(coefs), value = TRUE), function(term) {
    data.frame(
      model = label,
      contrast = sub("^group", "", term),
      log_or = coefs[term, 1],
      se = coefs[term, 2],
      or = exp(coefs[term, 1]),
      ci_low = exp(ci[term, 1]),
      ci_high = exp(ci[term, 2]),
      p = coefs[term, 4],
      mde80_log_or = 2.802 * coefs[term, 2],
      mde80_or = exp(2.802 * coefs[term, 2]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

res <- rbind(
  fit_contrast("logit_raw", "primary_raw_logit"),
  fit_contrast("logit_cc", "sensitivity_continuity_corrected")
)

pairwise <- pairwise.wilcox.test(
  sl$positive_pct,
  sl$group,
  p.adjust.method = "BH",
  exact = FALSE
)
pairwise_df <- as.data.frame(as.table(pairwise$p.value))
names(pairwise_df) <- c("group_a", "group_b", "bh_p")

write.csv(
  res,
  "INPUT_DIR/masl_contrast_results.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  pairwise_df,
  "INPUT_DIR/masl_contrast_wilcoxon.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
