d <- read.csv(
  "INPUT_DIR/celltype_pseudobulk.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
d <- d[trimws(d[["ct"]]) == "Hepatocyte", ]

sl <- read.csv(
  "INPUT_DIR/sc_nagq/sample_level_hepatocyte.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
d <- merge(d, sl[, c("sample", "median_nCount")], by = "sample", all.x = TRUE)
d$group <- factor(d$group, levels = c("MASL", "ctrl", "eMASH", "aMASH"))
d$log2cpm <- log2(d$cpm10k + 0.01)

fit_model <- function(label) {
  m <- glm(log2cpm ~ group + log(median_nCount), data = d)
  co <- summary(m)$coefficients
  ci <- confint.default(m)
  rows <- lapply(grep("^group", rownames(co), value = TRUE), function(term) {
    data.frame(
      model = label,
      contrast = sub("^group", "", term),
      beta_log2 = co[term, 1],
      se = co[term, 2],
      ci_low = ci[term, 1],
      ci_high = ci[term, 2],
      p = co[term, 4],
      fold_change = 2^co[term, 1],
      fold_low = 2^ci[term, 1],
      fold_high = 2^ci[term, 2],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

res <- fit_model("MASL_reference_primary_pseudobulk")

d_ctrl <- d
d_ctrl$group <- factor(d_ctrl$group, levels = c("ctrl", "MASL", "eMASH", "aMASH"))
m_ctrl <- glm(log2cpm ~ group + log(median_nCount), data = d_ctrl)
co_ctrl <- summary(m_ctrl)$coefficients
ci_ctrl <- confint.default(m_ctrl)
ctrl_rows <- lapply(grep("^group", rownames(co_ctrl), value = TRUE), function(term) {
  data.frame(
    model = "ctrl_reference_primary_pseudobulk",
    contrast = sub("^group", "", term),
    beta_log2 = co_ctrl[term, 1],
    se = co_ctrl[term, 2],
    ci_low = ci_ctrl[term, 1],
    ci_high = ci_ctrl[term, 2],
    p = co_ctrl[term, 4],
    fold_change = 2^co_ctrl[term, 1],
    fold_low = 2^ci_ctrl[term, 1],
    fold_high = 2^ci_ctrl[term, 2],
    stringsAsFactors = FALSE
  )
})
res <- rbind(res, do.call(rbind, ctrl_rows))

pairwise <- pairwise.wilcox.test(
  d$cpm10k,
  d$group,
  p.adjust.method = "BH",
  exact = TRUE
)
pairwise_df <- as.data.frame(as.table(pairwise$p.value))
names(pairwise_df) <- c("group_a", "group_b", "bh_p")

write.csv(
  res,
  "INPUT_DIR/pseudobulk_primary_results.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  pairwise_df,
  "INPUT_DIR/pseudobulk_primary_wilcoxon_exact.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
