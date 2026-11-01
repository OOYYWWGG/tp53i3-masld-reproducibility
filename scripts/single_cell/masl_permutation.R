suppressPackageStartupMessages(library(data.table))

pseudo <- read.csv(
  "INPUT_DIR/celltype_pseudobulk.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
pseudo <- pseudo[trimws(pseudo[["ct"]]) == "Hepatocyte", ]
sl <- read.csv(
  "INPUT_DIR/sc_nagq/sample_level_hepatocyte.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
d <- merge(
  pseudo[, c("sample", "group", "cpm10k")],
  sl[, c("sample", "median_nCount")],
  by = "sample"
)
d$group <- factor(d$group, levels = c("ctrl", "MASL", "eMASH", "aMASH"))
d$log_ncount <- log(d$median_nCount)

perm_p <- function(g1, g2, delta, B = 5000, seed = 20261020) {
  x <- d[d$group %in% c(g1, g2), ]
  x$y <- log2(x$cpm10k + delta)
  x$grp <- factor(ifelse(x$group == g1, 1, 0))
  fit <- lm(y ~ grp + log_ncount, data = x)
  obs <- unname(coef(fit)["grp1"])
  set.seed(seed)
  perm <- replicate(B, {
    xp <- x
    xp$grp <- sample(xp$grp)
    unname(coef(lm(y ~ grp + log_ncount, data = xp))["grp1"])
  })
  p_two <- (1 + sum(abs(perm) >= abs(obs))) / (B + 1)
  p_one <- (1 + sum(perm >= obs)) / (B + 1)
  data.frame(
    comparison = paste(g1, "vs", g2),
    delta = delta,
    n1 = sum(x$group == g1),
    n2 = sum(x$group == g2),
    log2_fold = obs,
    fold = 2^obs,
    p_two_sided = p_two,
    p_one_sided = p_one,
    stringsAsFactors = FALSE
  )
}

res <- rbind(
  perm_p("eMASH", "MASL", 0.001),
  perm_p("eMASH", "MASL", 0.01),
  perm_p("eMASH", "MASL", 0.1),
  perm_p("aMASH", "MASL", 0.001),
  perm_p("aMASH", "MASL", 0.01),
  perm_p("aMASH", "MASL", 0.1)
)

write.csv(
  res,
  "INPUT_DIR/masl_permutation_results.csv",
  row.names = FALSE
)
