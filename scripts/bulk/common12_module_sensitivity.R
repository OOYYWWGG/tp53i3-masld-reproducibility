suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(metafor)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

module_genes <- c(
  "CDKN1A", "MDM2", "BBC3", "GADD45A", "RRM2B", "SESN1",
  "ZMAT3", "TIGAR", "FDXR", "PHLDA3", "TP53INP1", "XPC", "DDB2"
)
common12 <- setdiff(module_genes, "TIGAR")
datasets <- c("GSE135251", "GSE48452", "GSE89632")

meta <- fread(
  file.path(base_dir, "bulk_module_scores.csv"),
  data.table = FALSE,
  check.names = FALSE
)

read_expr <- function(dataset) {
  x <- fread(
    file.path(base_dir, paste0(dataset, ".csv")),
    data.table = FALSE,
    check.names = FALSE
  )
  symbols <- as.character(x[[1]])
  matrix <- as.matrix(x[, -1, drop = FALSE])
  keep <- !is.na(symbols) & symbols != ""
  matrix <- matrix[keep, , drop = FALSE]
  symbols <- symbols[keep]
  matrix <- rowsum(matrix, group = symbols, reorder = FALSE)
  if (dataset == "GSE135251") {
    matrix <- sweep(matrix, 2, colSums(matrix), "/") * 1e6
    matrix <- log2(matrix + 1)
  }
  matrix
}

score_rows <- list()
for (dataset in datasets) {
  matrix <- read_expr(dataset)
  missing_genes <- setdiff(common12, rownames(matrix))
  if (length(missing_genes)) {
    stop(dataset, " missing common genes: ", paste(missing_genes, collapse = ","))
  }
  info <- meta[meta$dataset == dataset, , drop = FALSE]
  samples <- intersect(info$gsm, colnames(matrix))
  info <- info[match(samples, info$gsm), , drop = FALSE]
  matrix <- matrix[, samples, drop = FALSE]
  parameter <- ssgseaParam(
    matrix,
    list(P53_MODULE_COMMON12 = common12),
    normalize = TRUE,
    alpha = 0.25,
    verbose = FALSE
  )
  score <- as.numeric(gsva(parameter)["P53_MODULE_COMMON12", ])
  score_rows[[length(score_rows) + 1]] <- data.frame(
    dataset = dataset,
    gsm = samples,
    p53_common12_ssgsea = score,
    stringsAsFactors = FALSE
  )
}

scores <- do.call(rbind, score_rows)
scores$p53_common12_ssgsea_z <- ave(
  scores$p53_common12_ssgsea,
  scores$dataset,
  FUN = function(x) as.numeric(scale(x))
)
write.csv(
  scores,
  file.path(out_dir, "p53_common12_scores.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

merged <- merge(
  meta,
  scores[, c("dataset", "gsm", "p53_common12_ssgsea_z")],
  by = c("dataset", "gsm"),
  all.x = TRUE,
  sort = FALSE
)
write.csv(
  merged,
  file.path(out_dir, "bulk_module_scores_common12.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

spearman_row <- function(dataset, outcome) {
  x <- merged[merged$dataset == dataset, , drop = FALSE]
  ok <- is.finite(x$p53_common12_ssgsea_z) & is.finite(x[[outcome]])
  test <- suppressWarnings(cor.test(
    x$p53_common12_ssgsea_z[ok],
    x[[outcome]][ok],
    method = "spearman",
    exact = FALSE
  ))
  n <- sum(ok)
  z <- atanh(unname(test$estimate))
  se <- sqrt(1.06 / max(n - 3, 1))
  data.frame(
    dataset = dataset,
    outcome = outcome,
    n = n,
    rho = unname(test$estimate),
    p = test$p.value,
    lcl = tanh(z - 1.96 * se),
    ucl = tanh(z + 1.96 * se)
  )
}

correlations <- do.call(rbind, lapply(
  datasets,
  function(dataset) do.call(rbind, lapply(
    c("fib", "nas"),
    function(outcome) spearman_row(dataset, outcome)
  ))
))
write.csv(
  correlations,
  file.path(out_dir, "p53_common12_correlations.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

meta_rows <- list()
for (outcome in c("fib", "nas")) {
  x <- correlations[correlations$outcome == outcome, ]
  fit <- rma.uni(
    yi = atanh(x$rho),
    vi = 1.06 / (x$n - 3),
    method = "REML",
    test = "knha"
  )
  meta_rows[[length(meta_rows) + 1]] <- data.frame(
    outcome = outcome,
    k = nrow(x),
    pooled_rho = tanh(as.numeric(fit$b)),
    lcl = tanh(as.numeric(fit$ci.lb)),
    ucl = tanh(as.numeric(fit$ci.ub)),
    p = fit$pval,
    I2 = fit$I2
  )
}
meta_table <- do.call(rbind, meta_rows)
write.csv(
  meta_table,
  file.path(out_dir, "p53_common12_meta.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
