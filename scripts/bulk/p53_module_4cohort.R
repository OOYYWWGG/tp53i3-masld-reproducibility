suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(metafor)
})

P53 <- c(
  "CDKN1A", "MDM2", "BBC3", "GADD45A", "RRM2B", "SESN1", "ZMAT3",
  "TIGAR", "FDXR", "PHLDA3", "TP53INP1", "XPC", "DDB2"
)

score_file <- "INPUT_DIR/bulk_module_scores.csv"
scores <- read.csv(
  score_file, check.names = FALSE, stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

cohort_rows <- list()
for (ds in c("GSE135251", "GSE48452", "GSE89632")) {
  d <- scores[
    scores$dataset == ds &
      scores$group_unified %in% c("NAFL", "NASH") &
      !is.na(scores$fibrosis_stage_numeric) &
      scores$fibrosis_stage_numeric <= 3,
  ]
  ct <- cor.test(
    d$p53_ssgsea_z, d$fibrosis_stage_numeric,
    method = "spearman", exact = FALSE
  )
  cohort_rows[[ds]] <- data.frame(
    dataset = ds, n = nrow(d), rho = unname(ct$estimate),
    p = ct$p.value, stringsAsFactors = FALSE
  )
}

bulk <- "INPUT_DIR/composition_inputs/gse281797"
lines <- readLines(gzfile(file.path(bulk, "GSE281797_series_matrix.txt.gz")))
gv <- function(prefix) {
  z <- grep(paste0("^", prefix), lines, value = TRUE)
  gsub('^"|"$', "", strsplit(z[1], "\t", fixed = TRUE)[[1]][-1])
}
gc_ <- function(pattern) {
  z <- grep(pattern, lines, value = TRUE)
  v <- gsub('^"|"$', "", strsplit(z[1], "\t", fixed = TRUE)[[1]][-1])
  sub("^[^:]+: ", "", v)
}
meta <- data.frame(
  sample_id = sub("^HumanLiver_", "", gv("!Sample_title")),
  fibrosis = as.integer(gc_("fibrosis grade:")),
  diagnosis = gc_("diagnosis:"),
  stringsAsFactors = FALSE
)
counts <- fread(file.path(bulk, "GSE281797_GeneCount.txt.gz"), data.table = FALSE)
ids <- sub("\\..*$", "", as.character(counts[[1]]))
m <- as.matrix(counts[, -1, drop = FALSE])
rownames(m) <- ids
tpm <- fread(file.path(bulk, "GSE281797_tpm.txt.gz"), data.table = FALSE)[, 1:2]
colnames(tpm) <- c("ensembl", "symbol")
tpm$ensembl <- sub("\\..*$", "", tpm$ensembl)
map <- tpm[!duplicated(tpm$ensembl), ]
sym <- map$symbol[match(rownames(m), tpm$ensembl)]
ok <- !is.na(sym) & sym != ""
m <- rowsum(m[ok, , drop = FALSE], group = sym[ok], reorder = FALSE)
cpm <- sweep(m, 2, colSums(m), "/") * 1e6
logcpm <- log2(cpm + 1)
meta <- meta[match(colnames(logcpm), meta$sample_id), ]
meta <- meta[
  meta$diagnosis %in% c("MASL", "MASH") &
    !is.na(meta$fibrosis) & meta$fibrosis <= 3,
]
mat <- logcpm[, meta$sample_id, drop = FALSE]
sc <- gsva(
  ssgseaParam(mat, list(p53 = P53), normalize = TRUE, alpha = 0.25),
  verbose = FALSE
)
p53_score <- as.numeric(sc["p53", ])
ct <- cor.test(
  p53_score, meta$fibrosis, method = "spearman", exact = FALSE
)
cohort_rows[["GSE281797"]] <- data.frame(
  dataset = "GSE281797", n = nrow(meta), rho = unname(ct$estimate),
  p = ct$p.value, stringsAsFactors = FALSE
)

cohort <- do.call(rbind, cohort_rows)
mfit <- rma(
  yi = atanh(cohort$rho),
  vi = 1.06 / (cohort$n - 3),
  method = "REML",
  test = "knha"
)
res <- data.frame(
  k = mfit$k,
  pooled_rho = tanh(as.numeric(mfit$b)),
  ci_low = tanh(mfit$ci.lb),
  ci_high = tanh(mfit$ci.ub),
  p = mfit$pval,
  tau2 = mfit$tau2,
  I2 = mfit$I2,
  cohorts = paste(cohort$dataset, collapse = ",")
)

write.csv(cohort, "INPUT_DIR/p53_module_4cohort_cohorts.csv", row.names = FALSE)
write.csv(res, "INPUT_DIR/p53_module_4cohort_meta.csv", row.names = FALSE)
