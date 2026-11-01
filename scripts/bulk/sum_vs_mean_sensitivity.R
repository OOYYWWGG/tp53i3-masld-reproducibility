## ============================================================================
## 合并方式敏感性：sum（原方案） vs mean
## Run in a writable working directory containing the listed annotation tables
## and TP53I3_analysis_dataset.csv. All paths below are relative.
## ============================================================================
suppressPackageStartupMessages({
  library(GSVA)
})

P53 <- c("CDKN1A","MDM2","BBC3","GADD45A","RRM2B","SESN1","ZMAT3","TIGAR","FDXR","PHLDA3","TP53INP1","XPC","DDB2")
SETS <- c(GSE135251 = "GSE135251_全基因注释表.csv.gz",
          GSE48452  = "GSE48452_全基因注释表.csv.gz",
          GSE89632  = "GSE89632_全基因注释表.csv.gz")

clin <- read.csv("TP53I3_analysis_dataset.csv", check.names = FALSE,
                 stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
names(clin)[1] <- sub("^\ufeff", "", names(clin)[1])

out <- list()
for (ds in names(SETS)) {
  d <- read.csv(gzfile(SETS[[ds]]), check.names = FALSE, stringsAsFactors = FALSE)
  ## 与 4.2 打分脚本一致：不按 gene_type 过滤，只要 gene_symbol 非空
  d <- d[!is.na(d$gene_name) & d$gene_name != "", ]
  samples <- setdiff(colnames(d), c("input_id", "gene_name", "gene_type"))
  m <- as.matrix(d[, samples]); storage.mode(m) <- "double"
  g <- d$gene_name

  m_sum <- rowsum(m, group = g)
  cnt <- as.numeric(table(g)[rownames(m_sum)])
  m_mean <- m_sum / cnt

  cl <- clin[clin$dataset == ds & clin$group_unified %in% c("NAFL", "NASH"), ]
  keep <- intersect(samples, cl$gsm)
  cl <- cl[match(keep, cl$gsm), ]

  gs <- list(p53 = P53)
  run <- function(mat, tag) {
    sc <- gsva(ssgseaParam(mat[, keep, drop = FALSE], gs, normalize = TRUE, alpha = 0.25),
               verbose = FALSE)
    v <- as.numeric(sc["p53", ])
    data.frame(dataset = ds, merge = tag, n = length(keep),
               rho_fib = cor(v, cl$fibrosis_stage_numeric, method = "spearman", use = "complete.obs"),
               rho_nas = cor(v, cl$nas, method = "spearman", use = "complete.obs"),
               max_abs_diff_score = max(abs(v - as.numeric(
                 gsva(ssgseaParam(mat[, keep, drop = FALSE], gs, normalize = TRUE, alpha = 0.25),
                      verbose = FALSE)["p53", ]))),
               stringsAsFactors = FALSE)
  }
  tmp_sum <- run(m_sum, "sum")
  tmp_mean <- run(m_mean, "mean")
  tmp_sum$max_abs_diff_score <- NA
  tmp_mean$max_abs_diff_score <- NA
  s <- as.numeric(gsva(ssgseaParam(m_sum[, keep, drop = FALSE], gs, normalize = TRUE, alpha = 0.25), verbose = FALSE)["p53", ])
  mn <- as.numeric(gsva(ssgseaParam(m_mean[, keep, drop = FALSE], gs, normalize = TRUE, alpha = 0.25), verbose = FALSE)["p53", ])
  tmp_sum$max_abs_diff_score <- max(abs(s - mn))
  tmp_mean$max_abs_diff_score <- max(abs(s - mn))
  tmp_sum$score_cor <- cor(s, mn, method = "spearman")
  tmp_mean$score_cor <- cor(s, mn, method = "spearman")
  out[[ds]] <- rbind(tmp_sum, tmp_mean)
  cat(ds, "done; duplicated genes:", sum(cnt > 1), "\n")
}

res <- do.call(rbind, out)
write.csv(res, "sum_vs_mean_sensitivity.csv", row.names = FALSE)

## Intermediate, runtime-only: three-cohort Fisher-z meta summary. The manuscript
## reports only the per-dataset sum/mean score_cor values in
## sum_vs_mean_sensitivity.csv; this meta table is not a reported result and is
## written at runtime but not frozen in the package.
fz <- function(r) atanh(r)
meta <- aggregate(cbind(rho_fib, rho_nas) ~ merge, res, function(x) tanh(mean(fz(x))))
dir.create("results", showWarnings = FALSE, recursive = TRUE)
write.csv(meta, "results/sum_vs_mean_sensitivity_meta.csv", row.names = FALSE)
