## 5.3：MuSiC 估计与「经典 marker 打分」的交叉验证
##  marker：12 个肝细胞经典基因，队列内 z 后取均值
suppressPackageStartupMessages({ library(parallel) })

MARKERS <- c("ALB","APOA1","APOB","TTR","TF","SERPINA1","HNF4A","FGA","FGB","FGG","HP","APOC3")

dec <- read.csv("deconvolution_full_proportions.csv", check.names = FALSE,
                stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
names(dec)[1] <- sub("^\ufeff", "", names(dec)[1])
tp <- read.csv("TP53I3_analysis_dataset.csv", check.names = FALSE,
               stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
names(tp)[1] <- sub("^\ufeff", "", names(tp)[1])

FILES <- c(GSE135251 = "GSE135251_protein_coding.csv",
           GSE48452  = "GSE48452_protein_coding.csv",
           GSE89632  = "GSE89632_protein_coding.csv",
           GSE49541  = "GSE49541_protein_coding.csv")

out <- list()
for (ds in names(FILES)) {
  d <- read.csv(FILES[[ds]], check.names = FALSE, stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM")
  names(d)[1] <- sub("^\ufeff", "", names(d)[1])
  genes <- d[[1]]
  m <- as.matrix(d[, -1, drop = FALSE]); storage.mode(m) <- "double"
  m <- rowsum(m, group = genes)
  hit <- intersect(MARKERS, rownames(m))
  z <- t(scale(t(m[hit, , drop = FALSE])))
  score <- colMeans(z, na.rm = TRUE)

  dd <- dec[dec$dataset == ds, ]
  keep <- intersect(names(score), dd$sample)
  s <- data.frame(sample = keep, hep_score = score[keep])
  s <- merge(s, dd[, c("sample", "Hepatocyte", "Cholangiocyte")], by = "sample")
  s <- merge(s, tp[tp$dataset == ds, c("gsm", "tp53i3")], by.x = "sample", by.y = "gsm")

  sp <- function(x, y) {
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    c(rho = unname(ct$estimate), p = ct$p.value)
  }
  a <- sp(s$hep_score, s$Hepatocyte)
  b <- sp(s$hep_score, s$Cholangiocyte)
  c3 <- sp(s$hep_score, s$tp53i3)
  out[[ds]] <- data.frame(
    dataset = ds, n = nrow(s), markers = length(hit),
    hep_mean = mean(s$Hepatocyte), hep_min = min(s$Hepatocyte), hep_max = max(s$Hepatocyte),
    rho_music = a["rho"], p_music = a["p"],
    rho_chol = b["rho"], p_chol = b["p"],
    rho_tp53i3 = c3["rho"], p_tp53i3 = c3["p"],
    stringsAsFactors = FALSE)
  cat(sprintf("%-10s n=%3d markers=%2d | hep fraction %.3f (%.3f-%.3f) | rho(score, MuSiC)=%+.3f p=%.2g | rho(score, TP53I3)=%+.3f p=%.2g\n",
              ds, nrow(s), length(hit), mean(s$Hepatocyte), min(s$Hepatocyte), max(s$Hepatocyte),
              a["rho"], a["p"], c3["rho"], c3["p"]))
}

res <- do.call(rbind, out)
write.csv(res, "marker_crosscheck.csv", row.names = FALSE)
cat("\nwritten marker_crosscheck.csv\n")
