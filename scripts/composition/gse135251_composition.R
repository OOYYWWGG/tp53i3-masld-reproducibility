## ============================================================================
## GSE135251: TP53I3 ~ fibrosis/NAS + hepatocyte + ductular cells
## Input files: TP53I3_analysis_dataset.csv,
## deconvolution_full_proportions.csv, and ductular_score_bulk.csv
## under INPUT_DIR/composition_inputs.
## ============================================================================
suppressPackageStartupMessages(library(data.table))

BASE <- "INPUT_DIR/composition_inputs"
clin <- read.csv(file.path(BASE, "TP53I3_analysis_dataset.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
prop <- read.csv(file.path(BASE, "deconvolution_full_proportions.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
duct <- read.csv(file.path(BASE, "ductular_score_bulk.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")

d <- clin[clin$dataset == "GSE135251" & clin$group_unified %in% c("NAFL", "NASH"), ]
d$hep  <- prop$Hepatocyte[match(d$gsm, prop$sample)]
d$duct <- duct$ductular_score[match(d$gsm, duct$gsm)]

z <- function(x) as.numeric(scale(x))
d$y <- z(d$tp53i3_z); d$fib <- z(d$fibrosis_stage_numeric); d$nasz <- z(d$nas)
d$heps <- z(d$hep); d$ducts <- z(d$duct)

fit_show <- function(form, label) {
  m <- lm(form, data = d); co <- summary(m)$coefficients; ci <- confint(m)
  cat("\n==", label, "==\n")
  for (t in rownames(co)[-1])
    cat(sprintf("%-10s beta=%+.3f (%.3f, %.3f) P=%.3g\n", t, co[t,1], ci[t,1], ci[t,2], co[t,4]))
  cat("R2 =", round(summary(m)$r.squared, 3), " n =", nobs(m), "\n")
}
fit_show(y ~ fib, "fibrosis only")
fit_show(y ~ fib + heps, "fibrosis + hepatocyte")
fit_show(y ~ fib + heps + ducts, "fibrosis + hepatocyte + ductular")
fit_show(y ~ nasz, "NAS only")
fit_show(y ~ nasz + heps + ducts, "NAS + hepatocyte + ductular")
