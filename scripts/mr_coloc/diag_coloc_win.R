suppressPackageStartupMessages(library(coloc))

data_dir <- "data/mr"
out_dir <- "results/mr_coloc"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

e <- read.csv(
  file.path(data_dir, "exposure_TP53I3_eQTLGen_cis.csv"),
  stringsAsFactors = FALSE
)
o <- read.csv(
  file.path(data_dir, "outcome_GCST90054782_NAFLD_4761cases.csv"),
  stringsAsFactors = FALSE
)

m <- merge(e, o, by = "SNP")
m$pos <- m$pos.outcome
m$mafe <- pmin(m$eaf.exposure, 1 - m$eaf.exposure)
m$mafo <- pmin(m$eaf.outcome, 1 - m$eaf.outcome)

i <- which.min(m$pval.exposure)
lead <- m$pos[i]

for (hw in c(250000, 500000)) {
  d <- m[abs(m$pos - lead) <= hw, , drop = FALSE]
  D1 <- list(
    snp = d$SNP,
    position = d$pos,
    beta = d$beta.exposure,
    varbeta = d$se.exposure^2,
    type = "quant",
    N = 31684,
    MAF = d$mafe
  )
  D2 <- list(
    snp = d$SNP,
    position = d$pos,
    beta = d$beta.outcome,
    varbeta = d$se.outcome^2,
    type = "cc",
    N = 4761 + 373227,
    s = 4761 / (4761 + 373227),
    MAF = d$mafo
  )
  r <- tryCatch(coloc.abf(D1, D2, p12 = 1e-5), error = function(err) NULL)
  if (!is.null(r)) {
    coloc_summary <- as.data.frame(r$summary)
    coloc_summary$window_kb <- hw / 1000
    write.csv(
      coloc_summary,
      file.path(out_dir, sprintf("coloc_win_%dkb_summary.csv", hw / 1000)),
      row.names = FALSE
    )
  }
}
