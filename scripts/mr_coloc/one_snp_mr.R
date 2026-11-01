suppressPackageStartupMessages(library(TwoSampleMR))

out_dir <- "results/mr_coloc"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- data.frame(
  SNP = "rs11902756", beta.exposure = -0.255414, se.exposure = 0.0522525,
  pval.exposure = 2.36e-6, effect_allele.exposure = "G", other_allele.exposure = "A",
  eaf.exposure = 0.339, exposure = "x", id.exposure = "x",
  beta.outcome = 0.10, se.outcome = 0.05, pval.outcome = 0.05,
  effect_allele.outcome = "G", other_allele.outcome = "A", eaf.outcome = 0.34,
  outcome = "y", id.outcome = "y", mr_keep = TRUE, stringsAsFactors = FALSE
)

run_one <- function(tag, method_list) {
  result <- tryCatch(
    mr(dat, method_list = method_list),
    error = function(e) data.frame(error = conditionMessage(e))
  )
  write.csv(
    result,
    file.path(out_dir, paste0("one_snp_mr_", tag, ".csv")),
    row.names = FALSE
  )
  invisible(result)
}

run_one("ivw", "mr_ivw")
run_one("wald_ratio", "mr_wald_ratio")
run_one("ivw_wald_ratio", c("mr_ivw", "mr_wald_ratio"))
