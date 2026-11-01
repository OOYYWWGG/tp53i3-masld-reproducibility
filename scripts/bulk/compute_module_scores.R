suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else Sys.getenv("TEMP")
out_dir <- file.path(base_dir, "module_score_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_gencode_map <- function(path) {
  gtf <- data.table::fread(path, sep = "\t", header = FALSE,
                           select = c(3, 9), colClasses = "character",
                           showProgress = FALSE)
  gene_lines <- gtf[V3 == "gene"]
  gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", gene_lines$V9)
  gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", gene_lines$V9)
  map <- setNames(gene_name, sub("\\..*$", "", gene_id))
  map[!is.na(names(map)) & names(map) != ""]
}

collapse_rows <- function(mat, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  rowsum(mat, group = symbols, reorder = FALSE)
}

read_gse <- function(path, gtf) {
  counts <- data.table::fread(path, data.table = FALSE)
  ids <- sub("\\..*$", "", as.character(counts[[1]]))
  mat <- as.matrix(counts[, -1, drop = FALSE])
  rownames(mat) <- ids
  map <- read_gencode_map(gtf)
  symbols <- unname(map[rownames(mat)])
  mat <- collapse_rows(mat, symbols)
  cpm <- sweep(mat, 2, colSums(mat), "/") * 1e6
  log2(cpm + 1)
}

read_micro <- function(path) {
  expr <- data.table::fread(path, data.table = FALSE)
  symbols <- as.character(expr[[1]])
  mat <- as.matrix(expr[, -1, drop = FALSE])
  rownames(mat) <- symbols
  collapse_rows(mat, symbols)
}

clin <- read.csv(file.path(base_dir, "clinical.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
comp <- read.csv(file.path(base_dir, "composition.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
duct <- read.csv(file.path(base_dir, "ductular.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
names(clin)[1] <- "dataset"
names(comp)[1] <- "dataset"
names(duct)[1] <- "dataset"

d <- merge(clin,
           comp[, c("dataset", "gsm", "hepatocyte_fraction",
                    "cholangiocyte_fraction", "immune_fraction")],
           by = c("dataset", "gsm"), all.x = TRUE)
d <- merge(d, duct[, c("dataset", "gsm", "ductular_score")],
           by = c("dataset", "gsm"), all.x = TRUE)
tp53 <- read.csv(file.path(base_dir, "tp53_analysis.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
names(tp53)[1] <- "dataset"
tp53 <- tp53[, c("dataset", "gsm", "tp53i3")]
d <- merge(d, tp53, by = c("dataset", "gsm"), all.x = TRUE)
d <- d[d$dataset %in% c("GSE135251", "GSE48452", "GSE89632", "GSE49541"), ]
d <- d[!(d$group_unified %in% c("control", "obese")), ]
d$sex <- factor(ifelse(d$sex %in% c("male", "M", "1"), "male",
                       ifelse(d$sex %in% c("female", "F", "0"), "female", NA)))
d$age <- as.numeric(d$age_years)
d$bmi <- as.numeric(d$bmi_kg_m2)
d$nas <- as.numeric(d$nas)
d$fib <- as.numeric(d$fibrosis_stage_numeric)
d$fib_bin <- ifelse(is.na(d$fib), NA, ifelse(d$fib >= 2, 1, 0))
idx49541 <- d$dataset == "GSE49541"
d$fib_bin[idx49541] <- ifelse(d$group_unified[idx49541] == "NAFLD_advanced", 1, 0)
d$fib_ord <- d$fib
d$fib_ord[!is.na(d$fib_ord) & d$fib_ord %in% c(0.5, 3.5)] <- NA
d$nas_bin <- ifelse(is.na(d$nas), NA, ifelse(d$nas >= 5, 1, 0))

p53_genes <- c(
  "CDKN1A", "MDM2", "BBC3", "GADD45A", "RRM2B", "SESN1", "ZMAT3",
  "TIGAR", "FDXR", "PHLDA3", "TP53INP1", "XPC", "DDB2"
)
zonation_genes <- c(
  "CYP3A4", "CYP2E1", "GLUL", "SLCO1B3", "IGFBP1", "CYP2A6", "CYP2B6"
)
exploratory_genes <- c(
  "ALB", "APOA1", "TF", "TTR", "CYP2E1", "APOA2", "APOB", "SERPINA1",
  "CYP2A6", "TP53INP1", "ZMAT3", "DTNA", "IL32", "SPATA21", "BICC1",
  "SQSTM1", "ARHGAP27", "WDPCP", "APOC1", "AKR1C2"
)

logs <- list(
  GSE135251 = read_gse(
    file.path(base_dir, "gse135251.csv.gz"),
    file.path(base_dir, "gencode.v47.annotation.gtf.gz")
  ),
  GSE48452 = read_micro(file.path(base_dir, "gse48452.csv.gz")),
  GSE89632 = read_micro(file.path(base_dir, "gse89632.csv.gz")),
  GSE49541 = read_micro(file.path(base_dir, "gse49541_protein_coding.csv"))
)

score_rows <- list()
qc_rows <- list()
for (ds in names(logs)) {
  mat <- logs[[ds]]
  samples <- d$gsm[d$dataset == ds]
  missing_samples <- setdiff(samples, colnames(mat))
  if (length(missing_samples)) {
    stop(paste("Missing samples in", ds, ":", paste(missing_samples, collapse = ",")))
  }
  mat <- mat[, samples, drop = FALSE]

  gene_sets <- list(
    p53_module = intersect(p53_genes, rownames(mat)),
    zonation_module = intersect(zonation_genes, rownames(mat)),
    exploratory_signature = intersect(exploratory_genes, rownames(mat))
  )
  if (any(vapply(gene_sets, length, integer(1)) < 5)) {
    stop(paste("Too few genes in", ds))
  }

  param <- ssgseaParam(mat, gene_sets, normalize = TRUE, alpha = 0.25,
                       verbose = FALSE)
  ss <- gsva(param)

  mean_z <- function(genes) {
    g <- intersect(genes, rownames(mat))
    z <- t(scale(t(mat[g, , drop = FALSE])))
    z[is.na(z)] <- 0
    colMeans(z)
  }
  p53_mz <- mean_z(p53_genes)
  zon_mz <- mean_z(zonation_genes)
  expl_mz <- mean_z(exploratory_genes)

  score_rows[[ds]] <- data.frame(
    dataset = ds,
    gsm = colnames(mat),
    p53_ssgsea = as.numeric(ss["p53_module", ]),
    zonation_ssgsea = as.numeric(ss["zonation_module", ]),
    exploratory_ssgsea = as.numeric(ss["exploratory_signature", ]),
    p53_mean_z = as.numeric(p53_mz),
    zonation_mean_z = as.numeric(zon_mz),
    exploratory_mean_z = as.numeric(expl_mz),
    stringsAsFactors = FALSE
  )

  qc_rows[[ds]] <- data.frame(
    dataset = ds,
    n_samples = ncol(mat),
    n_genes_matrix = nrow(mat),
    p53_genes = length(gene_sets$p53_module),
    zonation_genes = length(gene_sets$zonation_module),
    exploratory_genes = length(gene_sets$exploratory_signature),
    p53_ssgsea_sd = sd(ss["p53_module", ]),
    zonation_ssgsea_sd = sd(ss["zonation_module", ]),
    cor_p53_zonation_ssgsea = cor(ss["p53_module", ],
                                  ss["zonation_module", ],
                                  method = "spearman"),
    cor_ssgsea_meanz_p53 = cor(ss["p53_module", ], p53_mz,
                               method = "spearman"),
    cor_ssgsea_meanz_zonation = cor(ss["zonation_module", ], zon_mz,
                                    method = "spearman"),
    stringsAsFactors = FALSE
  )
}

scores <- do.call(rbind, score_rows)
d <- merge(d, scores, by = c("dataset", "gsm"), all.x = TRUE)

for (ds in unique(d$dataset)) {
  i <- d$dataset == ds
  d$p53_ssgsea_z[i] <- as.numeric(scale(d$p53_ssgsea[i]))
  d$zonation_ssgsea_z[i] <- as.numeric(scale(d$zonation_ssgsea[i]))
  d$state_ssgsea_z[i] <- as.numeric(scale(
    d$p53_ssgsea_z[i] - d$zonation_ssgsea_z[i]
  ))
  d$state_mean_z[i] <- as.numeric(scale(
    d$p53_mean_z[i] - d$zonation_mean_z[i]
  ))
  d$exploratory_ssgsea_z[i] <- as.numeric(scale(d$exploratory_ssgsea[i]))
}

write.csv(d, file.path(out_dir, "bulk_module_scores.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

qc <- do.call(rbind, qc_rows)
write.csv(qc, file.path(out_dir, "bulk_module_score_QC.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("Scoring finished\n")
