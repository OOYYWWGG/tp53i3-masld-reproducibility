## ============================================================================
## GSE281797: TP53I3 ~ fibrosis/NAS + hepatocyte + ductular cells
## This script reads the prepared composition output for GSE281797.
## Output: GSE281797_composition_models.csv
## ============================================================================
suppressPackageStartupMessages(library(data.table))

BASE <- "INPUT_DIR/composition_inputs"
BULK <- file.path(BASE, "gse281797")

## ---------- Clinical data ----------
lines <- readLines(gzfile(file.path(BULK, "GSE281797_series_matrix.txt.gz")))
gv <- function(prefix) {
  z <- grep(paste0("^", prefix), lines, value = TRUE)
  if (!length(z)) return(NULL)
  gsub('^"|"$', "", strsplit(z[1], "\t", fixed = TRUE)[[1]][-1])
}
gc_ <- function(pattern) {
  z <- grep(pattern, lines, value = TRUE)
  if (!length(z)) return(rep(NA_character_, length(gv("!Sample_title"))))
  v <- gsub('^"|"$', "", strsplit(z[1], "\t", fixed = TRUE)[[1]][-1])
  sub("^[^:]+: ", "", v)
}
meta <- data.frame(
  sample_id = sub("^HumanLiver_", "", gv("!Sample_title")),
  sex = gc_("gender:"), age = as.numeric(gc_("age at_biopsy:")),
  bmi = as.numeric(gc_("body mass_index:")), diabetes = as.integer(gc_("diabetes:")),
  steatosis = as.integer(gc_("steatosis grade:")), fibrosis = as.integer(gc_("fibrosis grade:")),
  inflammation = as.integer(gc_("inflammation grade:")),
  ballooning = as.integer(gc_("ballooning grade:")),
  nas = as.integer(gc_("nafld activity_score:")),
  diagnosis = gc_("diagnosis:"), stringsAsFactors = FALSE
)
meta$sex <- factor(meta$sex, levels = c("Female", "Male"))
meta$diagnosis <- factor(meta$diagnosis, levels = c("No Pathology", "MASL", "MASH"))

## ---------- Expression (log2 CPM) ----------
counts <- fread(file.path(BULK, "GSE281797_GeneCount.txt.gz"), data.table = FALSE)
ids <- sub("\\..*$", "", as.character(counts[[1]]))
m <- as.matrix(counts[, -1, drop = FALSE]); rownames(m) <- ids
tpm <- fread(file.path(BULK, "GSE281797_tpm.txt.gz"), data.table = FALSE)[, 1:2]
colnames(tpm) <- c("ensembl", "symbol"); tpm$ensembl <- sub("\\..*$", "", tpm$ensembl)
map <- tpm[!duplicated(tpm$ensembl), ]
sym <- map$symbol[match(rownames(m), map$ensembl)]
ok <- !is.na(sym) & sym != ""
m <- rowsum(m[ok, , drop = FALSE], group = sym[ok], reorder = FALSE)
cpm <- sweep(m, 2, colSums(m), "/") * 1e6
logcpm <- log2(cpm + 1)

## ---------- Composition ----------
prop <- read.csv(file.path(BULK, "GSE281797_music_proportions.csv"),
                 check.names = FALSE, row.names = 1)
hep <- prop[, grep("Hepatocyte", colnames(prop), ignore.case = TRUE)[1]]

duct_genes <- c("KRT7", "KRT19", "EPCAM", "SOX9")
dmat <- logcpm[intersect(duct_genes, rownames(logcpm)), , drop = FALSE]
dz <- t(scale(t(dmat)))
duct <- colMeans(dz, na.rm = TRUE)

meta <- meta[match(colnames(logcpm), meta$sample_id), ]
d <- data.frame(sample = colnames(logcpm), tp53 = logcpm["TP53I3", ],
                hep = hep[match(colnames(logcpm), rownames(prop))],
                duct = duct[colnames(logcpm)], meta)
d <- d[d$diagnosis %in% c("MASL", "MASH"), ]
d <- d[complete.cases(d[, c("tp53", "fibrosis", "nas", "hep", "duct")]), ]

z <- function(x) as.numeric(scale(x))
d$y <- z(d$tp53); d$fib <- z(d$fibrosis); d$nasz <- z(d$nas)
d$heps <- z(d$hep); d$ducts <- z(d$duct)

fit_show <- function(form, label) {
  m <- lm(form, data = d); co <- summary(m)$coefficients; ci <- confint(m)
  rows <- lapply(rownames(co)[-1], function(t) data.frame(
    model = label, term = t, beta = co[t, 1], lcl = ci[t, 1], ucl = ci[t, 2],
    p = co[t, 4], r2 = summary(m)$r.squared, n = nobs(m)))
  cat("\n==", label, "==\n")
  for (t in rownames(co)[-1])
    cat(sprintf("%-10s beta=%+.3f (%.3f, %.3f) P=%.3g\n", t, co[t,1], ci[t,1], ci[t,2], co[t,4]))
  cat("R2 =", round(summary(m)$r.squared, 3), " n =", nobs(m), "\n")
  do.call(rbind, rows)
}

res <- rbind(
  fit_show(y ~ fib, "fibrosis only"),
  fit_show(y ~ fib + heps, "fibrosis + hepatocyte"),
  fit_show(y ~ fib + heps + ducts, "fibrosis + hepatocyte + ductular"),
  fit_show(y ~ nasz, "NAS only"),
  fit_show(y ~ nasz + heps + ducts, "NAS + hepatocyte + ductular")
)
write.csv(res, file.path(BULK, "GSE281797_composition_models.csv"), row.names = FALSE)
cat("\nwritten GSE281797_composition_models.csv\n")
