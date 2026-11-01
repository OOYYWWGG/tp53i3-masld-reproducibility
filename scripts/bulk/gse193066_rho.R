## Place the public GSE193066 input files under INPUT_DIR/gse193066.
base <- "INPUT_DIR/gse193066"
sm <- readLines(gzfile(file.path(base, "GSE193066_series_matrix.txt.gz")))
gr <- function(prefix) {
  z <- grep(paste0("^", prefix), sm, value = TRUE)
  gsub('^"|"$', "", strsplit(z[1], "\t", fixed = TRUE)[[1]][-1])
}
sf <- function(x, p) sub(p, "", x)
meta <- data.frame(
  sample = gr("!Sample_title"),
  fib = as.integer(sf(gr("!Sample_characteristics_ch1.*fibrosis stage:"), "^.*: ")),
  nas = as.integer(sf(gr("!Sample_characteristics_ch1.*nafld activity score:"), "^.*: ")),
  age = as.integer(sf(gr("!Sample_characteristics_ch1.*age:"), "^.*: ")),
  sex = sf(gr("!Sample_characteristics_ch1.*Sex:"), "^.*: "),
  stringsAsFactors = FALSE
)
g <- read.delim(
  gzfile(file.path(base, "GSE193066_NAFLD.HUn106.gct.gz")),
  skip = 2, header = TRUE, check.names = FALSE
)
tp <- g[g$gene_name == "TP53I3", , drop = FALSE]
meta <- meta[meta$sample %in% colnames(g), ]
meta$tp53 <- as.numeric(tp[1, meta$sample])
ct <- cor.test(meta$tp53, meta$nas, method = "spearman", exact = FALSE)
write.csv(meta, "INPUT_DIR/gse193066_meta.csv", row.names = FALSE)
