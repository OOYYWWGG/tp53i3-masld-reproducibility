args <- commandArgs(trailingOnly = TRUE)
raw_path <- args[1]
corrected_path <- args[2]
out_path <- args[3]

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

read_matrix <- function(path) {
  object <- readRDS(path)
  counts <- SeuratObject::LayerData(
    object,
    assay = "RNA",
    layer = "counts"
  )
  list(object = object, counts = counts)
}

aggregate_donor <- function(object, counts, matrix_type) {
  metadata <- object@meta.data
  tp53i3 <- as.numeric(counts["TP53I3", ])
  library_size <- as.numeric(Matrix::colSums(counts))
  donor <- data.frame(
    sample = as.character(metadata$sample),
    group = as.character(metadata$group),
    tp53i3_count = tp53i3,
    library_count = library_size,
    raw_median_nCount = as.numeric(metadata$nCount_RNA),
    stringsAsFactors = FALSE
  )
  rows <- lapply(split(donor, donor$sample), function(x) {
    if (length(unique(x$group)) != 1L) {
      stop("A donor was assigned to more than one group")
    }
    data.frame(
      matrix_type = matrix_type,
      sample = x$sample[1],
      group = unique(x$group),
      n_cells = nrow(x),
      sum_tp53i3 = sum(x$tp53i3_count),
      sum_library = sum(x$library_count),
      cpm10k = 10000 * sum(x$tp53i3_count) / sum(x$library_count),
      median_nCount = median(x$raw_median_nCount, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$log2cpm <- log2(result$cpm10k + 0.01)
  result
}

raw <- read_matrix(raw_path)
corrected <- read_matrix(corrected_path)
if (!identical(colnames(raw$counts), colnames(corrected$counts))) {
  stop("Raw and corrected matrices contain different nuclei")
}
if (!identical(
  as.character(raw$object@meta.data$sample),
  as.character(corrected$object@meta.data$sample)
)) {
  stop("Raw and corrected objects contain different donor assignments")
}

result <- rbind(
  aggregate_donor(raw$object, raw$counts, "matched_subset_raw"),
  aggregate_donor(
    corrected$object,
    corrected$counts,
    "matched_subset_decontx"
  )
)
result$group <- factor(
  result$group,
  levels = c("ctrl", "MASL", "eMASH", "aMASH")
)
result <- result[order(result$matrix_type, result$group, result$sample), ]

write.csv(
  result,
  out_path,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
