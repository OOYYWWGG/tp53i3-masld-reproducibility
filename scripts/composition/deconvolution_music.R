suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(MuSiC)
  library(data.table)
  library(ggplot2)
})
options(encoding = "UTF-8")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: deconvolution_music.R <deconvolution-input-directory>")
}
base_dir <- args[1]
sc_rda <- file.path(base_dir, "sc_obj", "sc_anno.rda")
out_dir <- file.path(base_dir, "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(12345)
max_cells_per_type <- as.integer(Sys.getenv("DECONV_MAX_CELLS", "1000"))

cat("Loading single-cell reference:", sc_rda, "\n")
load(sc_rda)
if (!exists("sc_anno")) stop("sc_anno object not found")

meta_all <- sc_anno@meta.data
if (!"cell_type" %in% colnames(meta_all)) stop("cell_type missing from metadata")
if (!"sample" %in% colnames(meta_all)) stop("sample missing from metadata")

cell_types_all <- sort(unique(as.character(meta_all$cell_type)))
reference_counts <- as.data.frame(table(meta_all$cell_type))
colnames(reference_counts) <- c("cell_type", "n_cells_reference")
write.csv(reference_counts, file.path(out_dir, "reference_celltype_counts.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cell_list <- split(colnames(sc_anno), meta_all$cell_type)
cells_keep <- unlist(lapply(cell_list, function(x) {
  sample(x, min(length(x), max_cells_per_type))
}), use.names = FALSE)

counts_all <- LayerData(sc_anno, assay = "RNA", layer = "counts")
counts_sub <- counts_all[, cells_keep, drop = FALSE]
meta_sub <- meta_all[cells_keep, c("cell_type", "sample", "group"), drop = FALSE]
rm(sc_anno, counts_all, meta_all, cell_list)
gc()

sce <- SingleCellExperiment(
  assays = list(counts = counts_sub),
  colData = meta_sub
)
rm(counts_sub); gc()

cat("Reference cells used:", ncol(sce), "\n")
cat("Reference cell types:", paste(cell_types_all, collapse = ", "), "\n")

read_gencode_map <- function(path) {
  gtf <- data.table::fread(
    path, sep = "\t", header = FALSE, select = c(3, 9),
    colClasses = "character", showProgress = FALSE
  )
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
  mat <- rowsum(mat, group = symbols, reorder = FALSE)
  mat
}

read_gse135251 <- function(path, gtf_path) {
  counts <- data.table::fread(path, data.table = FALSE)
  gene_ids <- sub("\\..*$", "", as.character(counts[[1]]))
  mat <- as.matrix(counts[, -1, drop = FALSE])
  rownames(mat) <- gene_ids
  gtf_map <- read_gencode_map(gtf_path)
  symbols <- unname(gtf_map[rownames(mat)])
  mat <- collapse_rows(mat, symbols)
  cpm <- sweep(mat, 2, colSums(mat), "/") * 1e6
  list(linear = cpm, log2 = log2(cpm + 1))
}

read_microarray <- function(path) {
  expr <- data.table::fread(path, data.table = FALSE)
  symbols <- as.character(expr[[1]])
  mat <- as.matrix(expr[, -1, drop = FALSE])
  rownames(mat) <- symbols
  mat <- collapse_rows(mat, symbols)
  linear <- 2^mat - 1
  linear[linear < 0] <- 0
  linear <- sweep(linear, 2, colSums(linear), "/") * 1e6
  list(linear = linear, log2 = log2(linear + 1))
}

cat("Loading bulk cohorts...\n")
bulk <- list(
  GSE135251 = read_gse135251(
    file.path(base_dir, "bulk", "GSE135251_counts.csv.gz"),
    file.path(base_dir, "gencode", "gencode.v47.annotation.gtf.gz")
  ),
  GSE48452 = read_microarray(
    file.path(base_dir, "chip", "GSE48452_expression.csv.gz")
  ),
  GSE89632 = read_microarray(
    file.path(base_dir, "chip", "GSE89632_expression.csv.gz")
  ),
  GSE49541 = read_microarray(
    file.path(base_dir, "gse49541", "GSE49541_protein_coding.csv")
  )
)

clinical <- read.csv(
  file.path(base_dir, "clinical", "clinical_metadata_master.csv"),
  fileEncoding = "UTF-8-BOM", check.names = FALSE, stringsAsFactors = FALSE
)

all_props <- list()

for (dataset in names(bulk)) {
  cat("\n===== MuSiC:", dataset, "=====\n")
  bulk_linear <- bulk[[dataset]]$linear
  bulk_log2 <- bulk[[dataset]]$log2
  common <- intersect(rownames(bulk_linear), rownames(sce))
  cat("Gene overlap:", length(common), "\n")
  if (length(common) < 1000) stop(paste("Too few overlapping genes:", dataset))

  bulk_mat <- bulk_linear[common, , drop = FALSE]
  fit <- music_prop(
    bulk.mtx = bulk_mat,
    sc.sce = sce,
    clusters = "cell_type",
    samples = "sample",
    select.ct = cell_types_all,
    verbose = TRUE
  )

  props <- as.data.frame(fit[["Est.prop.weighted"]])
  props$sample <- rownames(props)
  for (ct in cell_types_all) {
    if (!ct %in% colnames(props)) props[[ct]] <- 0
  }
  props <- props[, c("sample", cell_types_all), drop = FALSE]

  clinical_sub <- clinical[clinical$dataset == dataset, c("gsm", "group_unified"), drop = FALSE]
  props <- merge(
    data.frame(dataset = dataset, gc = props$sample, props[, -1, drop = FALSE],
               check.names = FALSE),
    clinical_sub,
    by.x = "gc", by.y = "gsm", all.x = TRUE
  )
  colnames(props)[colnames(props) == "gc"] <- "sample"
  all_props[[dataset]] <- props
}

full_props <- do.call(rbind, all_props)
write.csv(full_props, file.path(out_dir, "deconvolution_full_proportions.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

summary_spec <- list(
  Hepatocyte = c("Hepatocyte"),
  Cholangiocyte = c("Cholangiocyte"),
  Monocyte = c("Monocyte"),
  Macrophage = c("Macrophage"),
  `T/NK` = c("T_cell", "NK_cell"),
  `B/Plasma` = c("B_cell", "Plasma"),
  Endothelial = c("Endothelial"),
  `HSC/Stellate` = c("HSC_Stellate")
)

summary_props <- full_props[, c("dataset", "sample", "group_unified"), drop = FALSE]
for (target in names(summary_spec)) {
  components <- summary_spec[[target]]
  present <- intersect(components, colnames(full_props))
  summary_props[[target]] <- if (length(present)) {
    rowSums(full_props[, present, drop = FALSE])
  } else {
    0
  }
}
other_types <- setdiff(
  setdiff(cell_types_all, unlist(summary_spec)),
  character(0)
)
summary_props$Other <- rowSums(full_props[, other_types, drop = FALSE])
summary_props <- summary_props[, c("dataset", "sample", "group_unified",
                                   names(summary_spec), "Other")]
write.csv(summary_props, file.path(out_dir, "deconvolution_summary_proportions.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Cohort-level reasonableness summary
summary_long <- reshape(
  summary_props,
  varying = c(names(summary_spec), "Other"),
  v.names = "proportion",
  timevar = "cell_type",
  times = c(names(summary_spec), "Other"),
  direction = "long"
)
cohort_summary <- aggregate(
  proportion ~ dataset + cell_type, data = summary_long,
  FUN = function(x) c(
    n = length(x),
    mean = mean(x),
    median = median(x),
    min = min(x),
    max = max(x),
    n_negative = sum(x < -1e-6),
    n_gt_0.8 = sum(x > 0.8)
  )
)
cohort_summary <- data.frame(
  cohort_summary[, c("dataset", "cell_type")],
  as.data.frame(cohort_summary$proportion)
)
write.csv(cohort_summary, file.path(out_dir, "deconvolution_cohort_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Boxplot of summarized proportions
plot_data <- summary_long
plot_data$dataset <- factor(plot_data$dataset, levels = names(bulk))
plot_data$cell_type <- factor(plot_data$cell_type, levels = c(names(summary_spec), "Other"))
p <- ggplot(plot_data, aes(x = cell_type, y = proportion, fill = dataset)) +
  geom_boxplot(outlier.size = 0.4, width = 0.72) +
  facet_wrap(~dataset, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = "Estimated proportion", fill = "Cohort")
ggsave(file.path(out_dir, "deconv_proportion_boxplot.pdf"), p, width = 10, height = 7)
ggsave(file.path(out_dir, "deconv_proportion_boxplot.png"), p, width = 10, height = 7, dpi = 300)

neg_count <- sum(as.matrix(full_props[, cell_types_all, drop = FALSE]) < -1e-6)
row_sum_range <- range(rowSums(full_props[, cell_types_all, drop = FALSE]))

lines <- c(
  "Bulk 去卷积（MuSiC）QC 报告",
  "============================================================",
  paste("生成日期：", Sys.Date()),
  paste("MuSiC 版本：", as.character(packageVersion("MuSiC"))),
  paste("单细胞参考：GSE289173，使用细胞数=", ncol(sce)),
  paste("每类最多下采样：", max_cells_per_type),
  paste("参考细胞类型：", paste(cell_types_all, collapse = ", ")),
  "",
  "一、方法",
  "1. 以 GSE289173 sc_anno.rda 的 counts 为参考，按 cell_type 下采样。",
  "2. MuSiC music_prop，clusters=cell_type，samples=sample，估计每个 bulk 样本的细胞比例。",
  "3. GSE135251 使用 counts -> CPM（线性尺度）；三个芯片队列使用 log2 -> 2^x-1 -> 按样本总和缩放为 1e6。",
  "4. 比例报告 8 个目标类型 + Other（其余参考细胞类型）。",
  "",
  "二、合理性检查",
  paste("负值总数：", neg_count),
  paste("每样本比例总和范围：", paste(round(row_sum_range, 6), collapse = " - ")),
  paste("基因重叠数：", paste(names(bulk), sapply(names(bulk), function(d) length(intersect(rownames(bulk[[d]]$linear), rownames(sce)))), sep = "=", collapse = "; ")),
  "",
  "三、队列比例汇总（mean / median / min / max）",
  paste(
    cohort_summary$dataset, cohort_summary$cell_type,
    paste0(round(cohort_summary$mean, 4), " / ", round(cohort_summary$median, 4),
           " / ", round(cohort_summary$min, 4), " / ", round(cohort_summary$max, 4)),
    sep = "\t"
  ),
  "",
  "四、限制",
  "1. MuSiC 假设参考与 bulk 的细胞类型和表达程序可比；snRNA 参考与 bulk 组织仍有差异。",
  "2. 三个芯片队列由 log2 转线性后使用，和 RNA-seq 的绝对尺度不完全一致；比例应作相对估计。",
  "3. 肝 bulk 中 Hepatocyte 比例高是合理的；若某类型大量为零或极端，需要结合独立参考判断。"
)
writeLines(enc2utf8(lines), file.path(out_dir, "deconvolution_QC.txt"), useBytes = TRUE)

cat("DONE. Outputs:", out_dir, "\n")
sessionInfo()
