#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(lme4)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
preflight <- "--preflight" %in% args

input_path <- Sys.getenv("SC_INPUT", "/data/src/sc/sc_anno_score.rda")
out_dir <- Sys.getenv("OUT_DIR", "/data/out_singlecell")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf(...), "\n")

read_seurat <- function(path) {
  if (!file.exists(path)) stop("输入不存在: ", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (!inherits(obj, "Seurat")) stop("RDS 不是 Seurat 对象: ", class(obj)[1])
    return(obj)
  }
  if (grepl("\\.rda$", path, ignore.case = TRUE)) {
    e <- new.env()
    load(path, envir = e)
    objs <- ls(e)
    seurat_names <- objs[vapply(objs, function(nm) inherits(get(nm, envir = e), "Seurat"), logical(1))]
    if (!length(seurat_names)) stop("RDA 中没有 Seurat 对象；对象名: ", paste(objs, collapse = ", "))
    sizes <- vapply(seurat_names, function(nm) ncol(get(nm, envir = e)), numeric(1))
    pick <- seurat_names[which.max(sizes)]
    log_msg("RDA 中发现 Seurat 对象: %s；选择 %s（%d 个细胞）",
            paste(seurat_names, collapse = ", "), pick, as.integer(max(sizes)))
    return(get(pick, envir = e))
  }
  stop("只支持 .rds 或 .rda: ", path)
}

pick_col <- function(meta, candidates, label) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (!length(hit)) {
    stop("找不到 ", label, " 列。候选: ", paste(candidates, collapse = ", "),
         "；实际列: ", paste(colnames(meta), collapse = ", "))
  }
  hit[1]
}

norm_group <- function(x) {
  x <- as.character(x)
  x[tolower(x) %in% c("control", "ctrl", "healthy", "normal")] <- "ctrl"
  x[tolower(x) %in% c("masl", "nafl", "steatosis")] <- "MASL"
  x[tolower(x) %in% c("emash", "early_mash", "early mash")] <- "eMASH"
  x[tolower(x) %in% c("amash", "advanced_mash", "advanced mash")] <- "aMASH"
  x
}

get_counts <- function(obj, assay_name) {
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    z <- tryCatch(
      SeuratObject::LayerData(obj, assay = assay_name, layer = "counts"),
      error = function(e) NULL
    )
    if (!is.null(z)) return(z)
  }
  GetAssayData(obj, assay = assay_name, slot = "counts")
}

get_ncount <- function(obj, meta, counts) {
  if ("nCount_RNA" %in% colnames(meta)) return(as.numeric(meta$nCount_RNA))
  Matrix::colSums(counts)
}

sample_level <- function(idx, gene, counts, meta, sample_col, group_col, ncount) {
  idx <- idx[!is.na(idx)]
  if (!length(idx) || !(gene %in% rownames(counts))) return(NULL)
  pos <- as.numeric(counts[gene, idx] > 0)
  smp <- as.character(meta[[sample_col]][idx])
  grp <- norm_group(meta[[group_col]][idx])
  nc <- as.numeric(ncount[idx])
  d <- data.frame(sample = smp, group = grp, pos = pos, nCount = nc,
                  stringsAsFactors = FALSE)
  rows <- lapply(split(d, d$sample), function(z) {
    data.frame(
      sample = z$sample[1],
      group = z$group[1],
      n_cells = nrow(z),
      n_positive = sum(z$pos),
      n_negative = sum(1 - z$pos),
      positive_pct = 100 * mean(z$pos),
      median_nCount = median(z$nCount, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$gene <- gene
  out
}

fit_binomial <- function(d, label) {
  d <- d[d$n_cells >= 5, , drop = FALSE]
  if (nrow(d) < 8 || length(unique(d$group)) < 2) {
    return(list(model = NULL, coef = data.frame(), note = "样本不足"))
  }
  d$group <- factor(d$group, levels = c("ctrl", setdiff(unique(d$group), "ctrl")))
  d$depth_scaled <- d$median_nCount / median(d$median_nCount, na.rm = TRUE)
  fit <- try(
    glmer(
      cbind(n_positive, n_negative) ~ group + offset(log(depth_scaled)) + (1 | sample),
      data = d, family = binomial, nAGQ = 0,
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    fit <- glm(
      cbind(n_positive, n_negative) ~ group + offset(log(depth_scaled)),
      data = d, family = quasibinomial()
    )
  }
  co <- summary(fit)$coefficients
  terms <- grep("^group", rownames(co), value = TRUE)
  if (!length(terms)) return(list(model = fit, coef = data.frame(), note = "无 group 项"))
  out <- do.call(rbind, lapply(terms, function(term) {
    est <- co[term, "Estimate"]
    se <- co[term, "Std. Error"]
    data.frame(
      term = term,
      OR = exp(est),
      lcl = exp(est - qnorm(0.975) * se),
      ucl = exp(est + qnorm(0.975) * se),
      p = co[term, ncol(co)],
      n_samples = nrow(d),
      n_cells = sum(d$n_cells),
      stringsAsFactors = FALSE
    )
  }))
  list(model = fit, coef = out, note = label)
}

run_hurdle <- function(idx, gene, counts, meta, sample_col, group_col, ncount) {
  base <- sample_level(idx, gene, counts, meta, sample_col, group_col, ncount)
  if (is.null(base)) return(list(detection = data.frame(), expression = data.frame()))
  det <- fit_binomial(base, "detection")$coef
  if (nrow(det)) det$component <- "detection"

  pos_idx <- idx[counts[gene, idx] > 0]
  if (!length(pos_idx)) return(list(detection = det, expression = data.frame()))
  smp <- as.character(meta[[sample_col]][pos_idx])
  grp <- norm_group(meta[[group_col]][pos_idx])
  expr <- log1p(as.numeric(counts[gene, pos_idx]))
  ed <- data.frame(sample = smp, group = grp, expr = expr, stringsAsFactors = FALSE)
  erows <- lapply(split(ed, ed$sample), function(z) data.frame(
    sample = z$sample[1], group = z$group[1], mean_log1p = mean(z$expr),
    n_positive = nrow(z), stringsAsFactors = FALSE
  ))
  e <- do.call(rbind, erows)
  if (nrow(e) < 8 || length(unique(e$group)) < 2) {
    return(list(detection = det, expression = data.frame()))
  }
  e$group <- factor(e$group, levels = c("ctrl", setdiff(unique(e$group), "ctrl")))
  fit <- lm(mean_log1p ~ group, data = e)
  co <- summary(fit)$coefficients
  terms <- grep("^group", rownames(co), value = TRUE)
  eout <- do.call(rbind, lapply(terms, function(term) data.frame(
    term = term,
    beta = co[term, "Estimate"],
    lcl = co[term, "Estimate"] - qnorm(0.975) * co[term, "Std. Error"],
    ucl = co[term, "Estimate"] + qnorm(0.975) * co[term, "Std. Error"],
    p = co[term, ncol(co)],
    n_samples = nrow(e),
    stringsAsFactors = FALSE
  )))
  if (nrow(eout)) eout$component <- "expression_positive"
  list(detection = det, expression = eout)
}

obj <- read_seurat(input_path)
meta <- obj@meta.data
celltype_col <- pick_col(meta, c("cell_type", "celltype", "cell_type_final", "annotation", "celltype_final"), "cell type")
sample_col <- pick_col(meta, c("sample", "sample_id", "donor", "orig.ident"), "sample")
group_col <- pick_col(meta, c("group", "group_unified", "stage", "disease_group"), "group")

assay_name <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
counts <- get_counts(obj, assay_name)
ncount <- get_ncount(obj, meta, counts)
celltype <- as.character(meta[[celltype_col]])
group <- norm_group(meta[[group_col]])

celltype_group <- as.data.frame(table(celltype = celltype, group = group), stringsAsFactors = FALSE)
write.csv(celltype_group, file.path(out_dir, "celltype_group_counts.csv"), row.names = FALSE, fileEncoding = "UTF-8")

log_msg("输入: %s", input_path)
log_msg("assay: %s; 细胞数: %d; 基因数: %d", assay_name, ncol(counts), nrow(counts))
log_msg("cell type 列: %s；sample 列: %s；group 列: %s", celltype_col, sample_col, group_col)

if (preflight) {
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "run_session_info.txt"))
  log_msg("PREFLIGHT_ONLY；结果目录: %s", out_dir)
  quit(save = "no", status = 0)
}

p53_module <- c("CDKN1A", "MDM2", "BBC3", "GADD45A", "RRM2B", "SESN1",
                "ZMAT3", "TIGAR", "FDXR", "PHLDA3", "TP53INP1", "XPC", "DDB2")

is_hep <- grepl("hep", celltype, ignore.case = TRUE)
is_chol <- grepl("chol", celltype, ignore.case = TRUE)
if (!any(is_hep)) stop("没有识别到肝细胞；请检查 cell type 列")
if (!any(is_chol)) warning("没有识别到胆管细胞；胆管 stage 分析将跳过")

## expression-matched negative panel, selected in hepatocytes
hep_idx <- which(is_hep)
if (!("TP53I3" %in% rownames(counts))) stop("counts 中没有 TP53I3")
set.seed(20261018)
hep_for_panel <- if (length(hep_idx) > 8000) sample(hep_idx, 8000) else hep_idx
log_msg("选择阴性 panel：从 %d 个肝细胞中抽取 %d 个用于表达匹配",
        length(hep_idx), length(hep_for_panel))
det <- as.numeric(rowMeans(counts[, hep_for_panel, drop = FALSE] > 0))
mean_log <- as.numeric(rowMeans(log1p(counts[, hep_for_panel, drop = FALSE])))
all_cand <- data.frame(gene = rownames(counts), detection = det, mean_log = mean_log,
                       stringsAsFactors = FALSE)
exclude <- unique(c("TP53I3", p53_module, grep("^MT-", all_cand$gene, value = TRUE),
                    grep("^RP[SL]", all_cand$gene, value = TRUE), "MKI67", "TOP2A"))
cand <- all_cand[!all_cand$gene %in% exclude & is.finite(all_cand$detection) &
               is.finite(all_cand$mean_log) & all_cand$detection > 0.005 &
               all_cand$detection < 0.8, ]
tp <- all_cand[all_cand$gene == "TP53I3", ]
if (!nrow(tp)) stop("TP53I3 没有通过表达过滤")
tp_det <- as.numeric(tp$detection[1])
tp_mean <- as.numeric(tp$mean_log[1])
mean_range <- max(cand$mean_log, na.rm = TRUE) - min(cand$mean_log, na.rm = TRUE)
if (!is.finite(mean_range) || mean_range == 0) mean_range <- 1
cand$dist <- abs(log10((cand$detection + 1e-4) / (tp_det + 1e-4))) +
  abs(cand$mean_log - tp_mean) / mean_range
neg <- head(cand[order(cand$dist), ], 20)
log_msg("阴性 panel 距离范围：%s", paste(signif(range(neg$dist), 3), collapse = "–"))
write.csv(neg, file.path(out_dir, "negative_panel_genes.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## sample-level tables
hep_tab <- sample_level(hep_idx, "TP53I3", counts, meta, sample_col, group_col, ncount)
write.csv(hep_tab, file.path(out_dir, "sample_level_hepatocyte.csv"), row.names = FALSE, fileEncoding = "UTF-8")
hep_fit <- fit_binomial(hep_tab, "hepatocyte")
write.csv(hep_fit$coef, file.path(out_dir, "hepatocyte_glmm.csv"), row.names = FALSE, fileEncoding = "UTF-8")

chol_tab <- if (any(is_chol)) sample_level(which(is_chol), "TP53I3", counts, meta, sample_col, group_col, ncount) else NULL
if (!is.null(chol_tab)) write.csv(chol_tab, file.path(out_dir, "sample_level_cholangiocyte.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## bile duct stage GLMM
if (!is.null(chol_tab)) {
  log_msg("拟合胆管细胞 stage 模型")
  bd <- fit_binomial(chol_tab, "cholangiocyte")
  write.csv(bd$coef, file.path(out_dir, "bile_duct_glmm.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

## negative panel GLMM
neg_rows <- list()
for (i in seq_along(neg$gene)) {
  g <- neg$gene[i]
  log_msg("negative panel %d/%d: %s", i, length(neg$gene), g)
  tab <- sample_level(hep_idx, g, counts, meta, sample_col, group_col, ncount)
  if (is.null(tab)) next
  fit <- fit_binomial(tab, "negative_panel")
  if (nrow(fit$coef)) {
    fit$coef$gene <- g
    neg_rows[[length(neg_rows) + 1]] <- fit$coef
  }
}
neg_tab <- if (length(neg_rows)) do.call(rbind, neg_rows) else data.frame()
write.csv(neg_tab, file.path(out_dir, "negative_panel_glmm.csv"), row.names = FALSE, fileEncoding = "UTF-8")

neg_summary <- data.frame()
if (nrow(neg_tab)) {
  hep_or <- hep_fit$coef$OR[hep_fit$coef$term == "groupaMASH"]
  neg_or <- neg_tab$OR[neg_tab$term == "groupaMASH"]
  neg_summary <- data.frame(
    tp53i3_aMASH_OR = if (length(hep_or)) hep_or[1] else NA_real_,
    negative_panel_n = length(unique(neg_tab$gene)),
    negative_panel_median_OR = median(neg_or, na.rm = TRUE),
    negative_panel_min_OR = min(neg_or, na.rm = TRUE),
    negative_panel_max_OR = max(neg_or, na.rm = TRUE),
    negative_panel_percentile_below_TP53I3 = mean(neg_or < hep_or[1], na.rm = TRUE),
    negative_panel_sig_aMASH = sum(neg_tab$term == "groupaMASH" & neg_tab$p < 0.05, na.rm = TRUE),
    negative_panel_sig_eMASH = sum(neg_tab$term == "groupeMASH" & neg_tab$p < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  write.csv(neg_summary, file.path(out_dir, "negative_panel_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

## hurdle for TP53I3 in hepatocytes and cholangiocytes
log_msg("拟合 TP53I3 hurdle 模型")
hep_h <- run_hurdle(hep_idx, "TP53I3", counts, meta, sample_col, group_col, ncount)
chol_h <- if (any(is_chol)) run_hurdle(which(is_chol), "TP53I3", counts, meta, sample_col, group_col, ncount) else list(detection = data.frame(), expression = data.frame())

write.csv(rbind(
  if (nrow(hep_h$detection)) transform(hep_h$detection, cell_type = "Hepatocyte"),
  if (nrow(chol_h$detection)) transform(chol_h$detection, cell_type = "Cholangiocyte")
), file.path(out_dir, "hurdle_detection.csv"), row.names = FALSE, fileEncoding = "UTF-8")

write.csv(rbind(
  if (nrow(hep_h$expression)) transform(hep_h$expression, cell_type = "Hepatocyte"),
  if (nrow(chol_h$expression)) transform(chol_h$expression, cell_type = "Cholangiocyte")
), file.path(out_dir, "hurdle_expression.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## plots
if (!is.null(chol_tab)) {
  p1 <- ggplot(chol_tab, aes(x = group, y = positive_pct, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 1.5, alpha = 0.7) +
    labs(x = NULL, y = "TP53I3+ cholangiocytes (%)", title = "Cholangiocyte TP53I3 by stage") +
    theme_classic(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(out_dir, "fig_bile_duct_rate.pdf"), p1, width = 5.5, height = 4.2)
  ggsave(file.path(out_dir, "fig_bile_duct_rate.png"), p1, width = 5.5, height = 4.2, dpi = 300)
}

if (nrow(neg_tab)) {
  p2 <- ggplot(neg_tab, aes(x = OR, y = reorder(gene, OR))) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.2) +
    labs(x = "OR per stage group", y = NULL, title = "Expression-matched negative panel") +
    theme_classic(base_size = 10)
  ggsave(file.path(out_dir, "fig_negative_panel_effects.pdf"), p2, width = 7, height = 5)
  ggsave(file.path(out_dir, "fig_negative_panel_effects.png"), p2, width = 7, height = 5, dpi = 300)
}

## save a small object for reuse
keep_cells <- colnames(obj)[is_hep | is_chol]
keep_features <- unique(c("TP53I3", p53_module, neg$gene))
keep_features <- keep_features[keep_features %in% rownames(obj)]
small <- subset(obj, cells = keep_cells, features = keep_features)
saveRDS(small, file.path(out_dir, "sc_small_hep_chol.rds"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "run_session_info.txt"))
writeLines(c(
  paste("input:", input_path),
  paste("assay:", assay_name),
  paste("celltype column:", celltype_col),
  paste("sample column:", sample_col),
  paste("group column:", group_col),
  paste("hepatocytes:", sum(is_hep)),
  paste("cholangiocytes:", sum(is_chol)),
  paste("negative panel genes:", paste(neg$gene, collapse = ","))
), file.path(out_dir, "run_summary.txt"))

log_msg("DONE；输出目录: %s", out_dir)
