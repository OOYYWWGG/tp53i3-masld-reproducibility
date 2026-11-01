suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
})

args <- commandArgs(trailingOnly = TRUE)
incremental_dir <- args[1]
gse281797_dir <- args[2]
final_bulk_scores <- args[3]
out_dir <- args[4]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clean_names <- function(data) {
  names(data)[1] <- sub("^\ufeff", "", names(data)[1])
  data
}

read_csv_clean <- function(path) {
  clean_names(read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  ))
}

p53_genes <- c(
  "CDKN1A", "MDM2", "BBC3", "GADD45A", "RRM2B", "SESN1", "ZMAT3",
  "TIGAR", "FDXR", "PHLDA3", "TP53INP1", "XPC", "DDB2"
)

scores <- read_csv_clean(
  final_bulk_scores
)
scores <- scores[, c("dataset", "gsm", "p53_ssgsea_z")]
tp53 <- read_csv_clean(file.path(incremental_dir, "TP53I3_analysis_dataset.csv"))
clinical <- read_csv_clean(file.path(
  incremental_dir,
  "clinical_metadata_master.csv"
))
deconvolution <- read_csv_clean(file.path(
  incremental_dir,
  "deconvolution_full_proportions.csv"
))
ductular <- read_csv_clean(file.path(
  incremental_dir,
  "ductular_score_bulk.csv"
))

data <- merge(
  scores,
  tp53[, c("dataset", "gsm", "tp53i3", "nas", "group_unified")],
  by = c("dataset", "gsm")
)
data <- merge(
  data,
  clinical[, c(
    "dataset", "gsm", "age_years", "sex", "bmi_kg_m2", "diabetes",
    "alt_u_l", "ast_u_l"
  )],
  by = c("dataset", "gsm"),
  all.x = TRUE
)
data <- merge(
  data,
  deconvolution[, c("dataset", "sample", "Hepatocyte")],
  by.x = c("dataset", "gsm"),
  by.y = c("dataset", "sample"),
  all.x = TRUE
)
names(data)[names(data) == "Hepatocyte"] <- "hepatocyte_fraction"
data <- merge(
  data,
  ductular[, c("dataset", "gsm", "ductular_score")],
  by = c("dataset", "gsm"),
  all.x = TRUE
)

data <- data[
  !is.na(data$group_unified) &
    data$group_unified %in% c("NAFL", "NASH") &
    data$dataset %in% c("GSE135251", "GSE48452", "GSE89632"),
]
data$diabetes <- factor(
  ifelse(data$diabetes == "yes", "yes", "no"),
  levels = c("no", "yes")
)
data$sex <- factor(data$sex, levels = c("female", "male"))
data$nas_bin <- ifelse(is.na(data$nas), NA, as.integer(data$nas >= 5))

for (dataset in unique(data$dataset)) {
  rows <- data$dataset == dataset
  data$tp53i3_z[rows] <- as.numeric(scale(data$tp53i3[rows]))
  data$ductular_z[rows] <- as.numeric(scale(data$ductular_score[rows]))
  data$p53_ssgsea_z[rows] <- as.numeric(scale(data$p53_ssgsea_z[rows]))
}

read_series <- function(path) {
  lines <- readLines(gzfile(path))
  sample_values <- function(prefix) {
    selected <- grep(paste0("^", prefix), lines, value = TRUE)
    if (!length(selected)) return(NULL)
    gsub('^"|"$', "", strsplit(selected[1], "\t", fixed = TRUE)[[1]][-1])
  }
  characteristic <- function(pattern) {
    selected <- grep(pattern, lines, value = TRUE)
    if (!length(selected)) {
      return(rep(NA_character_, length(sample_values("!Sample_title"))))
    }
    values <- gsub(
      '^"|"$',
      "",
      strsplit(selected[1], "\t", fixed = TRUE)[[1]][-1]
    )
    sub("^[^:]+: ", "", values)
  }
  data.frame(
    sample_id = sub("^HumanLiver_", "", sample_values("!Sample_title")),
    sex = characteristic("gender:"),
    age_years = as.numeric(characteristic("age at_biopsy:")),
    bmi_kg_m2 = as.numeric(characteristic("body mass_index:")),
    diabetes = as.integer(characteristic("diabetes:")),
    steatosis = as.integer(characteristic("steatosis grade:")),
    fibrosis = as.integer(characteristic("fibrosis grade:")),
    inflammation = as.integer(characteristic("inflammation grade:")),
    ballooning = as.integer(characteristic("ballooning grade:")),
    nas = as.integer(characteristic("nafld activity_score:")),
    diagnosis = characteristic("diagnosis:"),
    stringsAsFactors = FALSE
  )
}

meta_281797 <- read_series(file.path(
  gse281797_dir,
  "GSE281797_series_matrix.txt.gz"
))
selected_samples_281797 <- meta_281797$sample_id[
  meta_281797$diagnosis %in% c("MASL", "MASH") &
    !is.na(meta_281797$fibrosis) &
    meta_281797$fibrosis <= 3
]

counts_281797 <- fread(
  file.path(gse281797_dir, "GSE281797_GeneCount.txt.gz"),
  data.table = FALSE
)
ensembl <- sub("\\..*$", "", as.character(counts_281797[[1]]))
count_matrix <- as.matrix(counts_281797[, -1, drop = FALSE])
rownames(count_matrix) <- ensembl
if (!all(selected_samples_281797 %in% colnames(count_matrix))) {
  stop("GSE281797 disease samples do not match the expression matrix")
}

tpm_map <- fread(
  file.path(gse281797_dir, "GSE281797_tpm.txt.gz"),
  data.table = FALSE
)[, 1:2]
colnames(tpm_map) <- c("ensembl", "symbol")
tpm_map$ensembl <- sub("\\..*$", "", tpm_map$ensembl)
tpm_map <- tpm_map[!duplicated(tpm_map$ensembl), ]
symbol <- tpm_map$symbol[match(rownames(count_matrix), tpm_map$ensembl)]
keep <- !is.na(symbol) & symbol != ""
count_matrix <- rowsum(
  count_matrix[keep, , drop = FALSE],
  group = symbol[keep],
  reorder = FALSE
)
if (any(!is.finite(colSums(count_matrix))) || any(colSums(count_matrix) <= 0)) {
  stop("GSE281797 contains a non-positive library size after gene mapping")
}
cpm <- sweep(count_matrix, 2, colSums(count_matrix), "/") * 1e6
logcpm <- log2(cpm + 1)
selected_logcpm <- logcpm[, selected_samples_281797, drop = FALSE]

module_matrix <- selected_logcpm[
  intersect(p53_genes, rownames(selected_logcpm)),
  ,
  drop = FALSE
]
if (nrow(module_matrix) != length(p53_genes)) {
  stop(
    "GSE281797 is missing p53 module genes: ",
    paste(setdiff(p53_genes, rownames(selected_logcpm)), collapse = ", ")
  )
}
if (any(!is.finite(module_matrix))) {
  stop("GSE281797 p53 module matrix contains non-finite values")
}
module_nonzero <- colSums(module_matrix > 0) >= 2
if (any(!module_nonzero)) {
  stop(
    "GSE281797 samples with fewer than two detected p53-module genes: ",
    paste(colnames(module_matrix)[!module_nonzero], collapse = ", ")
  )
}
module_score <- as.numeric(GSVA::gsva(
  GSVA::ssgseaParam(
    selected_logcpm,
    list(p53_module = p53_genes),
    normalize = TRUE,
    alpha = 0.25
  ),
  verbose = FALSE
)["p53_module", ])

ductular_genes <- c("KRT7", "KRT19", "EPCAM", "SOX9")
if (!all(ductular_genes %in% rownames(selected_logcpm))) {
  stop("GSE281797 is missing one or more ductular-response genes")
}
ductular_matrix <- selected_logcpm[ductular_genes, , drop = FALSE]
ductular_score <- colMeans(t(scale(t(ductular_matrix))), na.rm = TRUE)

proportion <- read.csv(
  file.path(gse281797_dir, "GSE281797_music_proportions.csv"),
  check.names = FALSE,
  row.names = 1
)
hepatocyte_fraction <- proportion[
  match(colnames(selected_logcpm), rownames(proportion)),
  grep("Hepatocyte", colnames(proportion), ignore.case = TRUE)[1]
]

gse281797 <- data.frame(
  dataset = "GSE281797",
  gsm = colnames(selected_logcpm),
  group_unified = meta_281797$diagnosis[
    match(colnames(selected_logcpm), meta_281797$sample_id)
  ],
  nas = meta_281797$nas[
    match(colnames(selected_logcpm), meta_281797$sample_id)
  ],
  tp53i3 = as.numeric(selected_logcpm["TP53I3", ]),
  p53_ssgsea_z = as.numeric(scale(module_score)),
  age_years = meta_281797$age_years[
    match(colnames(selected_logcpm), meta_281797$sample_id)
  ],
  sex = factor(
    tolower(meta_281797$sex[
      match(colnames(selected_logcpm), meta_281797$sample_id)
    ]),
    levels = c("female", "male")
  ),
  bmi_kg_m2 = meta_281797$bmi_kg_m2[
    match(colnames(selected_logcpm), meta_281797$sample_id)
  ],
  diabetes = factor(
    ifelse(
      meta_281797$diabetes[
        match(colnames(selected_logcpm), meta_281797$sample_id)
      ] == 1,
      "yes",
      "no"
    ),
    levels = c("no", "yes")
  ),
  hepatocyte_fraction = as.numeric(hepatocyte_fraction),
  ductular_score = as.numeric(ductular_score),
  stringsAsFactors = FALSE
)
gse281797 <- gse281797[
  gse281797$group_unified %in% c("MASL", "MASH"),
]
gse281797$nas_bin <- as.integer(gse281797$nas >= 5)
gse281797$tp53i3_z <- as.numeric(scale(gse281797$tp53i3))
gse281797$ductular_z <- as.numeric(scale(gse281797$ductular_score))

common_columns <- c(
  "dataset", "gsm", "group_unified", "nas", "nas_bin", "tp53i3",
  "tp53i3_z", "p53_ssgsea_z", "age_years", "sex", "bmi_kg_m2",
  "diabetes", "alt_u_l", "ast_u_l", "hepatocyte_fraction",
  "ductular_score", "ductular_z"
)
for (column in setdiff(common_columns, names(data))) {
  data[[column]] <- NA
}
for (column in setdiff(common_columns, names(gse281797))) {
  gse281797[[column]] <- NA
}
combined <- rbind(
  data[, common_columns],
  gse281797[, common_columns]
)
combined$dataset <- factor(
  combined$dataset,
  levels = c("GSE135251", "GSE48452", "GSE89632", "GSE281797")
)

complete <- complete.cases(combined[, c(
  "nas_bin", "tp53i3_z", "p53_ssgsea_z", "hepatocyte_fraction",
  "ductular_z"
)])

preflight <- do.call(rbind, lapply(
  levels(combined$dataset),
  function(dataset) {
    x <- combined[combined$dataset == dataset, ]
    data.frame(
      dataset = dataset,
      records = nrow(x),
      complete_core = sum(complete[combined$dataset == dataset]),
      events_nas_ge_5 = sum(x$nas_bin == 1, na.rm = TRUE),
      non_events_nas_lt_5 = sum(x$nas_bin == 0, na.rm = TRUE),
      missing_nas = sum(is.na(x$nas_bin)),
      missing_tp53i3 = sum(is.na(x$tp53i3_z)),
      missing_module = sum(is.na(x$p53_ssgsea_z)),
      missing_hepatocyte = sum(is.na(x$hepatocyte_fraction)),
      missing_ductular = sum(is.na(x$ductular_z)),
      missing_age = sum(is.na(x$age_years)),
      missing_sex = sum(is.na(x$sex)),
      missing_bmi = sum(is.na(x$bmi_kg_m2)),
      missing_diabetes = sum(is.na(x$diabetes)),
      stringsAsFactors = FALSE
    )
  }
))

write.csv(
  combined,
  file.path(out_dir, "incremental_model_input_with_GSE281797.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  preflight,
  file.path(out_dir, "incremental_model_preflight.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  data.frame(
    sample = gse281797$gsm,
    diagnosis = gse281797$group_unified,
    nas = gse281797$nas,
    nas_bin = gse281797$nas_bin,
    tp53i3_log2cpm = gse281797$tp53i3,
    p53_ssgsea_z = gse281797$p53_ssgsea_z,
    hepatocyte_fraction = gse281797$hepatocyte_fraction,
    ductular_score = gse281797$ductular_score,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "GSE281797_incremental_input.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
