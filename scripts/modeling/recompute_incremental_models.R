suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cv_seed <- 20261018L
bootstrap_seed <- 20261019L
bootstrap_B <- 500L

data <- read.csv(
  input_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
data$sex <- factor(data$sex, levels = c("female", "male"))
data$diabetes <- factor(data$diabetes, levels = c("no", "yes"))

configs <- list(
  list(
    id = "GSE89632_LOO",
    dataset = "GSE89632",
    role = "primary",
    m0 = "age_years + sex + bmi_kg_m2 + diabetes + alt_u_l + ast_u_l",
    cv = "LOO",
    K = NA_integer_,
    R = NA_integer_,
    seed = NA_integer_
  ),
  list(
    id = "GSE48452_5x20",
    dataset = "GSE48452",
    role = "primary",
    m0 = "age_years + sex + bmi_kg_m2",
    cv = "5x20",
    K = 5L,
    R = 20L,
    seed = cv_seed
  ),
  list(
    id = "GSE135251_5x20",
    dataset = "GSE135251",
    role = "supplementary",
    m0 = "1",
    cv = "5x20",
    K = 5L,
    R = 20L,
    seed = cv_seed
  ),
  list(
    id = "GSE281797_5x20",
    dataset = "GSE281797",
    role = "primary_event_limited",
    m0 = "age_years + sex + bmi_kg_m2 + diabetes",
    cv = "5x20",
    K = 5L,
    R = 20L,
    seed = cv_seed
  )
)

models <- list(
  M0 = function(m0) m0,
  M1 = function(m0) paste(m0, "+ tp53i3_z"),
  M2 = function(m0) paste(m0, "+ p53_ssgsea_z"),
  M12 = function(m0) paste(m0, "+ tp53i3_z + p53_ssgsea_z"),
  M1c = function(m0) paste(
    m0,
    "+ tp53i3_z + ductular_z + hepatocyte_z"
  ),
  M2c = function(m0) paste(
    m0,
    "+ p53_ssgsea_z + ductular_z + hepatocyte_z"
  ),
  M12c = function(m0) paste(
    m0,
    "+ tp53i3_z + p53_ssgsea_z + ductular_z + hepatocyte_z"
  )
)

for (dataset in unique(data$dataset)) {
  rows <- data$dataset == dataset
  data$tp53i3_z[rows] <- as.numeric(scale(data$tp53i3[rows]))
  data$p53_ssgsea_z[rows] <- as.numeric(scale(data$p53_ssgsea_z[rows]))
  data$ductular_z[rows] <- as.numeric(scale(data$ductular_score[rows]))
  data$hepatocyte_z[rows] <- as.numeric(scale(data$hepatocyte_fraction[rows]))
}

standardize_train_test <- function(train, test, variable_names) {
  for (variable in variable_names) {
    center <- mean(train[[variable]], na.rm = TRUE)
    spread <- sd(train[[variable]], na.rm = TRUE)
    if (!is.finite(spread) || spread == 0) spread <- 1
    train[[variable]] <- (train[[variable]] - center) / spread
    test[[variable]] <- (test[[variable]] - center) / spread
  }
  list(train = train, test = test)
}

fit_predict <- function(train, test, m0) {
  continuous <- c(
    "tp53i3_z", "p53_ssgsea_z", "ductular_z", "hepatocyte_z",
    "age_years", "bmi_kg_m2", "alt_u_l", "ast_u_l"
  )
  standardized <- standardize_train_test(train, test, continuous)
  train <- standardized$train
  test <- standardized$test

  output <- list()
  fit_status <- character(length(models))
  names(fit_status) <- names(models)
  for (model_name in names(models)) {
    formula <- as.formula(paste("nas_bin ~", models[[model_name]](m0)))
    fit <- try(
      suppressWarnings(glm(formula, data = train, family = binomial)),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      output[[model_name]] <- rep(NA_real_, nrow(test))
      fit_status[model_name] <- paste0(
        "fit_error: ",
        conditionMessage(attr(fit, "condition"))
      )
    } else {
      prediction <- try(
        predict(fit, newdata = test, type = "response"),
        silent = TRUE
      )
      if (inherits(prediction, "try-error")) {
        output[[model_name]] <- rep(NA_real_, nrow(test))
        fit_status[model_name] <- paste0(
          "prediction_error: ",
          conditionMessage(attr(prediction, "condition"))
        )
      } else {
        output[[model_name]] <- as.numeric(prediction)
        fit_status[model_name] <- if (isTRUE(fit$converged)) {
          "ok"
        } else {
          "not_converged_used"
        }
      }
    }
  }
  list(predictions = output, status = fit_status)
}

stratified_folds <- function(outcome, K) {
  folds <- integer(length(outcome))
  for (level in sort(unique(outcome))) {
    index <- which(outcome == level)
    index <- sample(index)
    folds[index] <- rep(seq_len(K), length.out = length(index))
  }
  folds
}

auc_value <- function(y, prediction) {
  valid <- !is.na(prediction)
  if (sum(valid) < 8 || length(unique(y[valid])) < 2) {
    return(NA_real_)
  }
  as.numeric(auc(roc(
    y[valid],
    prediction[valid],
    quiet = TRUE,
    direction = "<"
  )))
}

bootstrap_delta <- function(y, p0, p1, B, seed) {
  valid <- !is.na(p0) & !is.na(p1)
  y <- y[valid]
  p0 <- p0[valid]
  p1 <- p1[valid]
  point <- auc_value(y, p1) - auc_value(y, p0)
  controls <- which(y == 0)
  cases <- which(y == 1)
  set.seed(seed)
  draws <- numeric(0)
  for (index in seq_len(B)) {
    rows <- c(
      sample(controls, length(controls), replace = TRUE),
      sample(cases, length(cases), replace = TRUE)
    )
    auc0 <- auc_value(y[rows], p0[rows])
    auc1 <- auc_value(y[rows], p1[rows])
    if (is.finite(auc0) && is.finite(auc1)) {
      draws <- c(draws, auc1 - auc0)
    }
  }
  if (length(draws) < 50) {
    return(list(
      result = c(delta = point, lcl = NA_real_, ucl = NA_real_),
      draws = draws
    ))
  }
  interval <- quantile(draws, c(0.025, 0.975), names = FALSE)
  list(
    result = c(delta = point, lcl = interval[1], ucl = interval[2]),
    draws = draws
  )
}

comparisons <- list(
  "M12-M1" = c("M1", "M12"),
  "M12-M2" = c("M2", "M12"),
  "M2-M1" = c("M1", "M2"),
  "M2c-M1c" = c("M1c", "M2c"),
  "M12c-M1c" = c("M1c", "M12c")
)

all_delta <- list()
all_draws <- list()
all_auc <- list()
all_predictions <- list()
all_calibration <- list()
all_dca <- list()
fit_diagnostics <- list()

for (config in configs) {
  model_variables <- unique(c(
    "nas_bin", "tp53i3_z", "p53_ssgsea_z", "ductular_z", "hepatocyte_z",
    all.vars(as.formula(paste("nas_bin ~", config$m0)))
  ))
  x <- data[
    data$dataset == config$dataset &
      complete.cases(data[, model_variables, drop = FALSE]),
    ,
    drop = FALSE
  ]
  y <- x$nas_bin

  if (config$cv == "LOO") {
    prediction_matrix <- matrix(
      NA_real_,
      nrow = nrow(x),
      ncol = length(models),
      dimnames = list(NULL, names(models))
    )
    status_counts <- setNames(
      rep(0L, length(models)),
      names(models)
    )
    observed_statuses <- character(0)
    names(status_counts) <- names(models)
    for (row_index in seq_len(nrow(x))) {
      prediction <- fit_predict(
        x[-row_index, , drop = FALSE],
        x[row_index, , drop = FALSE],
        config$m0
      )
      observed_statuses <- c(observed_statuses, prediction$status)
      prediction_matrix[row_index, ] <- unlist(prediction$predictions)
      status_counts <- status_counts +
        as.integer(
          prediction$status %in% c("ok", "not_converged_used")
        )
    }
  } else {
    set.seed(config$seed)
    prediction_sum <- matrix(
      0,
      nrow = nrow(x),
      ncol = length(models),
      dimnames = list(NULL, names(models))
    )
    prediction_count <- matrix(
      0,
      nrow = nrow(x),
      ncol = length(models),
      dimnames = list(NULL, names(models))
    )
    status_counts <- setNames(
      rep(0L, length(models)),
      names(models)
    )
    observed_statuses <- character(0)
    for (repeat_index in seq_len(config$R)) {
      folds <- stratified_folds(y, config$K)
      for (fold in seq_len(config$K)) {
        test_rows <- which(folds == fold)
        prediction <- fit_predict(
          x[-test_rows, , drop = FALSE],
          x[test_rows, , drop = FALSE],
          config$m0
        )
        observed_statuses <- c(observed_statuses, prediction$status)
        valid <- prediction$status %in% c("ok", "not_converged_used")
        prediction_sum[test_rows, valid] <-
          prediction_sum[test_rows, valid, drop = FALSE] +
          matrix(
            unlist(prediction$predictions[valid]),
            nrow = length(test_rows),
            ncol = sum(valid)
          )
        prediction_count[test_rows, valid] <-
          prediction_count[test_rows, valid, drop = FALSE] + 1
        status_counts <- status_counts +
          as.integer(
            prediction$status %in% c("ok", "not_converged_used")
          ) *
          length(test_rows)
      }
    }
    prediction_matrix <- prediction_sum / prediction_count
    prediction_matrix[!is.finite(prediction_matrix)] <- NA_real_
  }

  fit_diagnostics[[config$id]] <- data.frame(
    analysis = config$id,
    dataset = config$dataset,
    n = nrow(x),
    events = sum(y == 1),
    model = names(models),
    successful_fold_predictions = as.integer(status_counts[names(models)]),
    expected_fold_predictions = if (config$cv == "LOO") {
      nrow(x)
    } else {
      config$R * nrow(x)
    },
    observed_statuses = paste(
      sort(unique(observed_statuses)),
      collapse = " | "
    ),
    stringsAsFactors = FALSE
  )

  prediction_table <- data.frame(
    analysis = config$id,
    dataset = config$dataset,
    gsm = x$gsm,
    y = y,
    prediction_matrix,
    check.names = FALSE
  )
  all_predictions[[config$id]] <- prediction_table

  for (model_name in names(models)) {
    all_auc[[length(all_auc) + 1]] <- data.frame(
      analysis = config$id,
      dataset = config$dataset,
      role = config$role,
      cv = config$cv,
      model = model_name,
      n = nrow(x),
      events = sum(y == 1),
      auc = auc_value(y, prediction_matrix[, model_name]),
      stringsAsFactors = FALSE
    )
  }

  for (comparison_index in seq_along(comparisons)) {
    comparison_name <- names(comparisons)[comparison_index]
    pair <- comparisons[[comparison_name]]
    boot <- bootstrap_delta(
      y,
      prediction_matrix[, pair[1]],
      prediction_matrix[, pair[2]],
      bootstrap_B,
      bootstrap_seed + comparison_index * 100L
    )
    half_width <- (boot$result[["ucl"]] - boot$result[["lcl"]]) / 2
    approximate_se <- half_width / 1.96
    minimum_detectable_effect <- 2.80 * approximate_se
    all_delta[[length(all_delta) + 1]] <- data.frame(
      analysis = config$id,
      dataset = config$dataset,
      role = config$role,
      cv = config$cv,
      n = nrow(x),
      events = sum(y == 1),
      comparison = comparison_name,
      model_00 = pair[1],
      model_01 = pair[2],
      delta_auc = boot$result[["delta"]],
      lcl_95 = boot$result[["lcl"]],
      ucl_95 = boot$result[["ucl"]],
      approximate_mde = minimum_detectable_effect,
      bootstrap_B = bootstrap_B,
      successful_bootstrap_draws = length(boot$draws),
      cv_seed = config$seed,
      bootstrap_seed = bootstrap_seed + comparison_index * 100L,
      stringsAsFactors = FALSE
    )
    all_draws[[length(all_draws) + 1]] <- data.frame(
      analysis = rep(config$id, length(boot$draws)),
      comparison = rep(comparison_name, length(boot$draws)),
      draw = seq_along(boot$draws),
      delta_auc = boot$draws,
      stringsAsFactors = FALSE
    )
  }

  for (model_name in c("M1", "M12")) {
    prediction <- pmin(
      pmax(prediction_matrix[, model_name], 1e-6),
      1 - 1e-6
    )
    if (sum(is.finite(prediction)) >= 8) {
      calibration_fit <- glm(
        y ~ qlogis(prediction),
        family = binomial
      )
      coefficient <- coef(summary(calibration_fit))
      all_calibration[[length(all_calibration) + 1]] <- data.frame(
        analysis = config$id,
        dataset = config$dataset,
        model = model_name,
        n = nrow(x),
        calibration_slope = coefficient["qlogis(prediction)", "Estimate"],
        calibration_intercept = coefficient["(Intercept)", "Estimate"],
        stringsAsFactors = FALSE
      )
    }
  }

  thresholds <- seq(0.05, 0.80, by = 0.01)
  prevalence <- mean(y)
  for (threshold in thresholds) {
    all_dca[[length(all_dca) + 1]] <- data.frame(
      analysis = config$id,
      dataset = config$dataset,
      model = "treat_all",
      threshold = threshold,
      net_benefit = prevalence -
        (1 - prevalence) * threshold / (1 - threshold),
      stringsAsFactors = FALSE
    )
    for (model_name in c("M1", "M12")) {
      prediction <- prediction_matrix[, model_name]
      true_positive <- mean(prediction >= threshold & y == 1, na.rm = TRUE)
      false_positive <- mean(prediction >= threshold & y == 0, na.rm = TRUE)
      net_benefit <- true_positive -
        false_positive * threshold / (1 - threshold)
      all_dca[[length(all_dca) + 1]] <- data.frame(
        analysis = config$id,
        dataset = config$dataset,
        model = model_name,
        threshold = threshold,
        net_benefit = net_benefit,
        stringsAsFactors = FALSE
      )
    }
  }
}

delta_table <- do.call(rbind, all_delta)
auc_table <- do.call(rbind, all_auc)
prediction_table <- do.call(rbind, all_predictions)
calibration_table <- if (length(all_calibration)) {
  do.call(rbind, all_calibration)
} else {
  data.frame()
}
dca_table <- do.call(rbind, all_dca)
fit_diagnostics_table <- do.call(rbind, fit_diagnostics)

write.csv(
  delta_table,
  file.path(out_dir, "incremental_delta_auc.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  do.call(rbind, all_draws),
  file.path(out_dir, "incremental_bootstrap_draws.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  auc_table,
  file.path(out_dir, "incremental_auc.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  prediction_table,
  file.path(out_dir, "incremental_predictions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  calibration_table,
  file.path(out_dir, "incremental_calibration.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  dca_table,
  file.path(out_dir, "incremental_dca.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  fit_diagnostics_table,
  file.path(out_dir, "incremental_fit_diagnostics.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

plot_delta <- delta_table[
  delta_table$comparison %in% c("M12-M1", "M2-M1", "M12c-M1c"),
]
plot_delta$cohort <- factor(
  plot_delta$dataset,
  levels = c("GSE135251", "GSE48452", "GSE89632", "GSE281797")
)
plot_delta$comparison_label <- factor(
  plot_delta$comparison,
  levels = c("M12-M1", "M2-M1", "M12c-M1c"),
  labels = c("M12 - M1", "M2 - M1", "M12c - M1c")
)

panel_a <- ggplot(
  plot_delta,
  aes(x = delta_auc, y = cohort, color = comparison_label)
) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_pointrange(
    aes(xmin = lcl_95, xmax = ucl_95),
    position = position_dodge(width = 0.62),
    linewidth = 0.35,
    size = 0.28
  ) +
  scale_color_manual(
    values = c("#B2182B", "#2166AC", "#7F7F7F"),
    name = NULL
  ) +
  labs(
    title = "A  Incremental value (bootstrap)",
    x = "Delta AUC (95% CI)",
    y = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7.2),
    plot.title = element_text(face = "bold", size = 8)
  )

roc_analyses <- c("GSE48452_5x20", "GSE89632_LOO", "GSE281797_5x20")
roc_rows <- list()
for (analysis_name in roc_analyses) {
  x <- prediction_table[prediction_table$analysis == analysis_name, ]
  for (model_name in c("M1", "M12")) {
    valid <- is.finite(x[[model_name]])
    if (sum(valid) < 8 || length(unique(x$y[valid])) < 2) next
    curve <- roc(
      x$y[valid],
      x[[model_name]][valid],
      quiet = TRUE,
      direction = "<"
    )
    roc_rows[[length(roc_rows) + 1]] <- data.frame(
      analysis = analysis_name,
      model = model_name,
      specificity = rev(curve$specificities),
      sensitivity = rev(curve$sensitivities)
    )
  }
}
roc_plot <- do.call(rbind, roc_rows)

panel_b <- ggplot(
  roc_plot,
  aes(x = 1 - specificity, y = sensitivity, color = model)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey75") +
  geom_path(linewidth = 0.45) +
  facet_wrap(
    ~analysis,
    labeller = as_labeller(c(
      GSE48452_5x20 = "GSE48452 (5x20)",
      GSE89632_LOO = "GSE89632 (LOO)",
      GSE281797_5x20 = "GSE281797 (5 events)"
    ))
  ) +
  scale_color_manual(values = c(M1 = "#2166AC", M12 = "#B2182B")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "B  ROC (out-of-fold)",
    x = "1 - specificity",
    y = "sensitivity",
    color = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7.2),
    panel.spacing.x = unit(0.45, "lines"),
    strip.text = element_text(size = 7),
    plot.title = element_text(face = "bold", size = 8)
  )

calibration_rows <- list()
for (analysis_name in roc_analyses) {
  x <- prediction_table[prediction_table$analysis == analysis_name, ]
  for (model_name in c("M1", "M12")) {
    prediction <- pmin(pmax(x[[model_name]], 1e-6), 1 - 1e-6)
    valid <- is.finite(prediction)
    prediction <- prediction[valid]
    outcome <- x$y[valid]
    breaks <- unique(quantile(prediction, probs = seq(0, 1, 0.2)))
    if (length(breaks) < 3) next
    bins <- cut(
      prediction,
      breaks = breaks,
      include.lowest = TRUE,
      labels = FALSE
    )
    for (bin in sort(unique(bins))) {
      rows <- which(bins == bin)
      calibration_rows[[length(calibration_rows) + 1]] <- data.frame(
        analysis = analysis_name,
        model = model_name,
        predicted = mean(prediction[rows]),
        observed = mean(outcome[rows]),
        bin = bin
      )
    }
  }
}
calibration_plot <- do.call(rbind, calibration_rows)

panel_c <- ggplot(
  calibration_plot,
  aes(x = predicted, y = observed, color = model)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey75") +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.1) +
  facet_wrap(
    ~analysis,
    labeller = as_labeller(c(
      GSE48452_5x20 = "GSE48452",
      GSE89632_LOO = "GSE89632",
      GSE281797_5x20 = "GSE281797"
    ))
  ) +
  scale_color_manual(values = c(M1 = "#2166AC", M12 = "#B2182B")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "C  Calibration (quintiles)",
    x = "Predicted probability",
    y = "Observed proportion",
    color = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7.2),
    panel.spacing.x = unit(0.45, "lines"),
    strip.text = element_text(size = 7),
    plot.title = element_text(face = "bold", size = 8)
  )

dca_plot <- dca_table[dca_table$analysis %in% roc_analyses, ]
panel_d <- ggplot(
  dca_plot,
  aes(x = threshold, y = net_benefit, color = model)
) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey55") +
  geom_line(linewidth = 0.45) +
  facet_wrap(
    ~analysis,
    labeller = as_labeller(c(
      GSE48452_5x20 = "GSE48452",
      GSE89632_LOO = "GSE89632",
      GSE281797_5x20 = "GSE281797"
    ))
  ) +
  scale_color_manual(
    values = c(
      M1 = "#2166AC",
      M12 = "#B2182B",
      treat_all = "#777777"
    )
  ) +
  labs(
    title = "D  Decision curve",
    x = "Threshold probability",
    y = "Net benefit",
    color = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7.2),
    panel.spacing.x = unit(0.45, "lines"),
    strip.text = element_text(size = 7),
    plot.title = element_text(face = "bold", size = 8)
  )

figure <- arrangeGrob(
  panel_a,
  panel_b,
  panel_c,
  panel_d,
  ncol = 2,
  heights = c(1.05, 1)
)
png(
  file.path(out_dir, "Figure5_incremental.png"),
  width = 7.4 * 600,
  height = 5.4 * 600,
  res = 600,
  bg = "white",
  type = "cairo"
)
grid::grid.draw(figure)
dev.off()
ggsave(
  file.path(out_dir, "Figure5_incremental.pdf"),
  figure,
  width = 7.4,
  height = 5.4,
  bg = "white"
)
ggsave(
  file.path(out_dir, "Figure5_incremental.svg"),
  figure,
  width = 7.4,
  height = 5.4,
  bg = "white"
)

writeLines(
  c(
    paste("R:", R.version.string),
    paste("pROC:", as.character(packageVersion("pROC"))),
    paste("input:", input_path),
    paste("cv_seed:", cv_seed),
    paste("bootstrap_seed:", bootstrap_seed),
    paste("bootstrap_B:", bootstrap_B),
    "GSE281797: 5-fold repeated 20 times; NAS >= 5 events = 5",
    "GSE281797 is event-limited and remains exploratory.",
    paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    capture.output(sessionInfo())
  ),
  file.path(out_dir, "run_metadata.txt"),
  useBytes = TRUE
)
