numeric_equal04 <- function(a, b, tolerance = 1e-11) {
  length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)) &&
    all(abs(a - b) <= tolerance * pmax(1, abs(a), abs(b)))
}
metrics04 <- function(y, predicted) {
  errors <- y - predicted
  c(rmse = sqrt(mean(errors^2)), mae = mean(abs(errors)),
    r2 = 1 - sum(errors^2) / sum((y - mean(y))^2), mse = mean(errors^2))
}
performance04 <- function(benchmarks, axes) {
  do.call(rbind, lapply(axes, function(a) {
    z <- benchmarks[benchmarks$axis == a, ]
    do.call(rbind, lapply(c("intercept_only", "arm_only", "full_ridge"), function(method) {
      m <- metrics04(z$observed_delta, z[[method]])
      data.frame(axis = a, method = method, n = nrow(z), oof_rmse = m[["rmse"]],
        oof_mae = m[["mae"]], oof_r2 = m[["r2"]], oof_mse = m[["mse"]])
    }))
  }))
}
residual_qc04 <- function(long, performance, axes) {
  do.call(rbind, lapply(axes, function(a) {
    z <- long[long$axis == a, ]
    r <- z$residual
    p <- performance[performance$axis == a, ]
    row <- data.frame(axis = a, n = nrow(z), missing = sum(is.na(r)),
      nonfinite = sum(!is.finite(r)), mean = mean(r), sd = sd(r), median = median(r),
      iqr = IQR(r), minimum = min(r), maximum = max(r), unique_values = length(unique(r)),
      sd_observed_delta = sd(z$observed_delta), sd_expected_delta = sd(z$expected_delta),
      cor_observed_expected = cor(z$observed_delta, z$expected_delta),
      cor_observed_residual = cor(z$observed_delta, r),
      cor_expected_residual = cor(z$expected_delta, r),
      full_ridge_rmse = p$oof_rmse[p$method == "full_ridge"],
      arm_only_rmse = p$oof_rmse[p$method == "arm_only"],
      intercept_only_rmse = p$oof_rmse[p$method == "intercept_only"])
    row$passed <- row$n == 150L && row$missing == 0L && row$nonfinite == 0L &&
      is.finite(row$sd) && row$sd > 0 && row$unique_values > 1L
    row
  }))
}
fold_diagnostics04 <- function(long, s) {
  do.call(rbind, lapply(s$axes, function(a) {
    do.call(rbind, lapply(seq_len(s$outer_folds), function(k) {
      z <- long[long$axis == a & long$outer_fold == k, ]
      m <- metrics04(z$observed_delta, z$expected_delta)
      data.frame(axis = a, outer_fold = k, n_test = nrow(z), test_rmse = m[["rmse"]],
        test_mae = m[["mae"]], mean_observed_delta = mean(z$observed_delta),
        mean_expected_delta = mean(z$expected_delta), mean_residual = mean(z$residual),
        sd_residual = sd(z$residual))
    }))
  }))
}
audit_saved04 <- function(x, s, folds) {
  long <- read_numeric04(stage04_long, c("patient_id", "arm", "axis"))
  wide <- read_numeric04(stage04_wide, c("patient_id", "arm"))
  scaling <- read_numeric04(stage04_scaling, "axis")
  bench <- read_numeric04("data/processed/residuals/neotrip_oof_benchmarks_v1.tsv",
    c("patient_id", "arm", "axis"))
  tuning <- read_numeric04("results/stage04_residuals/nested_tuning_results.tsv", c("axis", "selected_rule"))
  coefficients <- read_numeric04("results/stage04_residuals/fold_coefficients.tsv", c("axis", "term"))
  expected_long <- c("patient_id", "arm", "outer_fold", "axis", "baseline_score",
    "observed_delta", "expected_delta", "residual", "residual_z")
  expected_wide <- c("patient_id", "arm", "outer_fold", unlist(lapply(s$axes,
    function(a) paste0(c("observed_delta_", "expected_delta_", "residual_", "residual_z_"), a))))
  if (!identical(names(long), expected_long) || !identical(names(wide), expected_wide) ||
      nrow(long) != 900L || anyDuplicated(paste(long$patient_id, long$axis)) ||
      nrow(wide) != 150L || !identical(wide$patient_id, x$patient_id) ||
      !identical(wide$arm, x$arm) || !identical(as.integer(wide$outer_fold), folds$outer$outer_fold) ||
      !setequal(long$axis, s$axes) || !setequal(long$patient_id, x$patient_id) ||
      nrow(scaling) != 6L || !identical(scaling$axis, s$axes) || any(scaling$n != 150) ||
      any(!is.finite(as.matrix(long[, -(1:4)])))) stop04("Residual identity or schema invalid")
  if (nrow(bench) != 900L || anyDuplicated(paste(bench$patient_id, bench$axis)) ||
      nrow(tuning) != 60L || anyDuplicated(paste(tuning$axis, tuning$outer_fold)) ||
      nrow(coefficients) != 660L || anyDuplicated(paste(coefficients$axis, coefficients$outer_fold, coefficients$term)))
    stop04("Audit output identity invalid")
  design <- design04(x, s)
  separation <- list()
  for (k in seq_len(s$outer_folds)) {
    train_idx <- which(folds$outer$outer_fold != k)
    test_idx <- which(folds$outer$outer_fold == k)
    train <- x[train_idx, , drop = FALSE]
    test <- x[test_idx, , drop = FALSE]
    inner <- folds$inner[folds$inner$outer_fold == k, , drop = FALSE]
    inner <- inner[match(train$patient_id, inner$patient_id), , drop = FALSE]
    for (a in s$axes) {
      model <- readRDS(model_path04(a, k))
      validate_model04(model, s, train, test, inner, a, k)
      z <- long[long$outer_fold == k & long$axis == a, ]
      counts <- table(factor(z$patient_id, levels = test$patient_id))
      separation[[length(separation) + 1L]] <- data.frame(axis = a, outer_fold = k, n_test = nrow(test),
        test_in_training = length(intersect(test$patient_id, model$train_patient_ids)),
        test_in_inner = length(intersect(test$patient_id, model$inner_assignments$patient_id)),
        oof_count_violations = sum(counts != 1L) + sum(!z$patient_id %in% test$patient_id))
      if (!setequal(z$patient_id, test$patient_id)) stop04("OOF patient assignment failed")
      z <- z[match(test$patient_id, z$patient_id), ]
      pred <- as.numeric(predict(model$fit, newx = design[test_idx, , drop = FALSE], s = s$lambda_rule))
      if (!identical(z$arm, test$arm) || !numeric_equal04(z$expected_delta, pred) ||
          !numeric_equal04(z$baseline_score, test[[paste0("baseline_", a)]]) ||
          !numeric_equal04(z$observed_delta, test[[paste0("delta_", a)]]) ||
          !numeric_equal04(z$residual, z$observed_delta - z$expected_delta))
        stop04("Saved OOF predictions differ")
      b <- bench[bench$axis == a & bench$outer_fold == k, ]
      b <- b[match(test$patient_id, b$patient_id), ]
      arm_mean <- vapply(test$arm, function(arm) mean(train[[paste0("delta_", a)]][train$arm == arm]), numeric(1))
      if (anyNA(b$patient_id) || !identical(b$arm, test$arm) ||
          !numeric_equal04(b$observed_delta, z$observed_delta) ||
          !numeric_equal04(b$full_ridge, z$expected_delta) ||
          !numeric_equal04(b$intercept_only, rep(mean(train[[paste0("delta_", a)]]), nrow(test))) ||
          !numeric_equal04(b$arm_only, arm_mean)) stop04("OOF benchmark validation failed")
      t <- tuning[tuning$axis == a & tuning$outer_fold == k, ]
      if (nrow(t) != 1L || t$n_train != nrow(train) || t$n_test != nrow(test) ||
          t$selected_rule != "lambda.1se" ||
          !numeric_equal04(c(t$lambda_min, t$lambda_1se, t$inner_cv_mse_min, t$selected_lambda),
            c(model$fit$lambda.min, model$fit$lambda.1se, min(model$fit$cvm), model$fit$lambda.1se)))
        stop04("Tuning audit differs")
      co <- as.matrix(coef(model$fit, s = s$lambda_rule))
      saved <- coefficients[coefficients$axis == a & coefficients$outer_fold == k, ]
      if (!identical(saved$term, rownames(co)) || !numeric_equal04(saved$coefficient, as.numeric(co)) ||
          !numeric_equal04(saved$lambda, rep(model$fit$lambda.1se, nrow(co))))
        stop04("Coefficient audit differs")
    }
  }
  separation <- do.call(rbind, separation)
  write_tsv(separation, path("results/stage04_residuals/patient_separation_check.tsv"))
  if (any(as.matrix(separation[, -(1:3)]) != 0)) stop04("Patient separation violated")
  for (a in s$axes) {
    z <- long[long$axis == a, ]
    z <- z[match(x$patient_id, z$patient_id), ]
    sc <- scaling[scaling$axis == a, ]
    if (!numeric_equal04(c(sc$mean, sc$sd), c(mean(z$residual), sd(z$residual))) ||
        !numeric_equal04(z$residual_z, (z$residual - sc$mean) / sc$sd))
      stop04("Residual scaling differs")
    for (f in c("observed_delta", "expected_delta", "residual", "residual_z"))
      if (!numeric_equal04(wide[[paste0(f, "_", a)]], z[[f]])) stop04("Wide residual values differ")
  }
  list(long = long, wide = wide, scaling = scaling, benchmarks = bench, separation = separation)
}
all_names04 <- function(x) {
  unique(c(names(x), if (is.list(x) && !is.data.frame(x)) unlist(lapply(x, all_names04), use.names = FALSE)))
}
leakage04 <- function(s) {
  files <- unique(c(stage04_input, stage04_outer, stage04_inner, stage04_long, stage04_wide, stage04_scaling,
    "data/processed/residuals/neotrip_oof_benchmarks_v1.tsv",
    substring(list.files(path("results/stage04_residuals"), pattern = "[.]tsv$", full.names = TRUE),
      nchar(root) + 2L)))
  files <- files[!basename(files) %in% c("leakage_check.tsv", "stage04_output_checksums.tsv")]
  reports <- lapply(files, function(f) {
    fields <- names(read.delim(path(f), nrows = 0L, check.names = FALSE))
    data.frame(file = f, fields_scanned = length(fields), forbidden_fields = length(forbidden04(fields)),
      permitted_metadata = "", passed = !length(forbidden04(fields)))
  })
  fields <- c(setdiff(all_names04(s), c("outcomes", "outcome_information_used")), s$predictors, s$outcomes)
  permitted <- identical(s$outcome_information_used, FALSE) && identical(s$outcomes, paste0("delta_", s$axes))
  reports[[length(reports) + 1L]] <- data.frame(file = stage04_spec, fields_scanned = length(fields),
    forbidden_fields = length(forbidden04(fields)) + as.integer(!permitted),
    permitted_metadata = "outcomes: six prespecified molecular deltas; outcome_information_used: false (required protocol declarations)",
    passed = !length(forbidden04(fields)) && permitted)
  for (file in list.files(path("models/stage04"), pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)) {
    model <- readRDS(file)
    fields <- unique(c(all_names04(model), model$predictors, model$target, rownames(model$fit$glmnet.fit$beta)))
    reports[[length(reports) + 1L]] <- data.frame(file = substring(file, nchar(root) + 2L),
      fields_scanned = length(fields), forbidden_fields = length(forbidden04(fields)),
      permitted_metadata = "", passed = !length(forbidden04(fields)))
  }
  if (file.exists(path(stage04_freeze))) {
    f <- yaml::read_yaml(path(stage04_freeze))
    fields <- setdiff(names(f), "outcome_information_used")
    permitted <- identical(f$outcome_information_used, FALSE)
    reports[[length(reports) + 1L]] <- data.frame(file = stage04_freeze, fields_scanned = length(fields),
      forbidden_fields = length(forbidden04(fields)) + as.integer(!permitted),
      permitted_metadata = "outcome_information_used: false (required protocol declaration)",
      passed = !length(forbidden04(fields)) && permitted)
  }
  report <- do.call(rbind, reports)
  write_tsv(report, path("results/stage04_residuals/leakage_check.tsv"))
  if (any(!report$passed)) stop04("Outcome leakage detected")
  report
}
figures04 <- function(long, diagnostics, s) {
  output <- function(name, draw) {
    pdf(path("figures/stage04_residuals", paste0(name, ".pdf")), width = 12, height = 7)
    tryCatch(draw(), finally = dev.off())
    png(path("figures/stage04_residuals", paste0(name, ".png")), width = 2400, height = 1400, res = 200)
    tryCatch(draw(), finally = dev.off())
  }
  layout <- function(title) {
    par(mfrow = c(2, 3), mar = c(4.3, 4.5, 3, 1), oma = c(0, 0, 2, 0))
  }
  output("observed_expected_v1", function() {
    layout("Observed and cross-fitted expected molecular change")
    for (a in s$axes) {
      z <- long[long$axis == a, ]
      limits <- range(c(z$observed_delta, z$expected_delta))
      plot(z$expected_delta, z$observed_delta, main = gsub("_", " ", a),
        xlab = "Cross-fitted expected delta", ylab = "Observed delta",
        xlim = limits, ylim = limits, pch = 16, col = "#276A9890", cex = 0.65)
      abline(a = 0, b = 1, lty = 2, col = "grey50")
    }
    mtext("Observed and cross-fitted expected molecular change", outer = TRUE, cex = 1.1)
  })
  output("residual_distributions_v1", function() {
    layout("Cross-fitted residual distributions")
    for (a in s$axes) {
      r <- long$residual[long$axis == a]
      hist(r, breaks = "FD", main = gsub("_", " ", a),
        xlab = "Observed minus expected delta", ylab = "Patients", col = "#AACDE5", border = "white")
      abline(v = 0, col = "grey40", lty = 2)
    }
    mtext("Cross-fitted residual distributions", outer = TRUE, cex = 1.1)
  })
  output("outer_fold_errors_v1", function() {
    layout("Outer-fold prediction error")
    for (a in s$axes) {
      z <- diagnostics[diagnostics$axis == a, ]
      plot(z$outer_fold, z$test_rmse, type = "b", pch = 16, col = "#276A98",
        main = gsub("_", " ", a), xlab = "Outer fold", ylab = "Test RMSE",
        xaxt = "n", ylim = c(0, max(z$test_rmse) * 1.1))
      axis(1, at = seq_len(s$outer_folds))
    }
    mtext("Outer-fold prediction error", outer = TRUE, cex = 1.1)
  })
}
