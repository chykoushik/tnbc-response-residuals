design04 <- function(x, s) {
  if (!setequal(x$arm, names(s$arm_encoding))) stop04("Unexpected treatment arm")
  z <- x[, s$predictors, drop = FALSE]
  z$arm <- unname(unlist(s$arm_encoding)[x$arm])
  m <- as.matrix(z)
  storage.mode(m) <- "double"
  rownames(m) <- x$patient_id
  assert_fields04(colnames(m))
  if (any(!is.finite(m)) || ncol(m) != 10L) stop04("Predictor matrix invalid")
  m
}
model_path04 <- function(axis, fold) {
  path("models/stage04", axis, sprintf("outer_fold_%02d.rds", fold))
}
validate_model04 <- function(model, s, train, test, inner, axis, k) {
  if (!identical(model$axis, axis) || model$outer_fold != k ||
      !identical(model$train_patient_ids, train$patient_id) ||
      !identical(model$test_patient_ids, test$patient_id) ||
      length(intersect(model$train_patient_ids, model$test_patient_ids)) ||
      !identical(model$inner_assignments$patient_id, train$patient_id) ||
      !identical(as.integer(model$inner_assignments$inner_fold), as.integer(inner$inner_fold)) ||
      !identical(model$predictors, s$predictors) ||
      !identical(model$target, paste0("delta_", axis)) ||
      !identical(model$input_sha256, sha256(path(stage04_input))) ||
      !identical(model$model_spec_sha256, sha256(path(stage04_spec))) ||
      !identical(model$glmnet_version, s$glmnet_version) ||
      !identical(as.integer(model$fit$foldid), as.integer(inner$inner_fold)) ||
      model$fit$glmnet.fit$jerr != 0L) stop04("Saved training model invalid")
  assert_fields04(c(names(model), model$predictors, model$target, names(model$inner_assignments)))
  if (!is.finite(model$fit$lambda.1se) || model$fit$lambda.1se <= 0 ||
      !is.finite(model$fit$lambda.min)) stop04("Tuning selection invalid")
  invisible(TRUE)
}
fit_models04 <- function(x, s, folds) {
  design <- design04(x, s)
  rows <- tuning <- coefficients <- benchmark_rows <- list()
  training_warnings <- data.frame(axis = character(), outer_fold = integer(), warning = character())
  for (k in seq_len(s$outer_folds)) {
    test_idx <- which(folds$outer$outer_fold == k)
    train_idx <- which(folds$outer$outer_fold != k)
    train <- x[train_idx, , drop = FALSE]
    test <- x[test_idx, , drop = FALSE]
    inner <- folds$inner[folds$inner$outer_fold == k, , drop = FALSE]
    inner <- inner[match(train$patient_id, inner$patient_id), , drop = FALSE]
    for (a in s$axes) {
      file <- model_path04(a, k)
      y <- train[[paste0("delta_", a)]]
      if (file.exists(file)) {
        model <- readRDS(file)
      } else {
        fit <- withCallingHandlers(glmnet::cv.glmnet(x = design[train_idx, , drop = FALSE], y = y,
          family = s$model_family, alpha = s$alpha, intercept = s$intercept,
          standardize = s$standardize, type.measure = s$type_measure,
          nlambda = s$nlambda, lambda.min.ratio = s$lambda_min_ratio,
          nfolds = s$inner_folds, foldid = inner$inner_fold,
          grouped = s$cv_grouped, alignment = s$cv_alignment, parallel = s$cv_parallel, keep = TRUE),
          warning = function(w) {
            training_warnings <<- rbind(training_warnings,
              data.frame(axis = a, outer_fold = k, warning = conditionMessage(w)))
            invokeRestart("muffleWarning")
          })
        model <- list(axis = a, outer_fold = k, train_patient_ids = train$patient_id,
          test_patient_ids = test$patient_id, inner_assignments = inner,
          predictors = s$predictors, target = paste0("delta_", a),
          input_sha256 = sha256(path(stage04_input)), model_spec_sha256 = sha256(path(stage04_spec)),
          glmnet_version = s$glmnet_version, fit = fit)
        validate_model04(model, s, train, test, inner, a, k)
        dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
        saveRDS(model, file)
      }
      validate_model04(model, s, train, test, inner, a, k)
      fit <- model$fit
      predicted <- as.numeric(predict(fit, newx = design[test_idx, , drop = FALSE], s = s$lambda_rule))
      observed <- test[[paste0("delta_", a)]]
      rows[[length(rows) + 1L]] <- data.frame(patient_id = test$patient_id, arm = test$arm,
        outer_fold = k, axis = a, baseline_score = test[[paste0("baseline_", a)]],
        observed_delta = observed, expected_delta = predicted, residual = observed - predicted)
      tuning[[length(tuning) + 1L]] <- data.frame(axis = a, outer_fold = k,
        n_train = nrow(train), n_test = nrow(test), lambda_min = fit$lambda.min,
        lambda_1se = fit$lambda.1se, inner_cv_mse_min = min(fit$cvm),
        selected_lambda = fit$lambda.1se, selected_rule = s$lambda_rule)
      co <- as.matrix(coef(fit, s = s$lambda_rule))
      coefficients[[length(coefficients) + 1L]] <- data.frame(axis = a, outer_fold = k,
        term = rownames(co), coefficient = as.numeric(co), lambda = fit$lambda.1se)
      arm_train <- data.frame(delta = y, arm = unname(unlist(s$arm_encoding)[train$arm]))
      arm_fit <- lm(delta ~ arm, data = arm_train)
      arm_pred <- as.numeric(predict(arm_fit,
        newdata = data.frame(arm = unname(unlist(s$arm_encoding)[test$arm]))))
      benchmark_rows[[length(benchmark_rows) + 1L]] <- data.frame(patient_id = test$patient_id,
        arm = test$arm, outer_fold = k, axis = a, observed_delta = observed,
        intercept_only = mean(y), arm_only = arm_pred, full_ridge = predicted)
    }
  }
  write_tsv(training_warnings, path("results/stage04_residuals/training_warnings.tsv"))
  if (nrow(training_warnings)) stop04("Model training warnings require review", "MODIFY_REVIEW")
  long <- do.call(rbind, rows)
  long <- long[order(match(long$patient_id, x$patient_id), match(long$axis, s$axes)), , drop = FALSE]
  scaling <- do.call(rbind, lapply(s$axes, function(a) {
    r <- long$residual[long$axis == a]
    data.frame(axis = a, mean = mean(r), sd = sd(r), n = length(r))
  }))
  if (any(!is.finite(scaling$sd)) || any(scaling$sd <= 0)) stop04("Residual scaling invalid")
  j <- match(long$axis, scaling$axis)
  long$residual_z <- (long$residual - scaling$mean[j]) / scaling$sd[j]
  wide <- folds$outer
  for (a in s$axes) {
    z <- long[long$axis == a, , drop = FALSE]
    z <- z[match(wide$patient_id, z$patient_id), , drop = FALSE]
    for (f in c("observed_delta", "expected_delta", "residual", "residual_z"))
      wide[[paste0(f, "_", a)]] <- z[[f]]
  }
  for (z in list(long, wide, scaling)) assert_fields04(names(z))
  write_tsv(long, path(stage04_long))
  write_tsv(wide, path(stage04_wide))
  write_tsv(scaling, path(stage04_scaling))
  write_tsv(do.call(rbind, tuning), path("results/stage04_residuals/nested_tuning_results.tsv"))
  write_tsv(do.call(rbind, coefficients), path("results/stage04_residuals/fold_coefficients.tsv"))
  benchmarks <- do.call(rbind, benchmark_rows)
  write_tsv(benchmarks, path("data/processed/residuals/neotrip_oof_benchmarks_v1.tsv"))
  invisible(list(long = long, wide = wide, scaling = scaling, benchmarks = benchmarks))
}
read_numeric04 <- function(file, character_fields) {
  x <- read_tsv(path(file))
  assert_fields04(names(x))
  for (f in setdiff(names(x), character_fields)) x[[f]] <- as.numeric(x[[f]])
  x
}
