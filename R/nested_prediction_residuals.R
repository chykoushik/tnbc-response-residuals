split_path07 <- function(r, k) sprintf("repeat_%02d/fold_%02d", r, k)
numeric_expected07 <- function(x) {
  m <- x[, expected07, drop = FALSE]
  m$arm <- as.integer(m$arm == "CT/A")
  m <- as.matrix(m)
  assert07(!any(grepl("pcr|response|outcome|^delta_", colnames(m), ignore.case = TRUE)) &&
    all(is.finite(m)), "Expected-response design invalid")
  m
}
residual_split07 <- function(x, outer, crossfit, tuning, r, k) {
  folder <- path("models/stage07/expected_response", split_path07(r, k))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  checkpoint <- file.path(folder, "residual_results.rds")
  source_hash <- digest::digest(list(x, outer, crossfit, tuning, sha256(path(spec07_file))), algo = "sha256")
  if (file.exists(checkpoint)) {
    cached <- readRDS(checkpoint)
    assert07(identical(cached$source_hash, source_hash), "Residual checkpoint changed")
    return(cached)
  }
  train <- which(outer$outer_fold != k); test <- which(outer$outer_fold == k)
  crossfit <- crossfit[match(x$patient_id[train], crossfit$patient_id), ]
  assert07(identical(crossfit$patient_id, x$patient_id[train]) &&
    !any(crossfit$patient_id %in% x$patient_id[test]), "Cross-fit cohort contaminated")
  m <- numeric_expected07(x)
  rows <- audits <- scales <- fit_audits <- list()
  for (a in axes07) {
    target <- x[[paste0("delta_", a)]]
    expected <- rep(NA_real_, nrow(x)); seen <- integer(nrow(x)); fits <- list()
    for (cf in 0:5) {
      fit_ids <- if (cf == 0) train else train[crossfit$residual_fold != cf]
      predict_ids <- if (cf == 0) test else train[crossfit$residual_fold == cf]
      inner <- tuning[tuning$residual_fold == cf, ]
      inner <- inner[match(x$patient_id[fit_ids], inner$patient_id), ]
      overlap <- length(intersect(x$patient_id[test], x$patient_id[fit_ids])) +
        length(intersect(x$patient_id[test], inner$patient_id)) +
        length(intersect(fit_ids, predict_ids))
      assert07(overlap == 0 && !anyNA(inner$patient_id) &&
        identical(inner$patient_id, x$patient_id[fit_ids]) &&
        setequal(inner$tuning_fold, 1:5), "Expected-response patient leakage")
      warnings <- character()
      fit <- withCallingHandlers(glmnet::cv.glmnet(m[fit_ids, , drop = FALSE], target[fit_ids],
        family = "gaussian", alpha = 0, intercept = TRUE, standardize = TRUE,
        type.measure = "mse", nlambda = 100, lambda.min.ratio = 1e-4,
        foldid = inner$tuning_fold, nfolds = 5, grouped = TRUE, parallel = FALSE, keep = TRUE),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
      assert07(!length(warnings) && fit$glmnet.fit$jerr == 0,
        "Expected-response fit failed", "MODIFY_REVIEW")
      predicted <- as.numeric(predict(fit, newx = m[predict_ids, , drop = FALSE], s = "lambda.1se"))
      expected[predict_ids] <- predicted; seen[predict_ids] <- seen[predict_ids] + 1L
      fits[[as.character(cf)]] <- list(fit = fit, residual_fold = cf, training_ids = x$patient_id[fit_ids],
        prediction_ids = x$patient_id[predict_ids], tuning = inner, predictors = colnames(m),
        training_sha256 = digest::digest(list(m[fit_ids, , drop = FALSE], target[fit_ids]), algo = "sha256"),
        clinical_labels_used = FALSE)
      fit_row <- data.frame(repetition = r, outer_fold = k, axis = a, residual_fold = cf,
        n_fit = length(fit_ids), n_predict = length(predict_ids), selected_lambda = fit$lambda.1se,
        lambda_min = fit$lambda.min, patient_overlap = overlap, clinical_labels_used = FALSE,
        training_only_tuning = TRUE, warnings = length(warnings))
      names(fit_row)[1] <- "repeat"
      fit_audits[[paste(a, cf)]] <- fit_row
    }
    residual <- target - expected
    assert07(all(seen == 1L) && all(is.finite(residual)) && sd(residual[train]) > 0, "Prediction residual coverage failed")
    mean_train <- mean(residual[train]); sd_train <- sd(residual[train])
    z <- data.frame(repetition = r, outer_fold = k, patient_id = x$patient_id,
      set = ifelse(seq_len(nrow(x)) %in% train, "train", "test"), axis = a,
      observed_delta = target, expected_delta = expected, residual = residual,
      residual_z = (residual - mean_train) / sd_train)
    names(z)[1] <- "repeat"
    rows[[a]] <- z
    audit <- data.frame(repetition = r, outer_fold = k, axis = a, n_train = length(train), n_test = length(test),
      training_crossfit_complete = all(seen[train] == 1L), test_prediction_training_only = TRUE,
      patient_overlap = 0L, missing_residuals = sum(is.na(residual)), nonfinite_residuals = sum(!is.finite(residual)))
    names(audit)[1] <- "repeat"
    audits[[a]] <- audit
    scaling <- data.frame(repetition = r, outer_fold = k, axis = a,
      mean = mean_train, sd = sd_train, n_train = length(train), training_ids_sha256 = digest::digest(x$patient_id[train], algo = "sha256"))
    names(scaling)[1] <- "repeat"
    scales[[a]] <- scaling
    saveRDS(fits, file.path(folder, paste0(a, ".rds")))
  }
  result <- list(source_hash = source_hash, residuals = do.call(rbind, rows),
    audit = do.call(rbind, audits), scaling = do.call(rbind, scales), fits = do.call(rbind, fit_audits))
  saveRDS(result, checkpoint)
  result
}
verify_prediction_residuals07 <- function() {
  file <- path("metadata/prediction/stage07_residual_freeze_v1.yml")
  assert07(file.exists(file) && identical(readLines(paste0(file, ".sha256")), sha256(file)), "Prediction residual freeze changed")
  freeze <- yaml::read_yaml(file)
  for (i in seq_along(freeze$files$local_relative_path)) {
    assert07(sha256(path(freeze$files$local_relative_path[i])) == freeze$files$sha256[i],
      "Prediction residual bytes changed")
  }
  assert07(identical(freeze$clinical_labels_used, FALSE), "Clinical information used in residual construction")
}
