design07 <- function(x, residuals, r, k, model) {
  z <- residuals[residuals[["repeat"]] == r & residuals$outer_fold == k, ]
  for (a in axes07) {
    v <- z[z$axis == a, ]
    assert07(nrow(v) == 150 && !anyDuplicated(v$patient_id) &&
      setequal(v$patient_id, x$patient_id), "Residual feature identity failed")
    x[[paste0("residual_", a)]] <- v$residual_z[match(x$patient_id, v$patient_id)]
  }
  m <- x[, features07[[model]], drop = FALSE]
  m$arm <- as.integer(m$arm == "CT/A")
  m <- as.matrix(m)
  assert07(all(is.finite(m)), "Classifier design nonfinite")
  m
}
classifiers07 <- function(x, residuals, outer, inner) {
  scores <- tuning <- list()
  for (r in 1:10) for (k in 1:5) {
    o <- outer[outer[["repeat"]] == r, ]
    train <- which(o$outer_fold != k); test <- which(o$outer_fold == k)
    f <- inner[inner[["repeat"]] == r & inner$outer_fold == k, ]
    f <- f[match(x$patient_id[train], f$patient_id), ]
    assert07(identical(f$patient_id, x$patient_id[train]) &&
      !any(f$patient_id %in% x$patient_id[test]), "Classifier inner fold leakage")
    folder <- path("models/stage07", split_path07(r, k))
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    for (model in names(features07)) {
      m <- design07(x, residuals, r, k, model)
      warnings <- character()
      fit <- withCallingHandlers(glmnet::cv.glmnet(m[train, , drop = FALSE], x$pcr[train],
        family = "binomial", alpha = 0, standardize = TRUE, intercept = TRUE,
        type.measure = "deviance", nlambda = 100, lambda.min.ratio = 1e-4,
        foldid = f$inner_fold, nfolds = 5, grouped = TRUE, parallel = FALSE, keep = TRUE),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
      assert07(!length(warnings) && fit$glmnet.fit$jerr == 0, "Classifier fit failed", "MODIFY_REVIEW")
      probability <- as.numeric(predict(fit, newx = m[test, , drop = FALSE], s = "lambda.1se", type = "response"))
      assert07(all(is.finite(probability) & probability > 0 & probability < 1),
        "Classifier probability invalid", "MODIFY_REVIEW")
      saved <- list(fit = fit, feature_names = colnames(m), training_ids = x$patient_id[train],
        test_ids = x$patient_id[test], inner = f, standardization_training_only = TRUE,
        training_sha256 = digest::digest(list(m[train, , drop = FALSE], x$pcr[train]), algo = "sha256"))
      saveRDS(saved, file.path(folder, paste0(model, ".rds")))
      row <- data.frame(repetition = r, outer_fold = k, patient_id = x$patient_id[test],
        pcr = x$pcr[test], model = model, predicted_probability = probability)
      names(row)[1] <- "repeat"
      scores[[paste(r, k, model)]] <- row
      row <- data.frame(repetition = r, outer_fold = k, model = model, n_train = length(train),
        n_test = length(test), lambda_min = fit$lambda.min, lambda_1se = fit$lambda.1se,
        selected_lambda = fit$lambda.1se)
      names(row)[1] <- "repeat"
      tuning[[paste(r, k, model)]] <- row
    }
    cat("Classifier split", r, k, "\n")
  }
  list(scores = do.call(rbind, scores), tuning = do.call(rbind, tuning))
}
audit_models07 <- function(x, residuals, outer, scores) {
  assert07(requireNamespace("glmnet", quietly = TRUE), "glmnet unavailable", "MODIFY_REVIEW")
  blind <- blinded07(); em <- numeric_expected07(blind)
  audits <- reproductions <- list()
  for (r in 1:10) for (k in 1:5) {
    o <- outer[outer[["repeat"]] == r, ]
    train_ids <- x$patient_id[o$outer_fold != k]; test_ids <- x$patient_id[o$outer_fold == k]
    checks <- c(test_absent_expected_training = TRUE, test_absent_expected_tuning = TRUE,
      test_absent_classifier_training = TRUE, test_pcr_unused_expected = TRUE,
      test_pcr_unused_residuals = TRUE, training_residual_without_self = TRUE,
      scaling_outer_training_only = TRUE, classifier_preprocessing_outer_training_only = TRUE)
    sample_split <- (r == 1 && k == 1) || (r == 5 && k == 3) || (r == 10 && k == 5)
    for (a in axes07) {
      fits <- readRDS(path("models/stage07/expected_response", split_path07(r, k), paste0(a, ".rds")))
      z <- residuals[residuals[["repeat"]] == r & residuals$outer_fold == k & residuals$axis == a, ]
      tr <- z$patient_id %in% train_ids
      checks["scaling_outer_training_only"] <- checks["scaling_outer_training_only"] &&
        all(z$set == ifelse(tr, "train", "test")) &&
        max(abs(z$residual_z - (z$residual - mean(z$residual[tr])) / sd(z$residual[tr]))) < 1e-10
      for (f in fits) {
        checks["test_absent_expected_training"] <- checks["test_absent_expected_training"] &&
          !any(test_ids %in% f$training_ids) && all(f$training_ids %in% train_ids)
        checks["test_absent_expected_tuning"] <- checks["test_absent_expected_tuning"] &&
          identical(f$training_ids, f$tuning$patient_id) && !any(test_ids %in% f$tuning$patient_id)
        checks["test_pcr_unused_expected"] <- checks["test_pcr_unused_expected"] &&
          identical(f$clinical_labels_used, FALSE) && identical(f$predictors, expected07)
        checks["training_residual_without_self"] <- checks["training_residual_without_self"] &&
          !any(f$prediction_ids %in% f$training_ids)
        idx <- match(f$training_ids, blind$patient_id)
        checks["test_pcr_unused_residuals"] <- checks["test_pcr_unused_residuals"] &&
          identical(f$training_sha256, digest::digest(list(em[idx, , drop = FALSE],
            blind[[paste0("delta_", a)]][idx]), algo = "sha256"))
        if (sample_split) {
          idx <- match(f$prediction_ids, blind$patient_id)
          pred <- as.numeric(predict(f$fit, newx = em[idx, , drop = FALSE], s = "lambda.1se"))
          error <- max(abs(pred - z$expected_delta[match(f$prediction_ids, z$patient_id)]))
          row <- data.frame(repetition = r, outer_fold = k, type = "expected_response",
            model = paste(a, f$residual_fold, sep = "/"), max_absolute_difference = error, passed = error <= 1e-10)
          names(row)[1] <- "repeat"; reproductions[[length(reproductions) + 1L]] <- row
        }
      }
    }
    for (model in names(features07)) {
      f <- readRDS(path("models/stage07", split_path07(r, k), paste0(model, ".rds")))
      m <- design07(x, residuals, r, k, model)
      tr <- match(train_ids, x$patient_id)
      checks["test_absent_classifier_training"] <- checks["test_absent_classifier_training"] &&
        identical(f$training_ids, train_ids) && identical(f$inner$patient_id, train_ids) &&
        !any(test_ids %in% f$training_ids)
      checks["classifier_preprocessing_outer_training_only"] <-
        checks["classifier_preprocessing_outer_training_only"] &&
        identical(f$standardization_training_only, TRUE) && identical(f$feature_names, features07[[model]]) &&
        identical(f$training_sha256, digest::digest(list(m[tr, , drop = FALSE], x$pcr[tr]), algo = "sha256"))
      if (sample_split) {
        idx <- match(test_ids, x$patient_id)
        pred <- as.numeric(predict(f$fit, newx = m[idx, , drop = FALSE], s = "lambda.1se", type = "response"))
        s <- scores[scores[["repeat"]] == r & scores$outer_fold == k & scores$model == model, ]
        error <- max(abs(pred - s$predicted_probability[match(test_ids, s$patient_id)]))
        row <- data.frame(repetition = r, outer_fold = k, type = "classifier", model = model,
          max_absolute_difference = error, passed = error <= 1e-10)
        names(row)[1] <- "repeat"; reproductions[[length(reproductions) + 1L]] <- row
      }
    }
    row <- data.frame(repetition = r, outer_fold = k, as.list(checks), violations = sum(!checks))
    names(row)[1] <- "repeat"; audits[[length(audits) + 1L]] <- row
  }
  audit <- do.call(rbind, audits); reproduction <- do.call(rbind, reproductions)
  write07(audit, "results/stage07_prediction/leakage_audit.tsv")
  write07(reproduction, "results/stage07_prediction/model_reproduction_check.tsv")
  assert07(sum(audit$violations) == 0, "Prediction leakage audit failed")
  assert07(all(reproduction$passed), "Saved model reproduction failed", "MODIFY_REVIEW")
}
