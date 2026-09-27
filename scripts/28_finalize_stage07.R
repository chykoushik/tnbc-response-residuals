source("R/paths.R")
source("R/checksums.R")
source("R/prediction_folds.R")
source("R/nested_prediction_residuals.R")
source("R/prediction_models.R")
run07("28_finalize_stage07", {
  status <- list(stage = 7L, timestamp = utc_now(), stage06_verified = FALSE,
    patients = 0L, outer_folds = 5L, outer_repeats = 10L, outer_splits = 0L,
    prediction_models = 5L, oof_predictions = 0L, residual_pipeline_failures = 0L,
    classifier_failures = 0L, leakage_failures = 0L, reproduction_failures = 0L,
    all_metrics_complete = FALSE, P3_vs_P1_auc_mean_difference = NA_real_,
    P3_vs_P1_brier_mean_difference = NA_real_, P3_vs_P1_logloss_mean_difference = NA_real_,
    P4_vs_P2_auc_mean_difference = NA_real_, P4_vs_P2_brier_mean_difference = NA_real_,
    P4_vs_P2_logloss_mean_difference = NA_real_, external_transportability_status = "SKIP",
    decision_curve_status = "NOT_RUN", decision = "MODIFY_REVIEW", decision_reason = "Stage 7 incomplete")
  tryCatch({
    status$stage06_verified <- gate07()
    read_spec07(); verify_prediction_residuals07()
    failures <- list.files(path("results/stage07_prediction"), pattern = "_failure[.]json$", full.names = TRUE)
    if (length(failures)) {
      f <- jsonlite::read_json(failures[1], simplifyVector = TRUE)
      fail07(f$reason, f$decision)
    }
    x <- population07(); status$patients <- nrow(x)
    outer <- read_outer07()
    for (r in 1:10) {
      z <- outer[outer[["repeat"]] == r, ]
      assert07(identical(z$outer_fold, stratify07(x$patient_id, x$arm, 5L, 319650L + r)),
        "Outer folds do not reproduce arm-only seeds")
    }
    residuals <- read07("data/processed/prediction/stage07_nested_residuals_v1.tsv.gz")
    scores <- read07("data/processed/prediction/stage07_oof_predictions_v1.tsv.gz")
    audit_models07(x, residuals, outer, scores)
    status$oof_predictions <- nrow(scores)
    assert07(nrow(scores) == 7500 && !anyDuplicated(scores[, c("repeat", "model", "patient_id")]) &&
      all(table(scores[["repeat"]], scores$model) == 150) &&
      setequal(scores$model, names(features07)) &&
      all(is.finite(scores$predicted_probability) & scores$predicted_probability > 0 &
        scores$predicted_probability < 1), "OOF prediction coverage invalid", "MODIFY_REVIEW")
    for (r in 1:10) {
      z <- scores[scores[["repeat"]] == r, ]; o <- outer[outer[["repeat"]] == r, ]
      assert07(all(z$outer_fold == o$outer_fold[match(z$patient_id, o$patient_id)]) &&
        all(z$pcr == x$pcr[match(z$patient_id, x$patient_id)]), "OOF fold or outcome identity failed")
    }
    ra <- read07("results/stage07_prediction/prediction_residual_audit.tsv")
    status$residual_pipeline_failures <- sum(!ra$training_crossfit_complete |
      !ra$test_prediction_training_only | ra$patient_overlap != 0 |
      ra$missing_residuals != 0 | ra$nonfinite_residuals != 0)
    assert07(nrow(ra) == 300 && nrow(residuals) == 45000 &&
      !anyDuplicated(residuals[, c("repeat", "outer_fold", "axis", "patient_id")]) &&
      max(abs(residuals$residual - (residuals$observed_delta - residuals$expected_delta))) < 1e-10,
      "Residual coverage invalid", "MODIFY_REVIEW")
    tuning <- read07("results/stage07_prediction/classifier_tuning.tsv")
    status$outer_splits <- nrow(unique(tuning[, c("repeat", "outer_fold")]))
    status$classifier_failures <- sum(!is.finite(tuning$selected_lambda) |
      tuning$selected_lambda != tuning$lambda_1se)
    assert07(nrow(tuning) == 250 && status$outer_splits == 50 &&
      !anyDuplicated(tuning[, c("repeat", "outer_fold", "model")]), "Classifier fits incomplete", "MODIFY_REVIEW")
    la <- read07("results/stage07_prediction/leakage_audit.tsv")
    rp <- read07("results/stage07_prediction/model_reproduction_check.tsv")
    status$leakage_failures <- sum(la$violations)
    status$reproduction_failures <- sum(!rp$passed | rp$max_absolute_difference > 1e-10)
    assert07(nrow(la) == 50 && nrow(rp) == 123, "Audit coverage incomplete", "MODIFY_REVIEW")
    metrics <- read07("results/stage07_prediction/repeat_metrics.tsv")
    summary <- read07("results/stage07_prediction/prediction_summary.tsv")
    paired <- read07("results/stage07_prediction/paired_metric_summary.tsv")
    differences <- read07("results/stage07_prediction/paired_metric_differences.tsv")
    status$all_metrics_complete <- nrow(metrics) == 50 && nrow(summary) == 25 &&
      nrow(paired) == 15 && nrow(differences) == 50 &&
      all(is.finite(as.matrix(metrics[, -(1:2)]))) &&
      all(summary$range_label == "empirical repeat range")
    for (pair in c("P3_vs_P1", "P4_vs_P2")) for (metric in c("auc", "brier", "logloss")) {
      row <- paired$comparison == gsub("_vs_", " - ", pair) & paired$metric == paste0(metric, "_difference")
      assert07(sum(row) == 1, "Comparison result missing", "MODIFY_REVIEW")
      status[[paste0(pair, "_", metric, "_mean_difference")]] <- paired$mean_difference[row]
    }
    assert07(all(c(status$residual_pipeline_failures, status$classifier_failures,
      status$leakage_failures, status$reproduction_failures) == 0) && status$all_metrics_complete,
      "Required technical checks failed", "MODIFY_REVIEW")
    assert07(length(list.files(path("figures/stage07_prediction"), pattern = "[.](pdf|png)$")) == 6 &&
      read07("results/stage07_prediction/decision_curve_status.tsv")$status == "NOT_RUN",
      "Figure or decision curve status incomplete", "MODIFY_REVIEW")
    previous <- renv::lockfile_read(path("metadata/sensitivity/stage06_environment_v1.lock"))
    renv::snapshot(project = root, lockfile = path("metadata/prediction/stage07_environment_v1.lock"),
      packages = names(previous$Packages), prompt = FALSE)
    status$decision <- "PROCEED"
    status$decision_reason <- paste("All 50 outer splits, 250 classifiers and 7500 OOF predictions completed;",
      "residual, classifier, leakage and reproduction checks passed and all metrics are complete.",
      "Prediction is secondary; PROCEED reflects technical validity and does not require predictive improvement.",
      "External transportability remains SKIP; decision curves were not run.")
    verify_prior07()
  }, error = function(e) {
    status$decision <<- if (inherits(e, "stage07_error")) e$decision else "MODIFY_REVIEW"
    status$decision_reason <<- conditionMessage(e)
  })
  packages <- c("glmnet", "Matrix", "stats", "renv", "digest", "yaml", "jsonlite")
  versions <- data.frame(R_version = R.version.string, platform = R.version$platform,
    package = packages, package_version = vapply(packages, function(p) as.character(packageVersion(p)), ""))
  write07(versions, "metadata/prediction/software_versions_v1.tsv")
  old <- read_tsv(path("logs/software_versions.tsv"))
  write_tsv(rbind(old[!old$package %in% packages, ], versions), path("logs/software_versions.tsv"))
  jsonlite::write_json(status, path("results/stage07_prediction/stage07_status.json"),
    auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 16)
  write.table(data.frame(timestamp = status$timestamp, stage = 7L, decision = status$decision,
    reason = status$decision_reason, inputs = spec07_file, outputs = "results/stage07_prediction/stage07_status.json"),
    path("logs/decision_log.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
})
files <- unlist(lapply(path(dirs07), list.files, recursive = TRUE, full.names = TRUE), use.names = FALSE)
code <- c("R/prediction_folds.R", "R/nested_prediction_residuals.R", "R/prediction_models.R",
  "R/prediction_metrics.R", "R/prediction_calibration.R",
  "scripts/24_freeze_stage07_spec.R", "scripts/25_generate_stage07_folds.R",
  "scripts/26_generate_nested_prediction_residuals.R", "scripts/27_run_nested_prediction.R",
  "scripts/28_finalize_stage07.R")
logs <- list.files(path("logs"), pattern = "^(24_|25_|26_|27_|28_).*\\.log$", full.names = TRUE)
files <- sort(unique(c(files, path(code), logs)))
files <- files[file.exists(files) & basename(files) != "stage07_output_checksums.tsv"]
write_tsv(data.frame(local_relative_path = substring(files, nchar(root) + 2L),
  file_size_bytes = file.info(files)$size, sha256 = vapply(files, sha256, "")),
  path("results/stage07_prediction/stage07_output_checksums.tsv"))
cat("results/stage07_prediction/stage07_status.json\n")
