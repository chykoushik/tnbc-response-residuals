source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/prediction_folds.R"))
run07("24_freeze_stage07_spec", {
  gate07()
  s <- list(version = "v1", freeze_timestamp = utc_now(), master_seed = 319650L,
    outer_folds = 5L, outer_repeats = 10L, prediction_inner_folds = 5L,
    expected_response_crossfit_folds = 5L, expected_response_tuning_folds = 5L,
    classifier = list(method = "ridge_logistic", alpha = 0, lambda_rule = "lambda.1se",
      standardize = TRUE, intercept = TRUE, nlambda = 100L, lambda_min_ratio = 1e-4),
    expected_response = list(method = "gaussian_ridge", alpha = 0, lambda_rule = "lambda.1se",
      standardize = TRUE, intercept = TRUE, nlambda = 100L, lambda_min_ratio = 1e-4),
    primary_metrics = c("roc_auc", "brier_score", "log_loss", "calibration_intercept", "calibration_slope"),
    decision_curve_analysis = "SKIP",
    decision_curve_reason = "no external transportability and clinical utility remains hypothesis-generating",
    axes = axes07, features = features07, expected_response_predictors = expected07,
    outer_seed = "319650 + repeat_number", outer_stratification = "treatment_arm_only",
    split_seed = "319650 + 10000 * repeat + 100 * outer_fold",
    residual_crossfit_seed = "split_seed + 1",
    expected_tuning_seed = "split_seed + 10 + residual_crossfit_fold",
    full_training_expected_tuning_seed = "split_seed + 20",
    classifier_seed = "split_seed + 30", classifier_stratification = "pcr",
    RNGkind = c("Mersenne-Twister", "Inversion", "Rejection"),
    residual_feature_values = "residual_z using outer-training cross-fitted residual mean and sample SD",
    classifier_standardization = "glmnet training-only standardization in each fit",
    composition_units = "fraction; Stage 3 percentage values divided by 100 for classifier inputs",
    normal_epithelial_reference = TRUE, calibration_clip = 1e-6,
    auc_method = "Mann-Whitney rank statistic with average ranks for ties",
    summary_range = "2.5th and 97.5th empirical repeat percentiles, type 7; not formal confidence intervals",
    calibration_display = "Patient-level mean OOF probabilities; five equal-count bins per model for display only",
    reproduction_splits = c("repeat_01/fold_01", "repeat_05/fold_03", "repeat_10/fold_05"),
    glmnet_version = as.character(packageVersion("glmnet")))
  freeze07(spec07_file, s, yaml = TRUE)
  read_spec07()
})
cat("Prediction specification frozen\n")
