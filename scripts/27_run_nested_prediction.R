source("R/paths.R")
source("R/checksums.R")
source("R/prediction_folds.R")
source("R/nested_prediction_residuals.R")
source("R/prediction_models.R")
source("R/prediction_calibration.R")
source("R/prediction_metrics.R")
run07("27_run_nested_prediction", {
  gate07(); read_spec07(); verify_prediction_residuals07()
  x <- population07(); outer <- read_outer07()
  file <- "metadata/prediction/stage07_classifier_inner_folds_v1.tsv"
  freeze07(file, NULL); inner <- read07(file)
  residuals <- read07("data/processed/prediction/stage07_nested_residuals_v1.tsv.gz")
  result <- classifiers07(x, residuals, outer, inner)
  write07(result$scores, "data/processed/prediction/stage07_oof_predictions_v1.tsv.gz")
  write07(result$tuning, "results/stage07_prediction/classifier_tuning.tsv")
  scores <- read07("data/processed/prediction/stage07_oof_predictions_v1.tsv.gz")
  audit_models07(x, residuals, outer, scores)
  metrics07(scores, x)
  write07(data.frame(analysis = "decision_curve", status = "NOT_RUN",
    reason = "No external transportability; clinical utility remains hypothesis-generating"),
    "results/stage07_prediction/decision_curve_status.tsv")
})
cat("Nested prediction complete\n")
