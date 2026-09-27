source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/fold_generation.R"))
source(path("R/expected_response.R"))
source(path("R/residual_qc.R"))
run04("stage04_tests", {
  rejected <- function(code) tryCatch({force(code); FALSE}, stage04_error = function(e) TRUE)
  verify_stage03()
  x <- read_input04()
  s <- read_spec04()
  folds <- freeze_folds04(x, s)
  stopifnot(rejected(assert_fields04(c("patient_id", "Early.pCR"))))
  stopifnot(rejected(assert_fields04(c("patient_id", "response"))))
  stopifnot(rejected(assert_fields04(c("patient_id", "residual_disease"))))
  stopifnot(length(forbidden04(c("arm", "observed_delta", "residual"))) == 0L)
  stopifnot(identical(colnames(design04(x, s)), s$predictors))
  stopifnot(!any(grepl("^delta_|d1c2", colnames(design04(x, s)))))
  f <- assign_folds04(x$patient_id, x$arm, 10L, s$master_seed, s$rng_kind)
  reversed <- assign_folds04(rev(x$patient_id), rev(x$arm), 10L, s$master_seed, s$rng_kind)
  stopifnot(identical(f, rev(reversed)), identical(f, folds$outer$outer_fold))
  for (k in 1:10) {
    inner <- folds$inner[folds$inner$outer_fold == k, ]
    expected <- assign_folds04(inner$patient_id, inner$arm, 10L, s$master_seed + k, s$rng_kind)
    stopifnot(identical(expected, inner$inner_fold))
  }
  bad <- folds$inner
  bad$patient_id[1] <- folds$outer$patient_id[folds$outer$outer_fold == bad$outer_fold[1]][1]
  stopifnot(rejected(validate_folds04(x, s, folds$outer, bad)))
  bad_arm <- x
  bad_arm$arm[1] <- "unknown"
  stopifnot(rejected(design04(bad_arm, s)))
  audited <- audit_saved04(x, s, folds)
  performance <- performance04(audited$benchmarks, s$axes)
  constant <- audited$long
  constant$residual[constant$axis == s$axes[1]] <- 0
  qc <- suppressWarnings(residual_qc04(constant, performance, s$axes))
  stopifnot(!qc$passed[1])
  for (a in s$axes) {
    z <- audited$long$residual_z[audited$long$axis == a]
    stopifnot(abs(mean(z)) < 1e-10, abs(sd(z) - 1) < 1e-10)
  }
  checks <- c("forbidden_field_patterns_rejected", "allowed_molecular_fields",
    "exact_baseline_predictors", "postbaseline_predictors_excluded",
    "fold_seed_reproduction", "fold_order_invariance", "inner_seed_reproduction",
    "heldout_patient_contamination_rejected", "unexpected_arm_rejected",
    "all_saved_model_predictions_reproduced", "training_only_benchmarks_reproduced",
    "constant_residual_rejected", "standardized_residuals_valid", "prior_stage_bytes_preserved")
  write_tsv(data.frame(check = checks, passed = TRUE), path("results/stage04_residuals/stage04_regression_tests.tsv"))
})
cat("Checks passed\n")
