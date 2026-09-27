source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/outcome_manifest.R"))
source(path("R/primary_inference.R"))
source(path("R/stage05_v2.R"))
source(path("R/stage05_v2_validation.R"))
run_stage("stage05_v2_tests", {
  out <- validate_results_v2()
  s <- read_spec_v2()
  old <- read_spec05()
  stopifnot(identical(s$directions, old$directions), identical(s$base_terms, old$base_terms))
  stopifnot(identical(s$primary_terms, c(old$primary_terms, "delta_stromal_fraction_10")))
  stopifnot(!"delta_normal_epithelial_10" %in% s$primary_terms)
  x <- read_numeric05(analysis_v2, c("patient_id", "arm"))
  m <- model.matrix(reformulate(s$primary_terms), frame_v2(x, s$axes[1], s))
  duplicate <- cbind(m, duplicate = m[, 2])
  stopifnot(!design_check_v2(duplicate, "synthetic", "duplicate")$passed)
  constant <- cbind(m, constant = 0)
  stopifnot(!design_check_v2(constant, "synthetic", "constant")$passed)
  bad <- m
  bad[1, 2] <- Inf
  stopifnot(!design_check_v2(bad, "synthetic", "nonfinite")$passed)
  toy <- data.frame(y = c(0, 0, 0, 1, 1, 1), x = c(-3, -2, -1, 1, 2, 3))
  detected <- glm(y ~ x, data = toy, family = binomial("logit"), method = detectseparation::detect_separation)
  stopifnot(any(is.infinite(coef(detected))))
  wide <- read_numeric05(wide_file05, c("patient_id", "arm"))
  stopifnot(identical(x$patient_id, wide$patient_id))
  for (a in s$axes) for (f in paste0(c("residual_", "residual_z_"), a))
    stopifnot(isTRUE(all.equal(x[[f]], wide[[f]], tolerance = 1e-12)))
  check <- read_tsv(path("results/stage05_primary/composition_constraint_check_v2.tsv"))
  sample <- read_tsv(path("results/stage05_primary/composition_sample_sum_check_v2.tsv"))
  stopifnot(nrow(check) == 150L, all(check$passed == "TRUE"),
    nrow(sample) == 401L, all(sample$passed == "TRUE"))
  state <- jsonlite::read_json(path("metadata/outcome/stage05_execution_v2.json"), simplifyVector = TRUE)
  amendment <- yaml::read_yaml(path(amend_v2))
  stopifnot(amendment$freeze_timestamp <= state$started, s$freeze_timestamp <= state$started)
  tests <- c("prior_stage_and_v1_bytes_preserved", "amendment_frozen_before_models",
    "directions_unchanged", "exact_primary_and_base_terms", "reference_component_excluded",
    "duplicate_predictor_rejected", "constant_predictor_rejected", "nonfinite_predictor_rejected",
    "complete_separation_detected", "raw_and_scaled_residuals_unchanged",
    "sample_and_pair_composition_verified", "saved_model_likelihood_verified",
    "inference_reproduced", "bonferroni_verified", "familywise_intervals_verified")
  write_tsv(data.frame(check = tests, passed = TRUE),
    path("results/stage05_primary/validation_tests_v2.tsv"))
})
cat("V2 checks passed\n")
