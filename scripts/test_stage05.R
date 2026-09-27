source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/outcome_manifest.R"))
source(path("R/primary_inference.R"))
run_stage("stage05_tests", {
  verify_stage04_05()
  spec <- read_spec05()
  state <- jsonlite::read_json(path("metadata/outcome/stage05_execution_v1.json"), simplifyVector = TRUE)
  stopifnot(state$specification_verified, identical(state$specification_sha256, sha256(path(spec_file05))))
  stopifnot(state$outcome_access_started >= spec$freeze_timestamp)
  x <- read_numeric05(analysis_file05, c("patient_id", "arm"))
  input <- read_numeric05("data/processed/pathways/neotrip_model_inputs_blinded_v1.tsv", c("patient_id", "arm"))
  wide <- read_numeric05(wide_file05, c("patient_id", "arm"))
  stopifnot(identical(x$patient_id, wide$patient_id), identical(x$arm, wide$arm))
  i <- match(x$patient_id, input$patient_id)
  stopifnot(nrow(x) == 150L, sum(x$pcr == 1L) == 63L, sum(x$pcr == 0L) == 87L)
  for (a in spec$axes) {
    stopifnot(isTRUE(all.equal(x[[paste0("residual_", a)]], wide[[paste0("residual_", a)]], tolerance = 1e-12)))
    stopifnot(isTRUE(all.equal(x[[paste0("residual_z_", a)]], wide[[paste0("residual_z_", a)]], tolerance = 1e-12)))
    z <- x[[paste0("baseline_z_", a)]]
    stopifnot(abs(mean(z)) < 1e-10, abs(sd(z) - 1) < 1e-10)
  }
  for (f in c("delta_cancer_epithelial", "delta_immune_fraction", "delta_stromal_fraction")) {
    stopifnot(isTRUE(all.equal(x[[f]], input[[f]][i] / 100, tolerance = 1e-12)))
    stopifnot(isTRUE(all.equal(x[[paste0(f, "_10")]], input[[f]][i] / 10, tolerance = 1e-12)))
  }
  check <- read_numeric05("results/stage05_primary/composition_constraint_check.tsv", c("patient_id", "passed"))
  total <- rowSums(x[, c("delta_cancer_epithelial", "delta_immune_fraction", "delta_stromal_fraction")])
  stopifnot(isTRUE(all.equal(check$delta_sum, total, tolerance = 1e-12)))
  if (any(abs(total) > spec$composition_constraint_tolerance)) {
    stopifnot(identical(state$decision, "MODIFY_REVIEW"), identical(state$models_attempted, FALSE))
    stopifnot(!file.exists(path("results/stage05_primary/primary_models_v1.rds")))
    stopifnot(!file.exists(path("results/stage05_primary/primary_residual_associations.tsv")))
  }
  checks <- c("stage04_freeze_intact", "specification_precedes_outcome_access", "exact_patient_join",
    "outcome_counts_match", "raw_residuals_unchanged", "standardized_residuals_unchanged",
    "baseline_scaling_reproduced", "explicit_composition_unit_conversion",
    "ten_percentage_point_scaling", "composition_constraint_reproduced", "failed_constraint_prevents_models")
  write_tsv(data.frame(check = checks, passed = TRUE), path("results/stage05_primary/stage05_validation.tsv"))
})
cat("Checks passed\n")
