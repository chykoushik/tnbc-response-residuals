source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/outcome_manifest.R"))
source(path("R/primary_inference.R"))
source(path("R/stage05_v2.R"))
source(path("R/stage05_v2_validation.R"))
run_stage("17_finalize_stage05_v2", {
  status <- list(stage = 5L, version = "v2", timestamp = utc_now(), previous_decision = "MODIFY_REVIEW",
    methodology_amendment = amend_v2, amendment_reason = "four_component_composition",
    stage04_verified = FALSE, residual_freeze_verified = FALSE, four_component_constraint_verified = FALSE,
    patients = 0L, pcr = 0L, residual_disease = 0L, treatment_CT = 0L, treatment_CTA = 0L,
    primary_models = 0L, models_converged = 0L, rank_failures = NA_integer_,
    separation_failures = NA_integer_, missing_data_failures = NA_integer_, identity_failures = NA_integer_,
    direction_matching_pathways = NA_integer_, effect_size_gate_pathways = NA_integer_,
    multiplicity_gate_pathways = NA_integer_, both_effect_and_multiplicity_gate = NA_integer_,
    decision = "STOP_REVIEW", decision_reason = "Validation incomplete")
  tryCatch({
    verify_v2()
    status$stage04_verified <- TRUE
    status$residual_freeze_verified <- TRUE
    s <- read_spec_v2()
    state <- jsonlite::read_json(path("metadata/outcome/stage05_execution_v2.json"), simplifyVector = TRUE)
    if (!identical(state$specification_sha256, sha256(path(spec_v2))) ||
        !identical(state$amendment_sha256, sha256(path(amend_v2)))) stop05("V2 amendment changed after fitting")
    x <- read_numeric05(analysis_v2, c("patient_id", "arm"))
    status$patients <- nrow(x)
    status$pcr <- sum(x$pcr == 1L)
    status$residual_disease <- sum(x$pcr == 0L)
    status$treatment_CT <- sum(x$arm == "CT")
    status$treatment_CTA <- sum(x$arm == "CT/A")
    qc <- read_tsv(path("results/stage05_primary/outcome_qc_v2.tsv"))
    status$identity_failures <- as.integer(qc$observed[qc$check == "identity_failures"])
    status$missing_data_failures <- as.integer(qc$observed[qc$check == "missing_values"])
    if (any(qc$passed != "TRUE")) stop05("Outcome QC differs")
    paired <- read_tsv(path("results/stage05_primary/composition_constraint_check_v2.tsv"))
    sample <- read_tsv(path("results/stage05_primary/composition_sample_sum_check_v2.tsv"))
    status$four_component_constraint_verified <- nrow(paired) == 150L && all(paired$passed == "TRUE") &&
      nrow(sample) == 401L && all(sample$passed == "TRUE")
    if (!status$four_component_constraint_verified) stop05("Composition QC failed", "MODIFY_REVIEW")
    design <- read_tsv(path("results/stage05_primary/design_matrix_check_v2.tsv"))
    diagnostics <- read_tsv(path("results/stage05_primary/model_diagnostics_v2.tsv"))
    primary <- diagnostics[diagnostics$model == s$primary_model_name, ]
    status$primary_models <- nrow(primary)
    status$models_converged <- sum(primary$converged == "TRUE")
    status$rank_failures <- sum(as.integer(primary$rank) != as.integer(primary$number_parameters))
    status$separation_failures <- sum(primary$separation == "TRUE")
    if (state$decision != "PROCEED") stop05(state$reason, state$decision)
    if (status$primary_models != 6L || nrow(diagnostics) != 12L || any(diagnostics$passed != "TRUE") ||
        any(design$passed != "TRUE")) stop05("Model diagnostics failed", "MODIFY_REVIEW")
    out <- validate_results_v2()
    tests <- read_tsv(path("results/stage05_primary/validation_tests_v2.tsv"))
    if (any(tests$passed != "TRUE")) stop05("V2 regression checks failed")
    required <- c("figures/stage05_primary/primary_residual_forest_v2.pdf",
      "figures/stage05_primary/primary_residual_forest_v2.png",
      "figures/stage05_primary/composition_adjustment_comparison_v2.pdf",
      "figures/stage05_primary/composition_adjustment_comparison_v2.png",
      "tables/stage05_primary_residual_associations_v2.tsv")
    if (!all(file.exists(path(required)))) stop05("Required v2 outputs missing")
    status$direction_matching_pathways <- sum(out$primary$direction_matches)
    status$effect_size_gate_pathways <- sum(out$primary$effect_size_gate)
    status$multiplicity_gate_pathways <- sum(out$primary$multiplicity_gate)
    status$both_effect_and_multiplicity_gate <- sum(out$primary$effect_size_gate & out$primary$multiplicity_gate)
    old_lock <- renv::lockfile_read(path("metadata/modeling/stage04_renv.lock"))
    renv::snapshot(project = root, lockfile = path("metadata/outcome/stage05_environment_v2.lock"),
      packages = unique(c(names(old_lock$Packages), "detectseparation")), prompt = FALSE)
    verify_v2()
    status$decision <- "PROCEED"
    status$decision_reason <- paste("Four-component composition verified for 401 samples and all 150 pairs.",
      "Six prespecified primary models are full rank and converge without separation or non-finite estimates;",
      "12 model results reproduce. Stage 1-4 and Stage 5 v1 files remain unchanged.",
      "No pathway passes the Bonferroni multiplicity gate; PROCEED is based on technical validity.",
      "Stage 6 robustness and external-validation eligibility remain deferred.")
  }, error = function(e) {
    status$decision <<- if (inherits(e, "stage05_error")) e$decision else "STOP_REVIEW"
    status$decision_reason <<- conditionMessage(e)
  })
  packages <- unique(c("stats", "detectseparation",
    unlist(tools::package_dependencies("detectseparation", db = installed.packages(), recursive = TRUE), use.names = FALSE)))
  versions <- data.frame(R_version = R.version.string, platform = R.version$platform, package = packages,
    package_version = vapply(packages, function(p) as.character(packageVersion(p)), ""))
  write_tsv(versions, path("metadata/outcome/software_versions_v2.tsv"))
  current <- read_tsv(path("logs/software_versions.tsv"))
  write_tsv(rbind(current[!current$package %in% packages, ], versions), path("logs/software_versions.tsv"))
  jsonlite::write_json(status, path("results/stage05_primary/stage05_status_v2.json"),
    auto_unbox = TRUE, pretty = TRUE, na = "null")
  row <- data.frame(timestamp = status$timestamp, stage = 5L, decision = status$decision,
    reason = paste("v2:", status$decision_reason), inputs = spec_v2,
    outputs = "results/stage05_primary/stage05_status_v2.json")
  write.table(row, path("logs/decision_log.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE, append = TRUE)
})
files <- unlist(lapply(path(c(dirs05, "tables")), list.files, full.names = TRUE, recursive = TRUE), use.names = FALSE)
files <- files[grepl("_v2[.]", basename(files))]
files <- c(files, path(c("R/stage05_v2.R", "R/stage05_v2_validation.R",
  "scripts/15_amend_stage05_composition.R", "scripts/16_run_primary_outcome_v2.R",
  "scripts/17_finalize_stage05_v2.R", "scripts/test_stage05_v2.R", "scripts/install_stage05_v2.R",
  "logs/15_amend_stage05_composition.log", "logs/16_run_primary_outcome_v2.log",
  "logs/17_finalize_stage05_v2.log", "logs/stage05_v2_tests.log", "logs/stage05_v2_install.log")))
files <- sort(unique(files[file.exists(files) & basename(files) != "stage05_output_checksums_v2.tsv"]))
write_tsv(data.frame(local_relative_path = substring(files, nchar(root) + 2L),
  file_size_bytes = file.info(files)$size, sha256 = vapply(files, sha256, "")),
  path("results/stage05_primary/stage05_output_checksums_v2.tsv"))
cat("results/stage05_primary/stage05_status_v2.json\n")
