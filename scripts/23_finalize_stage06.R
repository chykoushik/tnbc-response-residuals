source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/stage06_common.R"))
source(path("R/stage06_validation.R"))
run_stage("23_finalize_stage06", {
  status <- list(stage = 6L, timestamp = utc_now(), stage05_verified = FALSE,
    five_sensitivities_completed = FALSE, permutations_completed = FALSE,
    permutations_per_pathway = 0L, influence_completed = FALSE, random_gene_controls_completed = FALSE,
    arm_interaction_completed = FALSE, arm_specific_residuals_completed = FALSE,
    complete_case_verified = FALSE, technical_failures = 0L, full_eligible_pathways = 0L,
    hypothesis_generating_pathways = 0L, not_robust_pathways = 0L,
    external_transportability_decision = "NOT_EVALUATED", decision = "MODIFY_REVIEW",
    decision_reason = "Stage 6 not complete")
  tryCatch({
    status$stage05_verified <- stage05_gate06()
    failures <- list.files(path("results/stage06_robustness"), pattern = "_failure[.]json$", full.names = TRUE)
    if (length(failures)) {
      f <- jsonlite::read_json(failures[1], simplifyVector = TRUE)
      fail06(f$reason, f$decision)
    }
    validate_outputs06()
    crossfit <- read06("results/stage06_robustness/crossfit_reproduction_validation.tsv")
    assert06(nrow(crossfit) == 630 && all(crossfit$passed), "Cross-fit verification incomplete")
    status$five_sensitivities_completed <- nrow(read06("results/stage06_robustness/five_primary_sensitivities.tsv")) == 30
    permutation <- read06("results/stage06_robustness/permutation_summary.tsv")
    status$permutations_completed <- nrow(permutation) == 6 && all(permutation$permutations == 10000)
    status$permutations_per_pathway <- min(permutation$permutations)
    status$influence_completed <- nrow(read06("results/stage06_robustness/influence_analysis.tsv")) == 900
    status$random_gene_controls_completed <- nrow(read06("results/stage06_robustness/random_gene_null_results.tsv.gz")) == 600
    status$arm_interaction_completed <- nrow(read06("results/stage06_robustness/outcome_arm_interaction.tsv")) == 6
    status$arm_specific_residuals_completed <- nrow(read06("results/stage06_robustness/arm_specific_residual_sensitivity.tsv")) == 6
    status$complete_case_verified <- identical(read06("results/stage06_robustness/complete_case_check.tsv")$result, "IDENTICAL_SAMPLE")
    completed <- unlist(status[c("five_sensitivities_completed", "permutations_completed", "influence_completed",
      "random_gene_controls_completed", "arm_interaction_completed", "arm_specific_residuals_completed", "complete_case_verified")])
    assert06(all(completed), "Required analyses incomplete", "MODIFY_REVIEW")
    eligibility <- read06("results/stage06_robustness/stage06_eligibility.tsv")
    status$full_eligible_pathways <- sum(eligibility$full_transportability_eligible)
    status$hypothesis_generating_pathways <- sum(eligibility$classification == "HYPOTHESIS_GENERATING")
    status$not_robust_pathways <- sum(eligibility$classification == "NOT_ROBUST")
    gate <- jsonlite::read_json(path("results/stage06_robustness/external_transportability_gate.json"), simplifyVector = TRUE)
    status$external_transportability_decision <- gate$decision
    assert06(length(list.files(path("figures/stage06_robustness"), pattern = "[.](pdf|png)$")) == 8,
      "Stage 6 figures incomplete", "MODIFY_REVIEW")
    previous_lock <- renv::lockfile_read(path("metadata/outcome/stage05_environment_v2.lock"))
    renv::snapshot(project = root, lockfile = path("metadata/sensitivity/stage06_environment_v1.lock"),
      packages = unique(c(names(previous_lock$Packages), "GSVA")), prompt = FALSE)
    verify_prior06()
    qualifies <- any(eligibility$primary_direction_gate & eligibility$sensitivity_4_of_5_gate)
    status$decision <- if (qualifies) "PROCEED" else "STOP_REVIEW"
    status$decision_reason <- if (qualifies) paste("All planned Stage 6 analyses completed with zero technical failures.",
      "At least one pathway matches its prespecified direction and is stable in at least four of five sensitivities.",
      "Primary multiplicity results remain unchanged; no pathway meets full external-transportability eligibility.",
      "External transportability is skipped; no external cohort or Stage 7 analysis was opened.") else
      "No pathway matches its prespecified direction and remains stable in at least four of five sensitivity analyses."
  }, error = function(e) {
    status$technical_failures <<- status$technical_failures + 1L
    status$decision <<- if (inherits(e, "stage06_error")) e$decision else "STOP_REVIEW"
    status$decision_reason <<- conditionMessage(e)
  })
  packages <- unique(c("GSVA", "glmnet", "singscore", "stats",
    unlist(tools::package_dependencies("GSVA", db = installed.packages(), recursive = TRUE), use.names = FALSE)))
  versions <- data.frame(R_version = R.version.string, platform = R.version$platform,
    package = packages, package_version = vapply(packages, function(p) as.character(packageVersion(p)), ""))
  write06(versions, "metadata/sensitivity/software_versions_v1.tsv")
  old <- read_tsv(path("logs/software_versions.tsv"))
  write_tsv(rbind(old[!old$package %in% packages, ], versions), path("logs/software_versions.tsv"))
  jsonlite::write_json(status, path("results/stage06_robustness/stage06_status.json"),
    auto_unbox = TRUE, pretty = TRUE, na = "null")
  write.table(data.frame(timestamp = status$timestamp, stage = 6L, decision = status$decision,
    reason = status$decision_reason, inputs = spec06_file, outputs = "results/stage06_robustness/stage06_status.json"),
    path("logs/decision_log.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
})
files <- unlist(lapply(path(dirs06), list.files, recursive = TRUE, full.names = TRUE), use.names = FALSE)
code <- c("R/stage06_common.R", "R/sensitivity_residuals.R", "R/ssgsea_scoring.R", "R/composition_sensitivity.R",
  "R/permutation_tests.R", "R/influence_analysis.R", "R/random_gene_nulls.R", "R/treatment_sensitivity.R",
  "R/eligibility_gate.R", "R/stage06_validation.R",
  "scripts/18_freeze_stage06_spec.R", "scripts/19_generate_blinded_sensitivity_residuals.R",
  "scripts/20_run_stage06_outcome_sensitivities.R", "scripts/21_run_stage06_negative_controls.R",
  "scripts/22_evaluate_stage06_eligibility.R", "scripts/23_finalize_stage06.R",
  "scripts/test_stage06.R", "scripts/install_stage06.R", "scripts/audit_stage06_crossfit.R")
logs <- list.files(path("logs"), pattern = "^(18_|19_|20_|21_|22_|23_|stage06_).*\\.log$", full.names = TRUE)
files <- sort(unique(c(files, path(code), logs)))
files <- files[file.exists(files) & basename(files) != "stage06_output_checksums.tsv"]
write_tsv(data.frame(local_relative_path = substring(files, nchar(root) + 2L),
  file_size_bytes = file.info(files)$size, sha256 = vapply(files, sha256, "")),
  path("results/stage06_robustness/stage06_output_checksums.tsv"))
cat("results/stage06_robustness/stage06_status.json\n")
