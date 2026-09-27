source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/outcome_manifest.R"))
source(path("R/primary_logistic.R"))
source(path("R/primary_inference.R"))
init05()
run_stage("14_finalize_stage05", {
  status <- list(stage = 5L, timestamp = utc_now(), stage04_verified = FALSE,
    residual_freeze_verified = FALSE, patients = 0L, pcr = 0L, residual_disease = 0L,
    treatment_CT = 0L, treatment_CTA = 0L, primary_models = 0L, models_converged = 0L,
    rank_failures = NA_integer_, missing_data_failures = NA_integer_, identity_failures = NA_integer_,
    direction_matching_pathways = NA_integer_, effect_size_gate_pathways = NA_integer_,
    multiplicity_gate_pathways = NA_integer_, both_effect_and_multiplicity_gate = NA_integer_,
    decision = "STOP_REVIEW", decision_reason = "Stage 5 validation incomplete")
  tryCatch({
    verify_stage04_05()
    status$stage04_verified <- TRUE
    status$residual_freeze_verified <- TRUE
    spec <- read_spec05()
    state <- jsonlite::read_json(path("metadata/outcome/stage05_execution_v1.json"), simplifyVector = TRUE)
    if (!identical(state$specification_sha256, sha256(path(spec_file05))))
      stop05("Primary specification changed after outcome access")
    data <- read_numeric05(analysis_file05, c("patient_id", "arm"))
    qc <- read_tsv(path("results/stage05_primary/outcome_qc.tsv"))
    missing <- read_tsv(path("results/stage05_primary/missing_data_check.tsv"))
    status$patients <- nrow(data)
    status$pcr <- sum(data$pcr == 1)
    status$residual_disease <- sum(data$pcr == 0)
    status$treatment_CT <- sum(data$arm == "CT")
    status$treatment_CTA <- sum(data$arm == "CT/A")
    status$identity_failures <- sum(as.integer(qc$failures))
    status$missing_data_failures <- sum(as.integer(missing$patient_count))
    if (status$identity_failures != 0L || status$patients != 150L || status$pcr != 63L ||
        status$residual_disease != 87L || status$treatment_CT != 79L || status$treatment_CTA != 71L)
      stop05("Outcome count or identity validation failed")
    if (status$missing_data_failures != 0L) stop05("Unexpected missingness", "MODIFY_REVIEW")
    if (state$decision != "PROCEED") stop05(state$reason, state$decision)
    composition_check05(data, spec)
    diagnostics <- read_tsv(path("results/stage05_primary/model_diagnostics.tsv"))
    primary <- diagnostics[diagnostics$model == "composition_adjusted_model", ]
    status$primary_models <- nrow(primary)
    status$models_converged <- sum(primary$converged == "TRUE")
    status$rank_failures <- sum(as.integer(primary$rank) != as.integer(primary$number_parameters))
    if (nrow(diagnostics) != 12L || status$primary_models != 6L || any(diagnostics$passed != "TRUE"))
      stop05("Model diagnostics failed", "MODIFY_REVIEW")
    models <- readRDS(path("results/stage05_primary/primary_models_v1.rds"))
    expected <- inference05(models, spec)
    stored <- read_tsv(path("results/stage05_primary/primary_residual_associations.tsv"))
    comparison <- read_tsv(path("results/stage05_primary/composition_adjustment_comparison.tsv"))
    compare <- function(saved, calculated) {
      if (!identical(names(saved), names(calculated)) || nrow(saved) != nrow(calculated)) return(FALSE)
      for (field in names(saved)) {
        if (is.numeric(calculated[[field]])) {
          if (!isTRUE(all.equal(as.numeric(saved[[field]]), calculated[[field]], tolerance = 1e-10,
              check.attributes = FALSE))) return(FALSE)
        } else if (!identical(saved[[field]], as.character(calculated[[field]]))) return(FALSE)
      }
      TRUE
    }
    if (!compare(stored, expected$primary) || !compare(comparison, expected$comparison))
      stop05("Stored inference cannot be reproduced")
    status$direction_matching_pathways <- sum(expected$primary$direction_matches)
    status$effect_size_gate_pathways <- sum(expected$primary$effect_size_gate)
    status$multiplicity_gate_pathways <- sum(expected$primary$multiplicity_gate)
    status$both_effect_and_multiplicity_gate <- sum(expected$primary$effect_size_gate & expected$primary$multiplicity_gate)
    status$decision <- "PROCEED"
    status$decision_reason <- "All six prespecified primary models pass technical validation; freeze intact and all 150 patients match. Results do not establish external-validation eligibility; Stage 6 robustness remains deferred."
  }, error = function(e) {
    status$decision <<- if (inherits(e, "stage05_error")) e$decision else "STOP_REVIEW"
    status$decision_reason <<- conditionMessage(e)
  })
  if (status$primary_models == 0L) {
    deferred <- c("results/stage05_primary/primary_residual_associations.tsv",
      "results/stage05_primary/composition_adjustment_comparison.tsv",
      "results/stage05_primary/model_diagnostics.tsv", "results/stage05_primary/primary_models_v1.rds",
      "tables/stage05_primary_residual_associations.tsv",
      "figures/stage05_primary/primary_residual_forest.pdf",
      "figures/stage05_primary/primary_residual_forest.png",
      "figures/stage05_primary/composition_adjustment_comparison.pdf",
      "figures/stage05_primary/composition_adjustment_comparison.png")
    write_tsv(data.frame(file = deferred, status = "NOT_RUN", reason = status$decision_reason),
      path("results/stage05_primary/deferred_outputs.tsv"))
  }
  packages <- c("stats", "graphics", "grDevices", "jsonlite", "yaml", "digest")
  existing <- read_tsv(path("logs/software_versions.tsv"))
  versions <- data.frame(R_version = R.version.string, platform = R.version$platform,
    package = packages, package_version = vapply(packages, function(p) as.character(packageVersion(p)), ""))
  write_tsv(rbind(existing[!existing$package %in% packages, ], versions), path("logs/software_versions.tsv"))
  jsonlite::write_json(status, path("results/stage05_primary/stage05_status.json"),
    auto_unbox = TRUE, pretty = TRUE, na = "null")
  entry <- data.frame(timestamp = status$timestamp, stage = 5L, decision = status$decision,
    reason = status$decision_reason,
    inputs = "Frozen Stage 4 residuals;Stage 1 patient manifest;Stage 3 blinded input",
    outputs = "results/stage05_primary/stage05_status.json")
  write.table(entry, path("logs/decision_log.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE, append = TRUE)
})
files <- c(unlist(lapply(path(dirs05), list.files, full.names = TRUE, recursive = TRUE), use.names = FALSE),
  path(c("R/outcome_manifest.R", "R/primary_logistic.R", "R/primary_inference.R", "R/primary_figures.R",
    "scripts/12_freeze_primary_outcome_spec.R", "scripts/13_run_primary_outcome_analysis.R",
    "scripts/14_finalize_stage05.R", "scripts/test_stage05.R", "scripts/audit_stage05_composition.R",
    "logs/12_freeze_primary_outcome_spec.log", "logs/13_run_primary_outcome_analysis.log",
    "logs/14_finalize_stage05.log", "logs/stage05_tests.log",
    "tables/stage05_primary_residual_associations.tsv")))
files <- sort(unique(files[file.exists(files) & basename(files) != "stage05_output_checksums.tsv"]))
write_tsv(data.frame(local_relative_path = substring(files, nchar(root) + 2L),
  file_size_bytes = file.info(files)$size, sha256 = vapply(files, sha256, "")),
  path("results/stage05_primary/stage05_output_checksums.tsv"))
cat("results/stage05_primary/stage05_status.json\n")
