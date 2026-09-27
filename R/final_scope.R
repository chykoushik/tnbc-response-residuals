freeze_scope08 <- function() {
  verify08()
  s <- json08("results/stage08_final_gate/stage08_status.json")
  assert08(s$decision == "PROCEED_NEOTRIP_ONLY", "Scope freeze requires successful final protocol gate")
  scope <- list(version = "v1", primary_dataset = "GSE319641", analytic_patients = 150L,
    confirmatory_pathways = 6L, primary_scoring = "singscore", primary_residual_model = "Stage4_v1",
    primary_outcome_model = "Stage5_v2", robustness_stage = "Stage6_v1", prediction_stage = "Stage7_v1",
    external_transportability = "SKIPPED", GSE122630 = "NOT_ANALYSED", decision_curve = "NOT_RUN",
    exploratory_hallmark = "NOT_RUN", single_cell_support = "PENDING_STAGE09", manuscript_status = "NOT_STARTED",
    project_decision = s$decision, claims_subject_to_stage09_audit = TRUE)
  file <- path("metadata/final_gate/analysis_scope_freeze_v1.yml")
  if (file.exists(file)) {
    assert08(identical(yaml::read_yaml(file), scope) &&
      identical(readLines(paste0(file, ".sha256")), sha256(file)), "Scope freeze differs")
  } else {
    yaml::write_yaml(scope, file); writeLines(sha256(file), paste0(file, ".sha256"))
  }
  prohibited <- c("independent TNBC replication", "external validation of residual phenotype",
    "validated clinical predictor", "clinically useful classifier", "significant confirmatory residual pathway",
    "causal treatment mechanism", "tumour-intrinsic interpretation without qualification")
  potential <- c("cross-fitted residual molecular response is technically measurable",
    "residuals contain structured biological variation",
    "some prespecified residual pathways show directionally stable hypothesis-generating pCR associations",
    "composition adjustment materially affects residual associations",
    "residual features do not show consistent incremental predictive improvement",
    "findings are single-cohort and require independent validation")
  boundaries <- data.frame(claim = c(prohibited, potential),
    status = c(rep("NOT_SUPPORTED", 7), rep("DEFER_STAGE09", 6)),
    reason = c("No independent cohort analysed.", "External transportability eligibility failed.",
      "No external validation or consistent incremental gain.", "Clinical utility not established; DCA not run.",
      "Zero pathways pass mandatory multiplicity correction.", "Observational associations do not establish causality.",
      "Bulk transcriptomic residuals remain composition-sensitive.",
      rep("Potentially supportable only with qualification and after the final Stage 9 audit; not a manuscript claim.", 6)))
  write08(boundaries, "metadata/final_gate/claim_boundaries_v1.tsv")
  packages <- c("stats", "digest", "jsonlite", "yaml")
  versions <- data.frame(R_version = R.version.string, platform = R.version$platform,
    package = packages, package_version = vapply(packages, function(p) as.character(packageVersion(p)), ""))
  write08(versions, "metadata/final_gate/software_versions_v1.tsv")
  old <- read08("logs/software_versions.tsv")
  write08(rbind(old[!old$package %in% packages, ], versions), "logs/software_versions.tsv")
  write.table(data.frame(timestamp = s$timestamp, stage = 8L, decision = s$decision,
    reason = s$decision_reason, inputs = "Frozen Stage 5 v2, Stage 6 and Stage 7 outputs",
    outputs = "results/stage08_final_gate/stage08_status.json"),
    path("logs/decision_log.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  verify08()
}
