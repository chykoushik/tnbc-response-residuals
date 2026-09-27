source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/stage06_common.R"))
source(path("R/composition_sensitivity.R"))
source(path("R/treatment_sensitivity.R"))
run06("20_run_stage06_outcome_sensitivities", {
  verify_prior06()
  context <- analysis_context06()
  sensitivity_outcomes06(context)
  treatment_models06(context)
  missing <- sum(!is.finite(as.matrix(context$x[, !names(context$x) %in% c("patient_id", "arm")])))
  assert06(missing == 0 && nrow(context$x) == 150, "Complete-case sample changed", "MODIFY_REVIEW")
  write06(data.frame(complete_case_n = 150L, primary_n = 150L, missing_values = missing, result = "IDENTICAL_SAMPLE"),
    "results/stage06_robustness/complete_case_check.tsv")
  write06(data.frame(analysis = "alternative_deconvolution", status = "NOT_RUN",
    reason = "no outcome-independent alternative method was prespecified and validated before pCR analysis"),
    "results/stage06_robustness/deconvolution_sensitivity_status.tsv")
})
cat("Outcome sensitivities complete\n")
