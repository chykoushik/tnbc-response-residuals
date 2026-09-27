source("R/paths.R")
source("R/checksums.R")
source("R/prediction_folds.R")
source("R/final_gate.R")
for (d in dirs08) dir.create(path(d), recursive = TRUE, showWarnings = FALSE)
run_stage("29_apply_final_protocol_gate", {
  status <- tryCatch(apply_gate08(), error = function(e) list(stage = 8L, timestamp = utc_now(),
    stage05_verified = FALSE, stage06_verified = FALSE, stage07_verified = FALSE,
    patients = NA_integer_, confirmatory_pathways = 6L, multiplicity_confirmed_pathways = NA_integer_,
    full_transportability_eligible_pathways = NA_integer_, GSE122630_status = "SKIP",
    consistent_incremental_prediction_gain = FALSE, fundamental_stop_condition = TRUE,
    decision = "STOP_REVIEW", decision_reason = conditionMessage(e)))
  jsonlite::write_json(status, path("results/stage08_final_gate/stage08_status.json"),
    auto_unbox = TRUE, pretty = TRUE, na = "null")
})
