source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/outcome_manifest.R"))
source(path("R/primary_inference.R"))
source(path("R/stage05_v2.R"))
run_stage("16_run_primary_outcome_v2", {
  state <- list(started = utc_now(), stage04_verified = FALSE, residual_freeze_verified = FALSE,
    four_component_constraint_verified = FALSE, models_attempted = FALSE,
    decision = "STOP_REVIEW", reason = "Not validated")
  tryCatch({
    verify_v2()
    state$stage04_verified <- TRUE
    state$residual_freeze_verified <- TRUE
    s <- read_spec_v2()
    state$specification_sha256 <- sha256(path(spec_v2))
    state$amendment_sha256 <- sha256(path(amend_v2))
    x <- prepare_v2(s)
    state$four_component_constraint_verified <- TRUE
    state$models_attempted <- TRUE
    fits <- fit_v2(x, s)
    out <- infer_v2(fits, s)
    save_inference_v2(out)
    figures_v2(out, s)
    state$decision <- "PROCEED"
    state$reason <- "Six primary and six base models pass technical checks; Stage 6 remains deferred."
  }, error = function(e) {
    state$decision <<- if (inherits(e, "stage05_error")) e$decision else "STOP_REVIEW"
    state$reason <<- conditionMessage(e)
  })
  state$completed <- utc_now()
  jsonlite::write_json(state, path("metadata/outcome/stage05_execution_v2.json"), auto_unbox = TRUE, pretty = TRUE)
})
cat("V2 analysis complete\n")
