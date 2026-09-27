source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/outcome_manifest.R"))
source(path("R/primary_logistic.R"))
source(path("R/primary_inference.R"))
source(path("R/primary_figures.R"))
init05()
run_stage("13_run_primary_outcome_analysis", {
  state <- list(started = utc_now(), stage04_verified = FALSE, residual_freeze_verified = FALSE,
    specification_verified = FALSE, models_attempted = FALSE, decision = "STOP_REVIEW", reason = "Not validated")
  tryCatch({
    verify_stage04_05()
    state$stage04_verified <- TRUE
    state$residual_freeze_verified <- TRUE
    spec <- read_spec05()
    state$specification_verified <- TRUE
    state$specification_sha256 <- sha256(path(spec_file05))
    state$outcome_access_started <- utc_now()
    data <- prepare_analysis05(spec)
    composition_check05(data, spec)
    state$models_attempted <- TRUE
    fits <- fit_logistic05(data, spec)
    inference <- inference05(fits$models, spec)
    save_inference05(inference)
    figures05(inference, spec)
    state$decision <- "PROCEED"
    state$reason <- "Six primary models and six base comparisons pass technical checks; robustness remains Stage 6."
  }, error = function(e) {
    state$decision <<- if (inherits(e, "stage05_error")) e$decision else "STOP_REVIEW"
    state$reason <<- conditionMessage(e)
  })
  state$completed <- utc_now()
  jsonlite::write_json(state, path("metadata/outcome/stage05_execution_v1.json"), auto_unbox = TRUE, pretty = TRUE)
})
cat("Outcome checks complete\n")
