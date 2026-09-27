source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/fold_generation.R"))
source(path("R/expected_response.R"))
run04("10_generate_crossfitted_residuals", {
  verify_stage03()
  x <- read_input04()
  s <- read_spec04()
  folds <- freeze_folds04(x, s)
  if (file.exists(path(stage04_freeze))) {
    verify_residual_freeze04()
  } else {
    fit_models04(x, s, folds)
  }
})
cat("Residual construction complete\n")
