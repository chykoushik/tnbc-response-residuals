source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/repository_audit.R"))
source(path("R/composition_checks.R"))
source(path("R/instaprism_runner.R"))
run_stage("04_composition", {
  authority <- verify_stage01()
  feasibility <- audit_instaprism()
  if (isTRUE(feasibility$reproducible_from_public_inputs)) {
    run_instaprism(authority)
  } else {
    write_tsv(data.frame(check = "reproduction_not_run", failures = NA_integer_, passed = FALSE,
      detail = feasibility$reason), path("results/stage02_composition/composition_qc.tsv"))
  }
})
cat("Composition step complete\n")
