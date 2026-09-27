source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/pathway_sets.R"))
stage03_run("06_freeze_pathways", {
  verify_stage02()
  blinded_identity()
  write_tsv(data.frame(layer = "broader_hallmark", status = "DEFERRED",
    reason = "Confirmatory residual analysis first"), path("metadata/pathways/exploratory_status.tsv"))
  freeze_panel()
})
cat("Panel frozen\n")
