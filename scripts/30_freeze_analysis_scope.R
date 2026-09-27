source("R/paths.R")
source("R/checksums.R")
source("R/prediction_folds.R")
source("R/final_gate.R")
source("R/final_scope.R")
run_stage("30_freeze_analysis_scope", freeze_scope08())
files <- unlist(lapply(path(dirs08), list.files, recursive = TRUE, full.names = TRUE), use.names = FALSE)
files <- sort(unique(c(files, path(c("R/final_gate.R", "R/final_scope.R",
  "scripts/29_apply_final_protocol_gate.R", "scripts/30_freeze_analysis_scope.R",
  "logs/29_apply_final_protocol_gate.log", "logs/30_freeze_analysis_scope.log")))))
files <- files[basename(files) != "stage08_output_checksums.tsv"]
write08(data.frame(local_relative_path = substring(files, nchar(root) + 2L),
  file_size_bytes = file.info(files)$size, sha256 = vapply(files, sha256, "")),
  "results/stage08_final_gate/stage08_output_checksums.tsv")
cat("results/stage08_final_gate/stage08_status.json\n")
