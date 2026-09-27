source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/fold_generation.R"))
source(path("R/expected_response.R"))
source(path("R/residual_qc.R"))
stage04_init()
finalize04 <- function() {
  if (file.exists(path(stage04_freeze)) && file.exists(path("results/stage04_residuals/stage04_status.json"))) {
    verify_stage03()
    verify_residual_freeze04()
    verify_hash_table04("results/stage04_residuals/stage04_output_checksums.tsv")
    return(invisible(NULL))
  }
  status <- list(stage = 4L, timestamp = utc_now(), stage03_verified = FALSE,
    model_version = "v1", patients = 0L, pathways = 6L, outer_folds = 10L, inner_folds = 10L,
    crossfitted_predictions = 0L, patient_separation_failures = NA_integer_,
    leakage_failures = NA_integer_, residual_qc_failures = NA_integer_,
    pathways_beating_intercept = 0L, pathways_beating_arm_only = 0L,
    all_residuals_frozen = FALSE, residual_freeze_sha256 = NA_character_,
    decision = "STOP_REVIEW", decision_reason = "Stage 4 validation incomplete")
  tryCatch({
    status$stage03_verified <- verify_stage03()
    failure_file <- path("results/stage04_residuals/stage04_failure.json")
    if (file.exists(failure_file)) {
      failure <- jsonlite::read_json(failure_file, simplifyVector = TRUE)
      stop04(failure$reason, failure$decision)
    }
    x <- read_input04()
    s <- read_spec04()
    folds <- freeze_folds04(x, s)
    audited <- audit_saved04(x, s, folds)
    status$patients <- nrow(audited$wide)
    status$crossfitted_predictions <- nrow(audited$long)
    status$patient_separation_failures <- sum(as.matrix(audited$separation[, -(1:3)]))
    performance <- performance04(audited$benchmarks, s$axes)
    qc <- residual_qc04(audited$long, performance, s$axes)
    diagnostics <- fold_diagnostics04(audited$long, s)
    write_tsv(performance, path("results/stage04_residuals/expected_response_performance.tsv"))
    write_tsv(qc, path("results/stage04_residuals/residual_qc.tsv"))
    write_tsv(diagnostics, path("results/stage04_residuals/fold_diagnostics.tsv"))
    status$residual_qc_failures <- sum(!qc$passed)
    if (status$residual_qc_failures > 0L) stop04("Residual QC failed")
    comparisons <- do.call(rbind, lapply(s$axes, function(a) {
      z <- performance[performance$axis == a, ]
      mse <- setNames(z$oof_mse, z$method)
      data.frame(axis = a, full_ridge_mse = mse[["full_ridge"]],
        intercept_only_mse = mse[["intercept_only"]], arm_only_mse = mse[["arm_only"]],
        beats_intercept = mse[["full_ridge"]] < mse[["intercept_only"]],
        beats_arm_only = mse[["full_ridge"]] < mse[["arm_only"]],
        retained_in_confirmatory_family = TRUE)
    }))
    write_tsv(comparisons, path("results/stage04_residuals/benchmark_comparison.tsv"))
    status$pathways_beating_intercept <- sum(comparisons$beats_intercept)
    status$pathways_beating_arm_only <- sum(comparisons$beats_arm_only)
    leakage <- leakage04(s)
    status$leakage_failures <- sum(leakage$forbidden_fields)
    figures04(audited$long, diagnostics, s)
    files <- c(stage04_input, "results/stage03_pathways/stage03_status.json")
    write_tsv(data.frame(local_relative_path = files, sha256 = vapply(path(files), sha256, "")),
      path("metadata/modeling/stage04_input_checksums.tsv"))
    provenance <- list(timestamp = utc_now(), R_version = R.version.string, glmnet_version = s$glmnet_version,
      input_file = stage04_input, input_sha256 = sha256(path(stage04_input)),
      models = 60L, predictions = 900L, predictor_scaling = "Training data only, internal to each glmnet fit",
      outer_predictions = "Only held-out patients; lambda.1se",
      lambda_grid = "Derived by cv.glmnet from each outer training set; inner fits use their own training scaling",
      benchmarks = c("Outer-training mean delta", "Outer-training delta ~ arm linear model"),
      residual_scaling = "Mean and sample SD of the 150 cross-fitted residuals per axis",
      clinical_labels_accessed = FALSE, prior_stage_validation = "Status and checksum manifests; raw sources hashed only",
      panel_changes = "none", model_redesign = "none")
    yaml::write_yaml(provenance, path("metadata/modeling/stage04_provenance_v1.yml"))
    old_lock <- renv::lockfile_read(path("renv.lock"))
    renv::snapshot(project = root, lockfile = path("metadata/modeling/stage04_renv.lock"),
      packages = unique(c(names(old_lock$Packages), "glmnet")), prompt = FALSE)
    verify_stage03()
    freeze <- list(freeze_timestamp = utc_now(), model_spec_sha256 = sha256(path(stage04_spec)),
      outer_folds_sha256 = sha256(path(stage04_outer)), inner_folds_sha256 = sha256(path(stage04_inner)),
      long_residual_file_sha256 = sha256(path(stage04_long)), wide_residual_file_sha256 = sha256(path(stage04_wide)),
      scaling_file_sha256 = sha256(path(stage04_scaling)), number_patients = 150L, number_pathways = 6L,
      outcome_information_used = FALSE)
    if (file.exists(path(stage04_freeze))) {
      verify_residual_freeze04()
    } else yaml::write_yaml(freeze, path(stage04_freeze))
    leakage <- leakage04(s)
    status$leakage_failures <- sum(leakage$forbidden_fields)
    verify_residual_freeze04()
    status$all_residuals_frozen <- TRUE
    status$residual_freeze_sha256 <- sha256(path(stage04_freeze))
    if (status$pathways_beating_arm_only > 0L) {
      status$decision <- "PROCEED"
      status$decision_reason <- paste("All 150 patients have six verified outer-fold predictions;",
        "patient separation, leakage and residual QC pass. Frozen v1 residuals retain all six pathways.",
        status$pathways_beating_arm_only, "pathways improve on arm-only OOF MSE and",
        status$pathways_beating_intercept, "improve on intercept-only OOF MSE.",
        "Clinical outcomes were not accessed; no model redesign or Stage 5 analysis was performed.")
    } else {
      status$decision <- "MODIFY_REVIEW"
      status$decision_reason <- paste("Technical residual construction succeeds, but pretreatment molecular/composition state",
        "provides no observed out-of-fold improvement beyond treatment context, requiring methodology review",
        "before the residual construct is exposed to pCR. All six pathways retained; no automatic redesign.")
    }
  }, error = function(e) {
    status$decision <<- if (inherits(e, "stage04_error")) e$decision else "STOP_REVIEW"
    status$decision_reason <<- conditionMessage(e)
  })
  packages <- c("glmnet", "shape")
  packages <- packages[vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  existing <- read_tsv(path("logs/software_versions.tsv"))
  versions <- data.frame(R_version = R.version.string, platform = R.version$platform,
    package = packages, package_version = vapply(packages, function(p) as.character(packageVersion(p)), ""))
  write_tsv(rbind(existing[!existing$package %in% packages, ], versions), path("logs/software_versions.tsv"))
  jsonlite::write_json(status, path("results/stage04_residuals/stage04_status.json"),
    auto_unbox = TRUE, pretty = TRUE, na = "null")
  entry <- data.frame(timestamp = status$timestamp, stage = 4L, decision = status$decision,
    reason = status$decision_reason, inputs = stage04_input,
    outputs = "results/stage04_residuals/stage04_status.json;metadata/modeling/residual_freeze_v1.yml")
  write.table(entry, path("logs/decision_log.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE, append = TRUE)
}
run_stage("11_finalize_stage04", finalize04())
files <- c(unlist(lapply(path(stage04_dirs), list.files, full.names = TRUE, recursive = TRUE), use.names = FALSE),
  path(c("R/fold_generation.R", "R/expected_response.R", "R/residual_qc.R",
    "scripts/09_freeze_residual_model.R", "scripts/10_generate_crossfitted_residuals.R",
    "scripts/11_finalize_stage04.R", "scripts/test_stage04.R", "scripts/install_stage04.R",
    "logs/stage04_install.log", "logs/09_freeze_residual_model.log",
    "logs/10_generate_crossfitted_residuals.log", "logs/11_finalize_stage04.log", "logs/stage04_tests.log")))
files <- sort(unique(files[file.exists(files) & basename(files) != "stage04_output_checksums.tsv"]))
write_tsv(data.frame(local_relative_path = substring(files, nchar(root) + 2L),
  file_size_bytes = file.info(files)$size, sha256 = vapply(files, sha256, "")),
  path("results/stage04_residuals/stage04_output_checksums.tsv"))
cat("results/stage04_residuals/stage04_status.json\n")
