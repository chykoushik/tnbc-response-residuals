source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/prediction_folds.R"))
source(path("R/nested_prediction_residuals.R"))
run07("26_generate_nested_prediction_residuals", {
  verify_prior07()
  read_spec07()
  if (file.exists(path("metadata/prediction/stage07_residual_freeze_v1.yml"))) {
    verify_prediction_residuals07()
  } else {
    x <- blinded07()
    outer <- read_outer07()
    files <- c("metadata/prediction/stage07_residual_crossfit_folds_v1.tsv",
      "metadata/prediction/stage07_expected_tuning_folds_v1.tsv")
    for (f in files) freeze07(f, NULL)
    cf <- read07(files[1], blinded = TRUE)
    tuning <- read07(files[2], blinded = TRUE)
    result <- list()
    for (r in 1:10) for (k in 1:5) {
      f <- outer[outer[["repeat"]] == r, ]
      c <- cf[cf[["repeat"]] == r & cf$outer_fold == k, ]
      t <- tuning[tuning[["repeat"]] == r & tuning$outer_fold == k, ]
      result[[paste(r, k)]] <- residual_split07(x, f, c, t, r, k)
      cat("Residual split", r, k, "\n"); flush.console()
    }
    write07(do.call(rbind, lapply(result, "[[", "residuals")), "data/processed/prediction/stage07_nested_residuals_v1.tsv.gz")
    write07(do.call(rbind, lapply(result, "[[", "audit")), "results/stage07_prediction/prediction_residual_audit.tsv")
    write07(do.call(rbind, lapply(result, "[[", "scaling")), "metadata/prediction/stage07_residual_scaling_v1.tsv")
    write07(do.call(rbind, lapply(result, "[[", "fits")), "results/stage07_prediction/expected_response_tuning.tsv")
    saved <- c("data/processed/prediction/stage07_nested_residuals_v1.tsv.gz",
      "metadata/prediction/stage07_residual_scaling_v1.tsv", outer07_file, files)
    freeze07("metadata/prediction/stage07_residual_freeze_v1.yml",
      list(timestamp = utc_now(), clinical_labels_used = FALSE,
        specification_sha256 = sha256(path(spec07_file)),
        files = data.frame(local_relative_path = saved, sha256 = unname(vapply(path(saved), sha256, "")))), yaml = TRUE)
  }
})
cat("Prediction residuals frozen\n")
