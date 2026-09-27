source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/prediction_folds.R"))
run07("25_generate_stage07_folds", {
  verify_prior07()
  read_spec07()
  x <- blinded07()
  assignments <- do.call(rbind, lapply(1:10, function(r) {
    z <- data.frame(repetition = r, patient_id = x$patient_id, arm = x$arm,
      outer_fold = stratify07(x$patient_id, x$arm, 5, 319650 + r))
    names(z)[1] <- "repeat"
    z
  }))
  freeze07(outer07_file, assignments)
  outer <- read_outer07()
  residual_crossfit <- expected_tuning <- list()
  for (r in 1:10) for (k in 1:5) {
    f <- outer[outer[["repeat"]] == r, ]
    train <- x[f$outer_fold != k, ]
    seed <- split_seed07(r, k)
    cf <- stratify07(train$patient_id, train$arm, 5, seed + 1)
    row <- data.frame(repetition = r, outer_fold = k, patient_id = train$patient_id,
      arm = train$arm, residual_fold = cf)
    names(row)[1] <- "repeat"
    residual_crossfit[[paste(r, k)]] <- row
    for (c in 0:5) {
      ids <- if (c == 0) seq_len(nrow(train)) else which(cf != c)
      tuned <- data.frame(repetition = r, outer_fold = k, residual_fold = c,
        patient_id = train$patient_id[ids], arm = train$arm[ids],
        tuning_fold = stratify07(train$patient_id[ids], train$arm[ids], 5,
          seed + if (c == 0) 20 else 10 + c))
      names(tuned)[1] <- "repeat"
      expected_tuning[[paste(r, k, c)]] <- tuned
    }
  }
  freeze07("metadata/prediction/stage07_residual_crossfit_folds_v1.tsv", do.call(rbind, residual_crossfit))
  freeze07("metadata/prediction/stage07_expected_tuning_folds_v1.tsv", do.call(rbind, expected_tuning))
  population <- population07()
  audits <- classifier_inner <- list()
  for (r in 1:10) for (k in 1:5) {
    f <- outer[outer[["repeat"]] == r, ]
    test <- which(f$outer_fold == k); train <- which(f$outer_fold != k)
    row <- data.frame(repetition = r, outer_fold = k, n_test = length(test), pcr_n = sum(population$pcr[test]),
      rd_n = sum(population$pcr[test] == 0), both_classes = length(unique(population$pcr[test])) == 2)
    names(row)[1] <- "repeat"
    audits[[paste(r, k)]] <- row
    row <- data.frame(repetition = r, outer_fold = k, patient_id = population$patient_id[train],
      inner_fold = stratify07(population$patient_id[train], population$pcr[train], 5, split_seed07(r, k) + 30))
    names(row)[1] <- "repeat"
    classifier_inner[[paste(r, k)]] <- row
  }
  audit <- do.call(rbind, audits)
  write07(audit, "results/stage07_prediction/outer_fold_outcome_audit.tsv")
  assert07(all(audit$both_classes), "Frozen outer test fold has one class", "MODIFY_REVIEW")
  freeze07("metadata/prediction/stage07_classifier_inner_folds_v1.tsv", do.call(rbind, classifier_inner))
  write07(population, "data/processed/prediction/stage07_analysis_input_v1.tsv")
})
cat("Prediction folds frozen\n")
