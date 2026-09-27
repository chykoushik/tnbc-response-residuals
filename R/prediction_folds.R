axes07 <- c("immune_interferon", "proliferation", "stromal_emt", "dna_repair", "hypoxia_stress", "apoptosis")
dirs07 <- c("metadata/prediction", "data/processed/prediction", "models/stage07",
  "results/stage07_prediction", "figures/stage07_prediction", "tables/stage07_prediction")
spec07_file <- "metadata/prediction/stage07_spec_v1.yml"
outer07_file <- "metadata/prediction/stage07_outer_folds_v1.tsv"
input07_file <- "data/processed/pathways/neotrip_model_inputs_blinded_v1.tsv"
expected07 <- c("arm", paste0("baseline_", axes07),
  "baseline_cancer_epithelial", "baseline_stromal_fraction", "baseline_immune_fraction")
benchmark07 <- c("arm", "baseline_cancer_epithelial", "baseline_immune_fraction", "baseline_stromal_fraction",
  "delta_cancer_epithelial", "delta_immune_fraction", "delta_stromal_fraction")
features07 <- list(P0 = benchmark07, P1 = c(benchmark07, paste0("baseline_", axes07)),
  P2 = c(benchmark07, paste0("baseline_", axes07), paste0("delta_", axes07)),
  P3 = c(benchmark07, paste0("baseline_", axes07), paste0("residual_", axes07)),
  P4 = c(benchmark07, paste0("baseline_", axes07), paste0("delta_", axes07), paste0("residual_", axes07)))
fail07 <- function(message, decision = "STOP_REVIEW") {
  stop(structure(list(message = message, call = NULL, decision = decision),
    class = c("stage07_error", "error", "condition")))
}
assert07 <- function(ok, message, decision = "STOP_REVIEW") if (!isTRUE(ok)) fail07(message, decision)
init07 <- function() for (d in dirs07) dir.create(path(d), recursive = TRUE, showWarnings = FALSE)
read07 <- function(file, blinded = FALSE) {
  header <- names(read.delim(path(file), nrows = 0, check.names = FALSE))
  if (blinded) assert07(!any(grepl("pcr|response|outcome|residual_disease", header, ignore.case = TRUE)),
    "Forbidden expected-response input field")
  read.delim(path(file), sep = "\t", quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
}
restricted07 <- function(file, fields) {
  header <- names(read.delim(path(file), nrows = 0, check.names = FALSE))
  assert07(all(fields %in% header), "Required identity fields missing")
  x <- read.delim(path(file), sep = "\t", quote = "", comment.char = "", check.names = FALSE,
    colClasses = ifelse(header %in% fields, "character", "NULL"), na.strings = NULL)
  x[, fields, drop = FALSE]
}
write07 <- function(x, file) {
  if (grepl("[.]gz$", file)) {
    con <- gzfile(path(file), "wt")
    tryCatch(write_tsv(x, con), finally = close(con))
  } else write_tsv(x, path(file))
}
verify_prior07 <- function() {
  files <- c("metadata/composition/stage01_input_checksums.tsv",
    "results/stage02_composition/stage02_output_checksums.tsv",
    "results/stage03_pathways/stage03_output_checksums.tsv",
    "results/stage04_residuals/stage04_output_checksums.tsv",
    "results/stage05_primary/stage05_output_checksums.tsv",
    "results/stage05_primary/stage05_output_checksums_v2.tsv",
    "results/stage06_robustness/stage06_output_checksums.tsv")
  for (file in files) {
    tab <- read_tsv(path(file)); names <- if ("path" %in% names(tab)) tab$path else tab$local_relative_path
    assert07(all(file.exists(path(names))) && identical(unname(vapply(path(names), sha256, "")), tab$sha256),
      "Prior-stage checksum mismatch")
  }
  invisible(TRUE)
}
gate07 <- function() {
  verify_prior07()
  s <- jsonlite::read_json(path("results/stage06_robustness/stage06_status.json"), simplifyVector = TRUE)
  gate <- jsonlite::read_json(path("results/stage06_robustness/external_transportability_gate.json"), simplifyVector = TRUE)
  assert07(s$decision == "PROCEED" && s$technical_failures == 0 && s$five_sensitivities_completed &&
    s$permutations_completed && s$influence_completed && s$random_gene_controls_completed &&
    s$complete_case_verified && s$external_transportability_decision == "SKIP" &&
    gate$eligible_count == 0 && gate$decision == "SKIP", "Stage 6 prerequisite failed")
  TRUE
}
freeze07 <- function(file, value, yaml = FALSE) {
  side <- paste0(file, ".sha256")
  if (file.exists(path(file)) || file.exists(path(side))) {
    assert07(file.exists(path(file)) && file.exists(path(side)) &&
      identical(readLines(path(side)), sha256(path(file))), "Frozen Stage 7 file changed")
    return(invisible(NULL))
  }
  assert07(!is.null(value), "Required frozen Stage 7 file missing")
  if (yaml) yaml::write_yaml(value, path(file)) else write07(value, file)
  writeLines(sha256(path(file)), path(side))
}
read_spec07 <- function() {
  assert07(file.exists(path(spec07_file)) &&
    identical(readLines(path(paste0(spec07_file, ".sha256"))), sha256(path(spec07_file))), "Prediction specification changed")
  s <- yaml::read_yaml(path(spec07_file))
  assert07(identical(s$features, features07) && s$outer_folds == 5L && s$outer_repeats == 10L &&
    identical(s$expected_response_predictors, expected07), "Prediction specification differs")
  s
}
stratify07 <- function(ids, strata, folds, seed) {
  assert07(!anyDuplicated(ids) && !anyNA(strata) && all(table(strata) >= folds), "Stratification invalid")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  assignments <- integer(length(ids))
  for (stratum in sort(unique(as.character(strata)), method = "radix")) {
    i <- which(as.character(strata) == stratum)
    i <- i[order(ids[i], method = "radix")]
    assignments[i] <- sample(rep(seq_len(folds), length.out = length(i)))
  }
  assignments
}
split_seed07 <- function(repetition, fold) 319650L + 10000L * repetition + 100L * fold
blinded07 <- function() {
  x <- read07(input07_file, blinded = TRUE)
  assert07(nrow(x) == 150 && !anyDuplicated(x$patient_id) && setequal(x$arm, c("CT", "CT/A")) &&
    all(is.finite(as.matrix(x[, -(1:2)]))), "Blinded prediction cohort invalid")
  x
}
read_outer07 <- function() {
  freeze07(outer07_file, NULL)
  folds <- read07(outer07_file, blinded = TRUE)
  x <- blinded07()
  assert07(nrow(folds) == 1500 && all(table(folds$patient_id) == 10L) &&
    !anyDuplicated(paste(folds[["repeat"]], folds$patient_id)), "Outer assignment coverage differs")
  for (r in 1:10) {
    f <- folds[folds[["repeat"]] == r, ]
    assert07(identical(f$patient_id, x$patient_id) && identical(f$arm, x$arm) &&
      setequal(f$outer_fold, 1:5), "Outer fold identity differs")
  }
  folds
}
population07 <- function() {
  blind <- blinded07()
  clinical <- restricted07("metadata/neotrip_pair_manifest_v1.tsv", c("patient_id", "arm", "pcr", "pcr_num", "paired"))
  clinical <- clinical[clinical$paired == "TRUE", ]
  stage4 <- restricted07("data/processed/residuals/neotrip_crossfitted_residuals_wide_v1.tsv", "patient_id")
  stage5 <- restricted07("data/processed/outcome/neotrip_primary_analysis_v2.tsv", c("patient_id", "arm", "pcr"))
  assert07(nrow(clinical) == 150 && nrow(stage4) == 150 && nrow(stage5) == 150 &&
    !anyDuplicated(clinical$patient_id) && !anyDuplicated(stage4$patient_id) &&
    !anyDuplicated(stage5$patient_id) && setequal(blind$patient_id, clinical$patient_id) &&
    setequal(blind$patient_id, stage4$patient_id) && setequal(blind$patient_id, stage5$patient_id),
    "Prediction cohort identity failed")
  clinical <- clinical[match(blind$patient_id, clinical$patient_id), ]
  stage5 <- stage5[match(blind$patient_id, stage5$patient_id), ]
  y <- unname(c(pCR = 1L, RD = 0L)[clinical$pcr])
  assert07(!anyNA(y) && all(clinical$pcr_num == as.character(y)) &&
    all(stage5$pcr == as.character(y)) && identical(clinical$arm, blind$arm) &&
    identical(stage5$arm, blind$arm), "Outcome coding mismatch")
  x <- blind
  x$pcr <- y
  for (f in c("baseline_cancer_epithelial", "baseline_immune_fraction", "baseline_stromal_fraction",
              "delta_cancer_epithelial", "delta_immune_fraction", "delta_stromal_fraction")) x[[f]] <- x[[f]] / 100
  counts <- c(patients = nrow(x), pcr = sum(y), RD = sum(y == 0),
    CT = sum(x$arm == "CT"), CTA = sum(x$arm == "CT/A"),
    missing = sum(!is.finite(as.matrix(x[, !names(x) %in% c("patient_id", "arm")]))))
  target <- c(150, 63, 87, 79, 71, 0)
  write07(data.frame(check = names(counts), observed = unname(counts), expected = target, passed = counts == target),
    "results/stage07_prediction/population_qc.tsv")
  assert07(all(counts == target), "Prediction population QC failed")
  x
}
run07 <- function(name, code) {
  init07()
  run_stage(name, {
    tryCatch({
      force(code)
      file <- path("results/stage07_prediction", paste0(name, "_failure.json"))
      if (file.exists(file)) unlink(file)
    }, error = function(e) {
      jsonlite::write_json(list(step = name, decision = if (inherits(e, "stage07_error")) e$decision else "STOP_REVIEW",
        reason = conditionMessage(e), timestamp = utc_now()),
        path("results/stage07_prediction", paste0(name, "_failure.json")), auto_unbox = TRUE, pretty = TRUE)
      stop(e)
    })
  })
}
