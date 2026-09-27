stage04_axes <- c("immune_interferon", "proliferation", "stromal_emt",
  "dna_repair", "hypoxia_stress", "apoptosis")
stage04_input <- "data/processed/pathways/neotrip_model_inputs_blinded_v1.tsv"
stage04_spec <- "metadata/modeling/stage04_model_spec_v1.yml"
stage04_outer <- "metadata/modeling/stage04_outer_folds_v1.tsv"
stage04_inner <- "metadata/modeling/stage04_inner_folds_v1.tsv"
stage04_long <- "data/processed/residuals/neotrip_crossfitted_residuals_long_v1.tsv"
stage04_wide <- "data/processed/residuals/neotrip_crossfitted_residuals_wide_v1.tsv"
stage04_scaling <- "metadata/modeling/stage04_residual_scaling_v1.tsv"
stage04_freeze <- "metadata/modeling/residual_freeze_v1.yml"
stage04_dirs <- c("metadata/modeling", "data/processed/residuals", "models/stage04",
  "results/stage04_residuals", "figures/stage04_residuals")
stage04_init <- function() {
  for (d in stage04_dirs) dir.create(path(d), recursive = TRUE, showWarnings = FALSE)
}
stop04 <- function(message, decision = "STOP_REVIEW") {
  stop(structure(list(message = message, call = NULL, decision = decision),
    class = c("stage04_error", "error", "condition")))
}
forbidden04 <- function(fields) {
  fields[grepl("pcr|early[.]pcr|early_pcr|response|outcome|residual_disease",
    fields, ignore.case = TRUE)]
}
assert_fields04 <- function(fields) {
  if (length(forbidden04(fields))) stop04("Forbidden outcome fields")
}
verify_hash_table04 <- function(file) {
  tab <- read_tsv(path(file))
  files <- if ("path" %in% names(tab)) tab$path else tab$local_relative_path
  if (any(!file.exists(path(files))) ||
      !identical(unname(vapply(path(files), sha256, "")), tab$sha256))
    stop04("Prior stage provenance changed")
}
verify_stage03 <- function() {
  s <- jsonlite::read_json(path("results/stage03_pathways/stage03_status.json"), simplifyVector = TRUE)
  if (!identical(s$decision, "PROCEED") || !isTRUE(s$all_sets_frozen) ||
      !isTRUE(s$all_sets_pass_coverage) || !isTRUE(s$expression_samples == 401) ||
      !isTRUE(s$paired_patients == 150) || !isTRUE(s$composition_linked_patients == 150) ||
      !isTRUE(s$score_qc_failures == 0) || !isTRUE(s$mapping_failures == 0) ||
      !isTRUE(s$forbidden_outcome_fields == 0)) stop04("Stage 3 prerequisite failed")
  for (f in c("metadata/composition/stage01_input_checksums.tsv",
      "results/stage02_composition/stage02_output_checksums.tsv",
      "results/stage03_pathways/stage03_output_checksums.tsv",
      "metadata/pathways/implementation_checksums.tsv")) verify_hash_table04(f)
  TRUE
}
read_input04 <- function() {
  header <- names(read.delim(path(stage04_input), nrows = 0L, check.names = FALSE))
  assert_fields04(header)
  expected <- c("patient_id", "arm", paste0("baseline_", stage04_axes), paste0("delta_", stage04_axes),
    "baseline_cancer_epithelial", "delta_cancer_epithelial",
    "baseline_stromal_fraction", "delta_stromal_fraction", "baseline_immune_fraction", "delta_immune_fraction")
  if (!identical(header, expected)) stop04("Blinded input schema changed")
  x <- read_tsv(path(stage04_input))
  if (nrow(x) != 150L || anyDuplicated(x$patient_id) || any(missing_value(x$patient_id)) ||
      !setequal(x$arm, c("CT", "CT/A"))) stop04("Blinded patient identity invalid")
  for (f in names(x)[-(1:2)]) x[[f]] <- as.numeric(x[[f]])
  if (any(!is.finite(as.matrix(x[, -(1:2)])))) stop04("Blinded input values invalid")
  x
}
specification04 <- function() {
  list(model_version = "v1", model_family = "gaussian", method = "ridge", alpha = 0,
    lambda_rule = "lambda.1se", outer_folds = 10L, inner_folds = 10L, master_seed = 319641L,
    stratification = "treatment_arm", standardize = TRUE, intercept = TRUE, nlambda = 100L,
    lambda_min_ratio = 1e-4, type_measure = "mse",
    predictors = c("arm", paste0("baseline_", stage04_axes),
      "baseline_cancer_epithelial", "baseline_stromal_fraction", "baseline_immune_fraction"),
    outcomes = paste0("delta_", stage04_axes), axes = stage04_axes,
    arm_encoding = list(CT = 0L, "CT/A" = 1L), interaction_terms = "none",
    outcome_information_used = FALSE, global_predictor_scaling = "none",
    lambda_min_sensitivity = "DEFERRED", feature_selection = "none",
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    fold_algorithm = "Sort patient IDs within CT then CT/A; independently shuffle balanced repeated fold labels",
    inner_seed_rule = "master_seed + outer_fold",
    benchmark_r2 = "1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2)",
    cv_grouped = TRUE, cv_alignment = "lambda", cv_parallel = FALSE,
    glmnet_version = as.character(packageVersion("glmnet")),
    documentation = "https://search.r-project.org/CRAN/refmans/glmnet/html/cv.glmnet.html")
}
sidecar04 <- function(file) sub("[.](yml|tsv)$", ".sha256", file)
verify_frozen04 <- function(file) {
  checksum <- path(sidecar04(file))
  if (!file.exists(path(file)) || !file.exists(checksum) ||
      !identical(readLines(checksum, warn = FALSE), sha256(path(file))))
    stop04("Frozen file checksum failed")
}
freeze_spec04 <- function() {
  if (file.exists(path(stage04_spec)) || file.exists(path(sidecar04(stage04_spec)))) {
    return(read_spec04())
  }
  spec <- specification04()
  yaml::write_yaml(spec, path(stage04_spec))
  writeLines(sha256(path(stage04_spec)), path(sidecar04(stage04_spec)))
  read_spec04()
}
read_spec04 <- function() {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop04("glmnet unavailable")
  verify_frozen04(stage04_spec)
  s <- yaml::read_yaml(path(stage04_spec))
  if (!isTRUE(all.equal(s, specification04(), check.attributes = FALSE)))
    stop04("Frozen specification differs", "MODIFY_REVIEW")
  assert_fields04(c(s$predictors, s$outcomes))
  s
}
assign_folds04 <- function(patient_id, arm, k, seed, rng_kind) {
  if (anyDuplicated(patient_id) || !setequal(arm, c("CT", "CT/A")) ||
      any(table(arm) < k)) stop04("Fold identity invalid")
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  do.call(RNGkind, as.list(rng_kind))
  set.seed(seed)
  folds <- integer(length(patient_id))
  for (a in c("CT", "CT/A")) {
    ids <- which(arm == a)
    ids <- ids[order(patient_id[ids], method = "radix")]
    folds[ids] <- sample(rep(seq_len(k), length.out = length(ids)))
  }
  folds
}
validate_folds04 <- function(x, s, outer, inner) {
  assert_fields04(names(outer))
  assert_fields04(names(inner))
  if (!identical(names(outer), c("patient_id", "arm", "outer_fold")) ||
      !identical(names(inner), c("outer_fold", "patient_id", "arm", "inner_fold")) ||
      nrow(outer) != 150L || anyDuplicated(outer$patient_id) ||
      !identical(outer$patient_id, x$patient_id) || !identical(outer$arm, x$arm) ||
      !setequal(outer$outer_fold, seq_len(s$outer_folds)) || anyNA(outer$outer_fold))
    stop04("Outer fold identity failed")
  reports <- lapply(seq_len(s$outer_folds), function(k) {
    train <- outer[outer$outer_fold != k, , drop = FALSE]
    test <- outer[outer$outer_fold == k, , drop = FALSE]
    z <- inner[inner$outer_fold == k, , drop = FALSE]
    failures <- c(test_in_training = length(intersect(train$patient_id, test$patient_id)),
      test_in_inner = length(intersect(test$patient_id, z$patient_id)),
      inner_missing_training = length(setdiff(train$patient_id, z$patient_id)),
      inner_extra_patients = length(setdiff(z$patient_id, train$patient_id)),
      inner_duplicate_patients = sum(duplicated(z$patient_id)),
      outer_arms = as.integer(!setequal(test$arm, c("CT", "CT/A"))),
      inner_fold_values = as.integer(!setequal(z$inner_fold, seq_len(s$inner_folds)) || anyNA(z$inner_fold)),
      inner_arm_mismatch = sum(z$arm != train$arm[match(z$patient_id, train$patient_id)]))
    inner_arms <- vapply(seq_len(s$inner_folds),
      function(j) setequal(z$arm[z$inner_fold == j], c("CT", "CT/A")), logical(1))
    failures <- c(failures, inner_arms = sum(!inner_arms))
    data.frame(outer_fold = k, check = names(failures), violations = as.integer(failures))
  })
  report <- do.call(rbind, reports)
  if (anyNA(report$violations) || any(report$violations != 0L) ||
      nrow(inner) != 1350L || !setequal(inner$outer_fold, seq_len(s$outer_folds)))
    stop04("Nested patient separation failed")
  report
}
freeze_folds04 <- function(x, s) {
  if (file.exists(path(stage04_outer)) || file.exists(path(sidecar04(stage04_outer)))) {
    verify_frozen04(stage04_outer)
    outer <- read_tsv(path(stage04_outer))
    outer$outer_fold <- as.integer(outer$outer_fold)
  } else {
    outer <- x[, c("patient_id", "arm")]
    outer$outer_fold <- assign_folds04(x$patient_id, x$arm, s$outer_folds, s$master_seed, s$rng_kind)
    write_tsv(outer, path(stage04_outer))
    writeLines(sha256(path(stage04_outer)), path(sidecar04(stage04_outer)))
  }
  if (file.exists(path(stage04_inner)) || file.exists(path(sidecar04(stage04_inner)))) {
    verify_frozen04(stage04_inner)
    inner <- read_tsv(path(stage04_inner))
    inner$outer_fold <- as.integer(inner$outer_fold)
    inner$inner_fold <- as.integer(inner$inner_fold)
  } else {
    inner <- do.call(rbind, lapply(seq_len(s$outer_folds), function(k) {
      train <- outer[outer$outer_fold != k, c("patient_id", "arm")]
      data.frame(outer_fold = k, train,
        inner_fold = assign_folds04(train$patient_id, train$arm, s$inner_folds, s$master_seed + k, s$rng_kind))
    }))
    write_tsv(inner, path(stage04_inner))
    writeLines(sha256(path(stage04_inner)), path(sidecar04(stage04_inner)))
  }
  report <- validate_folds04(x, s, outer, inner)
  report_file <- path("results/stage04_residuals/inner_fold_separation_check.tsv")
  if (!file.exists(path(stage04_freeze))) write_tsv(report, report_file)
  list(outer = outer, inner = inner)
}
verify_residual_freeze04 <- function() {
  f <- yaml::read_yaml(path(stage04_freeze))
  mapping <- c(model_spec_sha256 = stage04_spec, outer_folds_sha256 = stage04_outer,
    inner_folds_sha256 = stage04_inner, long_residual_file_sha256 = stage04_long,
    wide_residual_file_sha256 = stage04_wide, scaling_file_sha256 = stage04_scaling)
  for (field in names(mapping)) if (!identical(f[[field]], sha256(path(mapping[[field]]))))
    stop04("Frozen residual bytes changed")
  if (!identical(f$outcome_information_used, FALSE) || f$number_patients != 150L ||
      f$number_pathways != 6L) stop04("Residual freeze metadata invalid")
  invisible(f)
}
run04 <- function(name, code) {
  stage04_init()
  run_stage(name, {
    tryCatch({
      force(code)
      failure <- path("results/stage04_residuals/stage04_failure.json")
      if (file.exists(failure) && identical(jsonlite::read_json(failure)$step, name)) unlink(failure)
    }, error = function(e) {
      jsonlite::write_json(list(step = name, timestamp = utc_now(),
        decision = if (inherits(e, "stage04_error")) e$decision else "STOP_REVIEW",
        reason = conditionMessage(e)), path("results/stage04_residuals/stage04_failure.json"),
        auto_unbox = TRUE, pretty = TRUE)
      stop(e)
    })
  })
}
