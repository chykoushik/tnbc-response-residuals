axes06 <- c("immune_interferon", "proliferation", "stromal_emt", "dna_repair", "hypoxia_stress", "apoptosis")
variants06 <- c("S1_ssgsea", "S2_no_baseline_composition", "S3_lambda_min", "S4_arm_baseline_interaction", "S5_ilr_composition")
dirs06 <- c("metadata/sensitivity", "data/processed/sensitivity", "models/stage06",
  "results/stage06_robustness", "figures/stage06_robustness", "tables/stage06_robustness")
base_predictors06 <- c("arm", paste0("baseline_", axes06),
  "baseline_cancer_epithelial", "baseline_stromal_fraction", "baseline_immune_fraction")
spec06_file <- "metadata/sensitivity/stage06_spec_v1.yml"
freeze06_file <- "metadata/sensitivity/stage06_blinded_residual_freeze_v1.yml"
init06 <- function() for (d in dirs06) dir.create(path(d), recursive = TRUE, showWarnings = FALSE)
fail06 <- function(message, decision = "MODIFY_REVIEW") {
  stop(structure(list(message = message, call = NULL, decision = decision),
    class = c("stage06_error", "error", "condition")))
}
assert06 <- function(ok, message, decision = "STOP_REVIEW") if (!isTRUE(ok)) fail06(message, decision)
forbidden06 <- function(fields) fields[grepl("pcr|early[.]pcr|early_pcr|response|outcome|residual_disease", fields, ignore.case = TRUE)]
read06 <- function(file, blinded = FALSE) {
  header <- names(read.delim(path(file), nrows = 0L, check.names = FALSE))
  if (blinded) assert06(!length(forbidden06(header)), "Forbidden blinded input fields")
  read.delim(path(file), sep = "\t", quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
}
write06 <- function(x, file) {
  if (grepl("[.]gz$", file)) {
    con <- gzfile(path(file), "wt")
    tryCatch(write_tsv(x, con), finally = close(con))
  } else write_tsv(x, path(file))
}
hashes06 <- function(files) data.frame(local_relative_path = files, sha256 = unname(vapply(path(files), sha256, "")))
check_hashes06 <- function(tab) {
  fields <- if ("path" %in% names(tab)) tab$path else tab$local_relative_path
  assert06(all(file.exists(path(fields))), "Frozen files missing")
  assert06(identical(unname(vapply(path(fields), sha256, "")), as.character(tab$sha256)), "Frozen files changed")
}
verify_prior06 <- function() {
  for (f in c("metadata/composition/stage01_input_checksums.tsv",
    "results/stage02_composition/stage02_output_checksums.tsv",
    "results/stage03_pathways/stage03_output_checksums.tsv",
    "results/stage04_residuals/stage04_output_checksums.tsv",
    "results/stage05_primary/stage05_output_checksums.tsv",
    "results/stage05_primary/stage05_output_checksums_v2.tsv")) check_hashes06(read_tsv(path(f)))
  invisible(TRUE)
}
stage05_gate06 <- function() {
  verify_prior06()
  s <- jsonlite::read_json(path("results/stage05_primary/stage05_status_v2.json"), simplifyVector = TRUE)
  assert06(s$decision == "PROCEED" && s$patients == 150 && s$primary_models == 6 &&
    s$models_converged == 6 && s$rank_failures == 0 && s$separation_failures == 0 &&
    s$missing_data_failures == 0 && s$identity_failures == 0 &&
    isTRUE(s$residual_freeze_verified) && isTRUE(s$four_component_constraint_verified), "Stage 5 gate failed")
  TRUE
}
freeze_yaml06 <- function(value, file) {
  side <- sub("[.]yml$", ".sha256", file)
  if (file.exists(path(file))) {
    assert06(file.exists(path(side)) && identical(readLines(path(side)), sha256(path(file))), "Frozen metadata changed")
    return(yaml::read_yaml(path(file)))
  }
  yaml::write_yaml(value, path(file))
  writeLines(sha256(path(file)), path(side))
  value
}
read_spec06 <- function() {
  assert06(file.exists(path(spec06_file)) &&
    identical(readLines(path("metadata/sensitivity/stage06_spec_v1.sha256")), sha256(path(spec06_file))), "Stage 6 specification changed")
  s <- yaml::read_yaml(path(spec06_file))
  assert06(identical(s$primary_robustness_variants, variants06) && s$permutations == 10000 &&
    s$random_gene_sets_per_axis == 100, "Stage 6 specification differs")
  s
}
run06 <- function(name, code) {
  init06()
  run_stage(name, {
    tryCatch({
      force(code)
      failure <- path("results/stage06_robustness", paste0(name, "_failure.json"))
      if (file.exists(failure)) unlink(failure)
    }, error = function(e) {
      jsonlite::write_json(list(step = name, decision = if (inherits(e, "stage06_error")) e$decision else "STOP_REVIEW",
        reason = conditionMessage(e), timestamp = utc_now()),
        path("results/stage06_robustness", paste0(name, "_failure.json")), auto_unbox = TRUE, pretty = TRUE)
      stop(e)
    })
  })
}
blinded_context06 <- function() {
  files <- c("data/processed/pathways/neotrip_model_inputs_blinded_v1.tsv",
    "metadata/modeling/stage04_outer_folds_v1.tsv", "metadata/modeling/stage04_inner_folds_v1.tsv",
    "data/processed/composition/neotrip_pair_composition_change.tsv",
    "data/processed/composition/neotrip_composition.tsv")
  z <- lapply(files, read06, blinded = TRUE)
  x <- z[[1]]
  assert06(nrow(x) == 150 && !anyDuplicated(x$patient_id), "Blinded patient identity invalid")
  assert06(identical(x$patient_id, z[[2]]$patient_id) && identical(x$arm, z[[2]]$arm), "Outer fold identity differs")
  p <- z[[4]][match(x$patient_id, z[[4]]$patient_id), ]
  assert06(!anyNA(p$patient_id) && identical(p$patient_id, x$patient_id), "Pair identity differs")
  write06(data.frame(file = files, forbidden_fields = 0L, passed = TRUE), "results/stage06_robustness/blinded_input_leakage_check.tsv")
  list(x = x, outer = z[[2]], inner = z[[3]], pairs = p, identity = z[[5]][, 1:4])
}
verify_blind_freeze06 <- function() {
  read_spec06()
  assert06(file.exists(path(freeze06_file)), "Blinded residual freeze missing")
  assert06(identical(readLines(path(sub("[.]yml$", ".sha256", freeze06_file))), sha256(path(freeze06_file))), "Blinded freeze changed")
  f <- yaml::read_yaml(path(freeze06_file))
  check_hashes06(as.data.frame(f$files, stringsAsFactors = FALSE))
  assert06(identical(f$clinical_labels_accessed, FALSE), "Blinded construction leakage")
}
analysis_context06 <- function() {
  verify_blind_freeze06()
  x <- read06("data/processed/outcome/neotrip_primary_analysis_v2.tsv")
  primary <- read06("results/stage05_primary/primary_residual_associations_v2.tsv")
  assert06(nrow(x) == 150 && sum(x$pcr == 1) == 63 && sum(x$pcr == 0) == 87 &&
    !anyDuplicated(x$patient_id) && identical(primary$axis, axes06), "Outcome identity invalid")
  list(x = x, primary = primary)
}
glm_effect06 <- function(fit, term = "residual") {
  b <- unname(coef(fit)[term]); se <- sqrt(vcov(fit)[term, term])
  data.frame(beta = b, standard_error = se, odds_ratio = exp(b),
    ci95_lower = exp(b - qnorm(.975) * se), ci95_upper = exp(b + qnorm(.975) * se),
    p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE), direction = ifelse(b > 0, "positive", ifelse(b < 0, "negative", "zero")))
}
outcome_frame06 <- function(x, a) {
  data.frame(pcr = x$pcr, residual = x[[paste0("residual_z_", a)]],
    arm = as.integer(x$arm == "CT/A"), baseline = x[[paste0("baseline_z_", a)]],
    cancer = x$delta_cancer_epithelial_10, immune = x$delta_immune_fraction_10,
    stromal = x$delta_stromal_fraction_10)
}
fit_glm06 <- function(frame, formula = pcr ~ residual + arm + baseline + cancer + immune + stromal) {
  warnings <- character()
  fit <- withCallingHandlers(glm(formula, data = frame, family = binomial("logit"), na.action = na.fail),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  assert06(isTRUE(fit$converged) && fit$rank == length(coef(fit)) &&
    all(is.finite(coef(fit))) && all(is.finite(sqrt(diag(vcov(fit))))) && !length(warnings),
    "Logistic sensitivity failed", "MODIFY_REVIEW")
  fit
}
