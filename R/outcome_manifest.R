axes05 <- c("immune_interferon", "proliferation", "stromal_emt", "dna_repair", "hypoxia_stress", "apoptosis")
dirs05 <- c("metadata/outcome", "data/processed/outcome", "results/stage05_primary", "figures/stage05_primary")
spec_file05 <- "metadata/outcome/stage05_primary_spec_v1.yml"
freeze_file05 <- "metadata/modeling/residual_freeze_v1.yml"
long_file05 <- "data/processed/residuals/neotrip_crossfitted_residuals_long_v1.tsv"
wide_file05 <- "data/processed/residuals/neotrip_crossfitted_residuals_wide_v1.tsv"
analysis_file05 <- "data/processed/outcome/neotrip_primary_analysis_v1.tsv"
init05 <- function() for (d in dirs05) dir.create(path(d), recursive = TRUE, showWarnings = FALSE)
stop05 <- function(message, decision = "STOP_REVIEW") {
  stop(structure(list(message = message, call = NULL, decision = decision),
    class = c("stage05_error", "error", "condition")))
}
verify_table05 <- function(file) {
  tab <- read_tsv(path(file))
  files <- if ("path" %in% names(tab)) tab$path else tab$local_relative_path
  if (any(!file.exists(path(files))) || !identical(unname(vapply(path(files), sha256, "")), tab$sha256))
    stop05("Prior-stage checksums differ; no regeneration permitted")
}
verify_stage04_05 <- function() {
  s <- jsonlite::read_json(path("results/stage04_residuals/stage04_status.json"), simplifyVector = TRUE)
  if (!identical(s$decision, "PROCEED") || !isTRUE(s$patients == 150) || !isTRUE(s$pathways == 6) ||
      !isTRUE(s$crossfitted_predictions == 900) || !isTRUE(s$patient_separation_failures == 0) ||
      !isTRUE(s$leakage_failures == 0) || !isTRUE(s$residual_qc_failures == 0) ||
      !isTRUE(s$pathways_beating_arm_only == 6) || !isTRUE(s$all_residuals_frozen))
    stop05("Stage 4 prerequisite failed")
  if (!file.exists(path(freeze_file05))) stop05("Residual freeze missing")
  freeze <- yaml::read_yaml(path(freeze_file05))
  mapping <- c(model_spec_sha256 = "metadata/modeling/stage04_model_spec_v1.yml",
    outer_folds_sha256 = "metadata/modeling/stage04_outer_folds_v1.tsv",
    inner_folds_sha256 = "metadata/modeling/stage04_inner_folds_v1.tsv",
    long_residual_file_sha256 = long_file05, wide_residual_file_sha256 = wide_file05,
    scaling_file_sha256 = "metadata/modeling/stage04_residual_scaling_v1.tsv")
  report <- do.call(rbind, lapply(names(mapping), function(key) {
    current <- if (file.exists(path(mapping[[key]]))) sha256(path(mapping[[key]])) else NA_character_
    data.frame(file = mapping[[key]], expected_sha256 = freeze[[key]], current_sha256 = current,
      matched = identical(freeze[[key]], current), record = freeze_file05)
  }))
  report <- rbind(report, data.frame(file = freeze_file05, expected_sha256 = s$residual_freeze_sha256,
    current_sha256 = sha256(path(freeze_file05)), matched = identical(s$residual_freeze_sha256, sha256(path(freeze_file05))),
    record = "results/stage04_residuals/stage04_status.json"))
  write_tsv(report, path("results/stage05_primary/residual_freeze_validation.tsv"))
  if (any(!report$matched) || !identical(freeze$outcome_information_used, FALSE))
    stop05("Residual freeze invalid; no regeneration permitted")
  for (f in c("metadata/composition/stage01_input_checksums.tsv",
    "results/stage02_composition/stage02_output_checksums.tsv",
    "results/stage03_pathways/stage03_output_checksums.tsv",
    "results/stage04_residuals/stage04_output_checksums.tsv")) verify_table05(f)
  invisible(TRUE)
}
specification05 <- function() {
  list(version = "v1", outcome = "pCR", outcome_positive = "pCR", outcome_negative = "RD",
    model_family = "binomial_logit", primary_effect = "residual_z", effect_scale = "one_sd_residual",
    tests = "two_sided", multiplicity_family = "six_confirmatory_pathways",
    multiplicity_method = "bonferroni", familywise_alpha = 0.05,
    composition_adjusted_model = TRUE,
    composition_coordinates = c("delta_cancer_epithelial", "delta_immune_fraction"),
    omitted_composition_component = "delta_stromal_fraction",
    treatment_interaction = "none", subgroup_analysis = "none",
    axes = axes05, directions = as.list(setNames(c("positive", rep("negative", 4), "positive"), axes05)),
    direction_definition = list(positive = "OR > 1", negative = "OR < 1"),
    arm_encoding = list(CT = 0L, "CT/A" = 1L),
    primary_terms = c("residual_z_axis", "arm", "baseline_z_axis",
      "delta_cancer_epithelial_10", "delta_immune_fraction_10"),
    base_terms = c("residual_z_axis", "arm", "baseline_z_axis"),
    inference = "Wald normal approximation; two-sided",
    effect_size_threshold_positive = 1.5, effect_size_threshold_negative = 1 / 1.5,
    composition_constraint_tolerance = 1e-8,
    composition_input_units = "percentage_points",
    composition_analysis_units = "fraction_change",
    composition_conversion = "Stage 2 percentage-point delta / 100; modelling fraction delta / 0.10",
    composition_units_provenance = "metadata/composition/reproduction_provenance.json",
    baseline_scaling = "Mean and sample SD across all 150 patients, without outcome values",
    attenuation_zero_tolerance = sqrt(.Machine$double.eps),
    stage06_eligibility = "DEFERRED")
}
freeze_spec05 <- function() {
  sidecar <- path("metadata/outcome/stage05_primary_spec_v1.sha256")
  if (file.exists(path(spec_file05)) || file.exists(sidecar)) return(read_spec05())
  spec <- specification05()
  spec$freeze_timestamp <- utc_now()
  yaml::write_yaml(spec, path(spec_file05))
  writeLines(sha256(path(spec_file05)), sidecar)
  read_spec05()
}
read_spec05 <- function() {
  file <- path(spec_file05)
  sidecar <- path("metadata/outcome/stage05_primary_spec_v1.sha256")
  if (!file.exists(file) || !file.exists(sidecar) || !identical(sha256(file), readLines(sidecar)))
    stop05("Primary specification checksum failed")
  spec <- yaml::read_yaml(file)
  canonical <- yaml::yaml.load(yaml::as.yaml(specification05()))
  if (!isTRUE(all.equal(spec[names(canonical)], canonical, check.attributes = FALSE)))
    stop05("Frozen primary specification differs")
  spec
}
read_numeric05 <- function(file, text) {
  x <- read_tsv(path(file))
  for (field in setdiff(names(x), text)) x[[field]] <- as.numeric(x[[field]])
  x
}
read_manifest05 <- function() {
  read_spec05()
  file <- path("metadata/neotrip_pair_manifest_v1.tsv")
  header <- names(read.delim(file, nrows = 0L, check.names = FALSE))
  fields <- c("patient_id", "arm", "pcr", "pcr_num", "paired")
  if (!all(fields %in% header)) stop05("Outcome source schema differs")
  x <- read.delim(file, sep = "\t", quote = "", comment.char = "", check.names = FALSE,
    colClasses = ifelse(header %in% fields, "character", "NULL"), na.strings = NULL)
  x[x$paired == "TRUE", fields, drop = FALSE]
}
prepare_analysis05 <- function(spec) {
  wide <- read_numeric05(wide_file05, c("patient_id", "arm"))
  long <- read_numeric05(long_file05, c("patient_id", "arm", "axis"))
  input <- read_numeric05("data/processed/pathways/neotrip_model_inputs_blinded_v1.tsv", c("patient_id", "arm"))
  if (nrow(wide) != 150L || anyDuplicated(wide$patient_id) || nrow(input) != 150L ||
      anyDuplicated(input$patient_id) || !setequal(input$patient_id, wide$patient_id) ||
      nrow(long) != 900L || anyDuplicated(paste(long$patient_id, long$axis)))
    stop05("Frozen patient identity invalid")
  input <- input[match(wide$patient_id, input$patient_id), ]
  if (!identical(input$arm, wide$arm)) stop05("Frozen treatment arm mismatch")
  result <- wide[, c("patient_id", "arm")]
  scaling <- do.call(rbind, lapply(spec$axes, function(a) {
    z <- input[[paste0("baseline_", a)]]
    data.frame(axis = a, mean = mean(z), sd = sd(z), n = length(z))
  }))
  write_tsv(scaling, path("metadata/outcome/stage05_baseline_scaling_v1.tsv"))
  for (a in spec$axes) {
    result[[paste0("baseline_", a)]] <- input[[paste0("baseline_", a)]]
    result[[paste0("residual_", a)]] <- wide[[paste0("residual_", a)]]
    result[[paste0("residual_z_", a)]] <- wide[[paste0("residual_z_", a)]]
    sc <- scaling[scaling$axis == a, ]
    result[[paste0("baseline_z_", a)]] <- (result[[paste0("baseline_", a)]] - sc$mean) / sc$sd
  }
  provenance <- jsonlite::read_json(path(spec$composition_units_provenance), simplifyVector = TRUE)
  if (!identical(provenance$units$change, "percentage_points"))
    stop05("Composition source units differ", "MODIFY_REVIEW")
  composition <- c("delta_cancer_epithelial", "delta_immune_fraction", "delta_stromal_fraction")
  for (field in composition) {
    result[[field]] <- input[[field]] / 100
    result[[paste0(field, "_10")]] <- result[[field]] / 0.10
  }
  write_tsv(data.frame(variable = composition, source_units = "percentage_points",
    analysis_units = "fraction_change", source_to_fraction_divisor = 100,
    modelling_divisor = 0.10, coefficient_units = "10 percentage points"),
    path("metadata/outcome/stage05_composition_units_v1.tsv"))
  manifest <- read_manifest05()
  labels <- c(pCR = 1L, RD = 0L)
  mapped <- unname(labels[manifest$pcr])
  i <- match(wide$patient_id, manifest$patient_id)
  qc <- c(unique_patients = as.integer(nrow(manifest) != 150L) + sum(duplicated(manifest$patient_id)),
    pcr_count = as.integer(sum(manifest$pcr == "pCR") != 63L),
    rd_count = as.integer(sum(manifest$pcr == "RD") != 87L),
    CT_count = as.integer(sum(manifest$arm == "CT") != 79L),
    CTA_count = as.integer(sum(manifest$arm == "CT/A") != 71L),
    missing_labels = sum(missing_value(manifest$pcr)),
    binary_coding = sum(is.na(mapped)) + sum(!manifest$pcr_num %in% c("0", "1")),
    label_numeric_agreement = sum(as.character(mapped) != manifest$pcr_num, na.rm = TRUE),
    patient_identity = sum(is.na(i)) + sum(!manifest$patient_id %in% wide$patient_id),
    treatment_identity = sum(manifest$arm[i] != wide$arm, na.rm = TRUE),
    residual_patient_axis = sum(table(factor(long$patient_id, levels = wide$patient_id),
      factor(long$axis, levels = spec$axes)) != 1L))
  write_tsv(data.frame(check = names(qc), failures = as.integer(qc), passed = qc == 0L),
    path("results/stage05_primary/outcome_qc.tsv"))
  if (anyNA(qc) || any(qc != 0L)) stop05("Outcome coding or patient identity failed")
  result$pcr <- mapped[i]
  result <- result[, c("patient_id", "arm", "pcr", setdiff(names(result), c("patient_id", "arm", "pcr")))]
  missing <- do.call(rbind, lapply(names(result), function(f) {
    bad <- if (is.numeric(result[[f]])) !is.finite(result[[f]]) else missing_value(result[[f]])
    data.frame(variable = f, patient_count = sum(bad), patient_ids = paste(result$patient_id[bad], collapse = ";"))
  }))
  write_tsv(missing, path("results/stage05_primary/missing_data_check.tsv"))
  write_tsv(result, path(analysis_file05))
  if (any(missing$patient_count > 0L)) stop05("Unexpected missing analysis values", "MODIFY_REVIEW")
  result
}
composition_check05 <- function(x, spec) {
  coordinates <- c("delta_cancer_epithelial", "delta_immune_fraction", "delta_stromal_fraction")
  total <- rowSums(x[, coordinates])
  report <- data.frame(patient_id = x$patient_id, x[, coordinates], delta_sum = total,
    absolute_error = abs(total), tolerance = spec$composition_constraint_tolerance,
    passed = is.finite(total) & abs(total) <= spec$composition_constraint_tolerance)
  write_tsv(report, path("results/stage05_primary/composition_constraint_check.tsv"))
  if (any(!report$passed)) stop05(paste(sum(!report$passed),
    "of 150 patients fail the required three-component zero-sum constraint after explicit percentage-point to fraction conversion;",
    "maximum absolute sum =", format(max(abs(total)), digits = 10),
    ". No logistic models were fitted; the two-coordinate parameterization requires methodology review."),
    "MODIFY_REVIEW")
  report
}
