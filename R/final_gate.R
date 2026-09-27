dirs08 <- c("metadata/final_gate", "results/stage08_final_gate", "tables/stage08_final_gate")
read08 <- function(file) read_tsv(path(file))
json08 <- function(file) jsonlite::read_json(path(file), simplifyVector = TRUE)
write08 <- function(x, file) write_tsv(x, path(file))
assert08 <- function(ok, reason) if (!isTRUE(ok)) stop(reason, call. = FALSE)
verify08 <- function() {
  verify_prior07()
  manifest <- read08("results/stage07_prediction/stage07_output_checksums.tsv")
  assert08(identical(unname(vapply(path(manifest$local_relative_path), sha256, "")), manifest$sha256),
    "Stage 7 checksum mismatch")
  s4 <- json08("results/stage04_residuals/stage04_status.json")
  s5 <- json08("results/stage05_primary/stage05_status_v2.json")
  s6 <- json08("results/stage06_robustness/stage06_status.json")
  s7 <- json08("results/stage07_prediction/stage07_status.json")
  gate <- json08("results/stage06_robustness/external_transportability_gate.json")
  freeze <- read08("results/stage05_primary/residual_freeze_validation_v2.tsv")
  assert08(s4$all_residuals_frozen && s4$residual_freeze_sha256 ==
    sha256(path("metadata/modeling/residual_freeze_v1.yml")) &&
    all(freeze$matched == "TRUE") && identical(unname(vapply(path(freeze$file), sha256, "")),
      freeze$expected_sha256) && s5$residual_freeze_verified, "Residual freeze invalid")
  assert08(s5$decision == "PROCEED" && s6$decision == "PROCEED" && s7$decision == "PROCEED" &&
    s6$technical_failures == 0 && s7$leakage_failures == 0 && s7$reproduction_failures == 0 &&
    gate$eligible_count == 0 && gate$decision == "SKIP", "Stage 8 prerequisite failed")
  list(s4 = s4, s5 = s5, s6 = s6, s7 = s7)
}
apply_gate08 <- function() {
  s <- verify08()
  p <- read08("results/stage05_primary/primary_residual_associations_v2.tsv")
  e <- read08("results/stage06_robustness/stage06_eligibility.tsv")
  stability <- read08("results/stage06_robustness/five_sensitivity_stability.tsv")
  perm <- read08("results/stage06_robustness/permutation_summary.tsv")
  influence <- read08("results/stage06_robustness/influence_summary.tsv")
  random <- read08("results/stage06_robustness/random_gene_null_summary.tsv")
  for (tab in list(p, e, perm, influence, random))
    assert08(nrow(tab) == 6 && !anyDuplicated(tab$axis) && setequal(tab$axis, axes07), "Pathway identity mismatch")
  align <- function(tab) tab[match(p$axis, tab$axis), , drop = FALSE]
  e <- align(e); perm <- align(perm); influence <- align(influence); random <- align(random)
  assert08(identical(p$prespecified_direction, e$prespecified_direction) &&
    identical(p$odds_ratio, e$primary_OR) && identical(p$p_bonferroni, e$primary_bonferroni_p) &&
    identical(p$direction_matches, e$primary_direction_gate) &&
    identical(p$multiplicity_gate, e$primary_multiplicity_gate) &&
    identical(p$effect_size_gate, e$primary_effect_size_gate) &&
    identical(perm$permutation_p, e$permutation_p) &&
    identical(perm$permutation_gate, e$permutation_gate) &&
    identical(influence$influence_gate, e$influence_gate), "Frozen evidence sources disagree")
  for (a in p$axis) {
    z <- stability[stability$axis == a, ]
    assert08(nrow(z) == 5 && all(z$stable_count == e$sensitivity_stable_count[e$axis == a]),
      "Sensitivity summary mismatch")
  }
  assert08(all(p$multiplicity_gate == "FALSE") && all(e$full_transportability_eligible == "FALSE") &&
    all(e$classification %in% c("HYPOTHESIS_GENERATING", "NOT_ROBUST")) &&
    all(e$primary_direction_gate[e$classification == "HYPOTHESIS_GENERATING"] == "TRUE"),
    "Frozen confirmatory or direction gate contradiction")
  confirm <- data.frame(axis = p$axis, prespecified_direction = p$prespecified_direction,
    primary_OR = p$odds_ratio, primary_CI95 = paste(p$ci95_lower, p$ci95_upper, sep = " to "),
    primary_bonferroni_p = p$p_bonferroni, direction_matches = p$direction_matches,
    effect_size_gate = p$effect_size_gate, multiplicity_gate = p$multiplicity_gate)
  robust <- data.frame(axis = p$axis, sensitivity_stable_count = e$sensitivity_stable_count,
    permutation_gate = e$permutation_gate, influence_gate = e$influence_gate,
    random_gene_null_percentile = random$primary_percentile,
    full_transportability_eligible = e$full_transportability_eligible, stage06_classification = e$classification)
  write08(confirm, "tables/stage08_final_gate/confirmatory_summary.tsv")
  write08(robust, "tables/stage08_final_gate/robustness_summary.tsv")
  final <- confirm; names(final)[names(final) == "primary_bonferroni_p"] <- "bonferroni_p"
  final$sensitivity_stable_count <- e$sensitivity_stable_count
  final$permutation_p <- perm$permutation_p
  final$permutation_gate <- e$permutation_gate; final$influence_gate <- e$influence_gate
  final$random_gene_null_percentile <- random$primary_percentile
  final$full_transportability_eligible <- e$full_transportability_eligible
  final$final_evidence_level <- e$classification
  write08(final, "tables/stage08_final_gate/final_pathway_evidence.tsv")
  inventory <- list.files(root, recursive = TRUE, all.files = TRUE, include.dirs = FALSE)
  matches <- inventory[grepl("GSE122630", inventory, ignore.case = TRUE)]
  presence <- if (length(matches)) "PRESENT_NOT_ANALYSED" else "ABSENT_NOT_DOWNLOADED"
  write08(data.frame(status = presence, matching_path = if (length(matches)) matches else "",
    inspection = "Filename/path inventory only; no candidate file contents opened"),
    "results/stage08_final_gate/GSE122630_presence_audit.tsv")
  jsonlite::write_json(list(eligible_pathways = character(), eligible_count = 0L,
    GSE122630_status = "SKIP", reason = "No confirmatory pathway satisfied the complete prespecified eligibility rule.",
    data_presence = presence), path("results/stage08_final_gate/transportability_decision.json"),
    auto_unbox = TRUE, pretty = TRUE)
  prediction <- read08("results/stage07_prediction/prediction_summary.tsv")
  paired <- read08("results/stage07_prediction/paired_metric_summary.tsv")
  assert08(nrow(prediction) == 25 && nrow(paired) == 15, "Prediction summaries incomplete")
  increments <- lapply(c("P3 - P1", "P4 - P2"), function(comparison) {
    z <- paired[paired$comparison == comparison, ]
    z <- z[match(c("auc_difference", "brier_difference", "logloss_difference"), z$metric), ]
    assert08(nrow(z) == 3 && !anyNA(z), "Prediction comparison missing")
    data.frame(comparison = comparison, mean_auc_difference = z$mean_difference[1],
      mean_brier_difference = z$mean_difference[2], mean_logloss_difference = z$mean_difference[3],
      auc_repeats_favoring_residuals = z$repeats_favoring_later[1],
      brier_repeats_favoring_residuals = z$repeats_favoring_later[2],
      logloss_repeats_favoring_residuals = z$repeats_favoring_later[3])
  })
  increments <- do.call(rbind, increments)
  gain <- as.numeric(increments$mean_auc_difference) > 0 &
    as.numeric(increments$mean_brier_difference) < 0 & as.numeric(increments$mean_logloss_difference) < 0
  assert08(!any(gain), "Frozen prediction gain differs from expected protocol assessment")
  write08(increments, "tables/stage08_final_gate/prediction_increment_summary.tsv")
  jsonlite::write_json(list(prediction_secondary = TRUE, consistent_incremental_gain = any(gain),
    decision_curve_status = "NOT_RUN", interpretation = "No residual-addition comparison improves all three mean metrics simultaneously; predictive improvement and clinical utility are not supported."),
    path("results/stage08_final_gate/prediction_decision.json"), auto_unbox = TRUE, pretty = TRUE)
  composition <- json08("results/stage02_composition/stage02_status.json")
  stable <- any(e$primary_direction_gate == "TRUE" & e$sensitivity_4_of_5_gate == "TRUE")
  assert08(composition$stage01_verified && composition$instaprism_reproduced &&
    composition$composition_qc_failures == 0 && s$s4$patients == 150 && s$s5$patients == 150 &&
    s$s5$models_converged == 6 && s$s4$patient_separation_failures == 0 &&
    s$s4$leakage_failures == 0 && stable, "Fundamental continuation condition failed")
  assessment <- data.frame(
    criterion = c("public_primary_dataset", "paired_cohort_reconstructed", "composition_reproducible",
      "residual_construct_novel", "crossfitting_valid", "outcome_leakage_absent",
      "confirmatory_multiplicity_success", "robust_hypothesis_generating_signal",
      "incremental_prediction_success", "external_transportability_eligible",
      "external_transportability_completed", "fundamental_stop_condition"),
    status = c("PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "FAIL", "PASS", "FAIL", "FAIL", "SKIP", "PASS"),
    evidence_file = c(rep("results/stage02_composition/stage02_status.json", 4),
      "results/stage04_residuals/stage04_status.json", "results/stage07_prediction/leakage_audit.tsv",
      "results/stage05_primary/primary_residual_associations_v2.tsv",
      "results/stage06_robustness/stage06_eligibility.tsv", "results/stage07_prediction/paired_metric_summary.tsv",
      rep("results/stage06_robustness/external_transportability_gate.json", 2),
      "results/stage08_final_gate/stage08_status.json"),
    reason = c("Public NeoTRIP inputs verified by the frozen reconstruction audit.",
      "Exactly 150 paired analytic patients reconstructed.", "Published InstaPrism composition reproduced; QC passed.",
      "No substantial residual-method overlap identified in inspected code and abstract; novelty remains qualified because publisher full text/supplements were unavailable.",
      "Frozen cross-fitted residuals and separation checks valid.", "Frozen biological and prediction leakage checks report zero failures.",
      "Zero pathways pass the frozen Stage 5 multiplicity gate.", "At least one direction-matching association is stable in at least four of five sensitivities.",
      "Neither residual addition improves all three mean predictive metrics.", "Zero pathways satisfy the complete eligibility rule, including mandatory multiplicity.",
      "External transportability skipped by the frozen gate.", "PASS denotes absence of a fundamental STOP condition; all continuation prerequisites hold."))
  write08(assessment, "results/stage08_final_gate/protocol_assessment.tsv")
  list(stage = 8L, timestamp = utc_now(), stage05_verified = TRUE, stage06_verified = TRUE,
    stage07_verified = TRUE, patients = 150L, confirmatory_pathways = 6L,
    multiplicity_confirmed_pathways = 0L, full_transportability_eligible_pathways = 0L,
    GSE122630_status = "SKIP", consistent_incremental_prediction_gain = FALSE,
    fundamental_stop_condition = FALSE, decision = "PROCEED_NEOTRIP_ONLY",
    decision_reason = paste("Technical reconstruction, composition, residuals and primary models are valid;",
      "directionally stable hypothesis-generating associations remain. No multiplicity-confirmed pathway or",
      "external-transportability eligibility exists. Continue with single-cohort scope; no positive confirmatory",
      "result, predictive improvement, clinical utility or high-impact external-validation framing is supported."))
}
