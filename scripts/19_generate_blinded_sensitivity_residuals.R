source("D:/TNBC-ResponseResiduals/R/paths.R")
setwd(root)
source(path("R/checksums.R"))
source(path("R/stage06_common.R"))
source(path("R/sensitivity_residuals.R"))
source(path("R/ssgsea_scoring.R"))
source(path("R/treatment_sensitivity.R"))
source(path("R/random_gene_nulls.R"))
run06("19_generate_blinded_sensitivity_residuals", {
  verify_prior06()
  read_spec06()
  if (file.exists(path(freeze06_file))) {
    verify_blind_freeze06()
  } else {
    context <- blinded_context06()
    context$arm_inner <- arm_inner06(context)
    m <- expression06()
    ss <- ssgsea06(m, context)
    ss_pairs <- pair_scores06(ss, context)
    write06(ss_pairs, "data/processed/sensitivity/ssgsea_blinded_pairs_v1.tsv")
    construct_variant06(ss_pairs, context, "S1_ssgsea")
    construct_variant06(context$x, context, "S2_no_baseline_composition",
      predictors = c("arm", paste0("baseline_", axes06)))
    construct_variant06(context$x, context, "S3_lambda_min", rule = "lambda.min")
    construct_variant06(context$x, context, "S4_arm_baseline_interaction", interaction = TRUE)
    construct_variant06(context$x, context, "arm_specific",
      predictors = setdiff(base_predictors06, "arm"), arm_specific = TRUE)
    random_residuals06(m, context)
    residual_files <- c(paste0("data/processed/sensitivity/residuals_", tolower(variants06[1:4]), ".tsv"),
      "data/processed/sensitivity/residuals_arm_specific.tsv",
      "data/processed/sensitivity/random_null_residuals_v1.tsv.gz")
    qc <- do.call(rbind, lapply(residual_files, function(f) {
      z <- read06(f, blinded = TRUE)
      group <- if ("null_set" %in% names(z)) paste(z$axis, z$null_set) else z$axis
      counts <- table(group)
      data.frame(file = f, patients = length(unique(z$patient_id)), predictions = nrow(z),
        missing = sum(is.na(z$residual_z)), nonfinite = sum(!is.finite(z$residual_z)),
        duplicate_keys = sum(duplicated(paste(z$patient_id, group))), forbidden_fields = length(forbidden06(names(z))),
        passed = all(counts == 150) && all(is.finite(z$residual_z)) && !anyDuplicated(paste(z$patient_id, group)))
    }))
    write06(qc, "results/stage06_robustness/blinded_residual_qc.tsv")
    assert06(all(qc$passed), "Blinded residual QC failed")
    frozen <- c(residual_files, "data/processed/sensitivity/random_null_residuals_raw_v1.tsv.gz",
      "metadata/sensitivity/random_matched_gene_sets_v1.tsv.gz",
      "data/processed/sensitivity/ssgsea_confirmatory_scores_v1.tsv.gz",
      "metadata/sensitivity/arm_specific_inner_folds_v1.tsv")
    freeze_yaml06(list(timestamp = utc_now(), clinical_labels_accessed = FALSE,
      specification_sha256 = sha256(path(spec06_file)), files = hashes06(frozen)), freeze06_file)
  }
})
cat("Blinded residuals frozen\n")
