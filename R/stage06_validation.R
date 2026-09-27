validate_crossfit06 <- function() {
  context <- blinded_context06()
  ss_input <- NULL
  reports <- list()
  model_files <- list.files(path("models/stage06"), pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
  model_files <- model_files[basename(model_files) %in% paste0(axes06, ".rds")]
  for (file in model_files) {
    label <- basename(dirname(file)); a <- sub("[.]rds$", "", basename(file))
    x <- context$x
    arm_specific <- label == "arm_specific"
    if (label == "S1_ssgsea") {
      if (is.null(ss_input)) {
        raw <- read06("data/raw/geo/GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz", blinded = TRUE)
        expression <- as.matrix(raw[, -(1:3)])
        rownames(expression) <- raw$HGNC_approved_symbol
        panel <- read06("metadata/pathways/confirmatory_gene_sets_v1.tsv", blinded = TRUE)
        genes <- setNames(lapply(axes06, function(axis) panel$gene_symbol[panel$axis == axis]), axes06)
        scores <- GSVA::gsva(GSVA::ssgseaParam(expression, genes, alpha = .25,
          normalize = TRUE, minSize = 1, maxSize = Inf), verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
        stored <- read06("data/processed/sensitivity/ssgsea_confirmatory_scores_v1.tsv.gz", blinded = TRUE)
        assert06(max(abs(t(scores[axes06, ]) - as.matrix(stored[, axes06]))) < 1e-12,
          "ssGSEA score reproduction differs")
        b <- match(context$pairs$baseline_sample_id, colnames(scores))
        d <- match(context$pairs$d1c2_sample_id, colnames(scores))
        ss_input <- context$x
        for (axis in axes06) {
          ss_input[[paste0("baseline_", axis)]] <- unname(scores[axis, b])
          ss_input[[paste0("delta_", axis)]] <- unname(scores[axis, d] - scores[axis, b])
        }
      }
      x <- ss_input
    }
    if (startsWith(label, "random_")) {
      key <- sub("^random_", "", label)
      checkpoint <- readRDS(path("data/processed/sensitivity/null_checkpoints", paste0(key, ".rds")))
      z <- checkpoint$residuals
      assert06(identical(z$patient_id, x$patient_id), "Null checkpoint patient mismatch")
      x[[paste0("baseline_", a)]] <- z$baseline_score
      x[[paste0("delta_", a)]] <- z$observed_delta
    } else {
      z <- read06(paste0("data/processed/sensitivity/residuals_", tolower(label), ".tsv"), blinded = TRUE)
      z <- z[z$axis == a, ]; z <- z[match(x$patient_id, z$patient_id), ]
    }
    saved <- readRDS(file)
    coverage <- integer(nrow(x))
    for (name in names(saved)) {
      fit <- saved[[name]]
      train <- match(fit$train_patient_ids, x$patient_id); test <- match(fit$test_patient_ids, x$patient_id)
      assert06(!anyNA(c(train, test)) && !length(intersect(train, test)) &&
        !length(intersect(fit$test_patient_ids, fit$inner_folds$patient_id)), "Saved model patient leakage")
      fold <- unique(context$outer$outer_fold[test])
      assert06(length(fold) == 1 && all(context$outer$outer_fold[train] != fold), "Saved outer fold differs")
      expected_train <- which(context$outer$outer_fold != fold &
        if (arm_specific) x$arm == x$arm[test[1]] else rep(TRUE, nrow(x)))
      assert06(identical(train, expected_train), "Saved training cohort differs")
      m <- x[, intersect(fit$predictors, names(x)), drop = FALSE]
      if ("arm" %in% names(m)) m$arm <- as.integer(x$arm == "CT/A")
      if ("arm_baseline_interaction" %in% fit$predictors)
        m$arm_baseline_interaction <- as.integer(x$arm == "CT/A") * x[[paste0("baseline_", a)]]
      m <- as.matrix(m[, fit$predictors, drop = FALSE])
      actual_hash <- digest::digest(list(m[train, , drop = FALSE], x[[paste0("delta_", a)]][train]), algo = "sha256")
      assert06(identical(actual_hash, fit$training_data_sha256), paste("Saved training values differ", label, a, name))
      co <- as.matrix(coef(fit$fit, s = fit$selected_rule))
      predicted <- as.numeric(cbind(1, m[test, , drop = FALSE]) %*% co)
      assert06(max(abs(predicted - z$expected_delta[test])) < 1e-10, "Stored OOF predictions differ")
      assert06(max(abs(z$observed_delta[test] - z$expected_delta[test] - z$residual[test])) < 1e-10, "Residual arithmetic differs")
      coverage[test] <- coverage[test] + 1L
    }
    assert06(all(coverage == 1), "OOF coverage differs")
    reports[[file]] <- data.frame(variant = label, axis = a, patients = nrow(z),
      folds = length(saved), separation_failures = 0L, prediction_failures = 0L, passed = TRUE)
  }
  report <- do.call(rbind, reports)
  assert06(nrow(report) == 630L, "Cross-fit model count differs")
  write06(report, "results/stage06_robustness/crossfit_reproduction_validation.tsv")
  report
}
validate_outputs06 <- function() {
  verify_prior06()
  verify_blind_freeze06()
  s <- read_spec06()
  five <- read06("results/stage06_robustness/five_sensitivity_stability.tsv")
  primary <- read06("results/stage05_primary/primary_residual_associations_v2.tsv")
  p <- primary[match(five$axis, primary$axis), ]
  assert06(nrow(five) == 30 && !anyDuplicated(paste(five$axis, five$sensitivity_id)), "Sensitivity count differs")
  expected <- sign(five$beta) == sign(p$beta) & five$beta >= log(p$ci95_lower) & five$beta <= log(p$ci95_upper)
  assert06(identical(expected, five$stable), "Frozen stability rule differs")
  null <- read06("results/stage06_robustness/permutation_nulls.tsv.gz")
  perm <- read06("results/stage06_robustness/permutation_summary.tsv")
  assert06(nrow(null) == 60000 && all(table(null$axis) == 10000) && all(null$passed), "Permutation count failed")
  indices <- readRDS(path("metadata/sensitivity/permutation_indices_v1.rds"))
  assert06(identical(dim(indices), c(150L, 10000L)) &&
    all(apply(indices, 2, function(i) identical(sort(i), 1:150))), "Invalid permutation labels")
  for (a in axes06) {
    z <- null$directional_statistic[null$axis == a]; row <- perm[perm$axis == a, ]
    assert06(abs(row$permutation_p - (1 + sum(z >= row$observed_T)) / 10001) < 1e-12 &&
      abs(row$permutation_95th_percentile - quantile(z, .95, names = FALSE, type = 7)) < 1e-12,
      "Permutation inference differs")
  }
  influence <- read06("results/stage06_robustness/influence_analysis.tsv")
  assert06(nrow(influence) == 900 && !anyDuplicated(paste(influence$axis, influence$omitted_patient_id)),
    "Influence coverage differs")
  expected_change <- abs(influence$beta_leave_one_out - influence$beta_full) / abs(influence$beta_full)
  assert06(max(abs(expected_change - influence$relative_effect_change), na.rm = TRUE) < 1e-10,
    "Influence arithmetic differs")
  genes <- read06("metadata/sensitivity/random_matched_gene_sets_v1.tsv.gz", blinded = TRUE)
  bins <- read06("metadata/sensitivity/random_gene_expression_bins_v1.tsv.gz", blinded = TRUE)
  panel <- read06("metadata/pathways/confirmatory_gene_sets_v1.tsv", blinded = TRUE)
  assert06(!any(genes$gene_symbol %in% panel$gene_symbol), "Confirmatory genes in null pool")
  for (a in axes06) {
    target <- panel$gene_symbol[panel$axis == a]
    target <- target[target %in% bins$gene_symbol]
    target_bins <- table(factor(bins$expression_bin[match(target, bins$gene_symbol)], levels = 1:20))
    for (j in 1:100) {
      z <- genes[genes$axis == a & genes$null_set == j, ]
      assert06(nrow(z) == length(target) && !anyDuplicated(z$gene_symbol) &&
        identical(table(factor(z$expression_bin, levels = 1:20)), target_bins), "Null expression matching differs")
    }
  }
  results <- read06("results/stage06_robustness/random_gene_null_results.tsv.gz")
  assert06(nrow(results) == 600 && all(table(results$axis) == 100) && all(is.finite(results$beta)), "Null result count differs")
  eligibility <- read06("results/stage06_robustness/stage06_eligibility.tsv")
  assert06(identical(eligibility$primary_multiplicity_gate, primary$multiplicity_gate) &&
    !any(eligibility$full_transportability_eligible & !primary$multiplicity_gate), "Primary multiplicity bypassed")
  checks <- c("prior_stage_files_unchanged", "blinded_freeze_intact", "five_variant_rule_verified",
    "permutation_assignments_valid", "permutation_pvalues_reproduced", "influence_refits_complete",
    "influence_arithmetic_verified", "null_genes_exclude_confirmatory", "null_sizes_and_bins_match",
    "six_hundred_null_results", "primary_multiplicity_preserved")
  write06(data.frame(check = checks, passed = TRUE), "results/stage06_robustness/stage06_validation.tsv")
  invisible(TRUE)
}
