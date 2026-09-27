ilr06 <- function(context) {
  full <- read06("data/processed/composition/neotrip_instaprism_full.tsv.gz", blinded = TRUE)
  fractions <- cbind(cancer_epithelial = full$Cancer_Epithelial, normal_epithelial = full$Normal_Epithelial,
    immune = rowSums(full[, c("B-cells", "T-cells", "Myeloid")]),
    stromal = rowSums(full[, c("Endothelial", "CAFs", "PVL")]))
  assert06(all(is.finite(fractions) & fractions > 0) && all(abs(rowSums(fractions) - 1) < 1e-8),
    "ILR requires positive complete composition", "MODIFY_REVIEW")
  basis <- contr.helmert(4)
  basis <- sweep(basis, 2, sqrt(colSums(basis^2)), "/")
  rownames(basis) <- colnames(fractions); colnames(basis) <- paste0("ilr", 1:3)
  assert06(max(abs(crossprod(basis) - diag(3))) < 1e-12 && max(abs(colSums(basis))) < 1e-12,
    "ILR basis not orthonormal")
  write06(data.frame(component = rownames(basis), basis), "metadata/sensitivity/stage06_ilr_basis_v1.tsv")
  scores <- log(fractions) %*% basis
  pairs <- read06("data/processed/composition/neotrip_pair_composition_change.tsv", blinded = TRUE)
  pairs <- pairs[match(context$x$patient_id, pairs$patient_id), ]
  b <- match(pairs$baseline_sample_id, full$sample_id); d <- match(pairs$d1c2_sample_id, full$sample_id)
  assert06(!anyNA(c(b, d)) && identical(full$patient_id[b], context$x$patient_id) &&
    identical(full$patient_id[d], context$x$patient_id), "ILR pair identity invalid")
  delta <- scores[d, , drop = FALSE] - scores[b, , drop = FALSE]
  colnames(delta) <- paste0("delta_ilr", 1:3)
  output <- data.frame(patient_id = context$x$patient_id, delta)
  write06(output, "data/processed/sensitivity/ilr_composition_changes_v1.tsv")
  write06(data.frame(samples = nrow(full), components = 4L, minimum_fraction = min(fractions),
    zero_fractions = sum(fractions == 0), pseudocount_used = FALSE,
    maximum_sum_error = max(abs(rowSums(fractions) - 1)), basis_verified = TRUE),
    "results/stage06_robustness/ilr_composition_qc.tsv")
  output
}
sensitivity_outcomes06 <- function(context) {
  x <- context$x
  effects <- list(); fits <- list()
  ids <- c(variants06[1:4], "arm_specific")
  for (id in ids) {
    z <- read06(paste0("data/processed/sensitivity/residuals_", tolower(id), ".tsv"), blinded = TRUE)
    for (a in axes06) {
      r <- z[z$axis == a, ]; r <- r[match(x$patient_id, r$patient_id), ]
      assert06(identical(r$patient_id, x$patient_id) && identical(r$arm, x$arm), "Sensitivity patient identity differs")
      frame <- outcome_frame06(x, a)
      frame$residual <- r$residual_z
      if (id == "S1_ssgsea") frame$baseline <- as.numeric(scale(r$baseline_score))
      fit <- fit_glm06(frame)
      key <- paste(id, a, sep = "__")
      fits[[key]] <- fit
      effects[[key]] <- data.frame(axis = a, sensitivity_id = id, glm_effect06(fit))
    }
  }
  ilr <- ilr06(context)
  comparisons <- list()
  for (a in axes06) {
    frame <- outcome_frame06(x, a)
    frame <- cbind(frame, ilr[, -1])
    forms <- list(A = pcr ~ residual + arm + baseline,
      B = pcr ~ residual + arm + baseline + cancer,
      C = pcr ~ residual + arm + baseline + cancer + immune + stromal,
      D = pcr ~ residual + arm + baseline + delta_ilr1 + delta_ilr2 + delta_ilr3)
    for (name in names(forms)) {
      fit <- fit_glm06(frame, forms[[name]])
      comparisons[[paste(a, name)]] <- data.frame(axis = a, model = name, glm_effect06(fit))
      fits[[paste("composition", name, a, sep = "__")]] <- fit
      if (name == "D") effects[[paste("S5", a)]] <- data.frame(axis = a,
        sensitivity_id = "S5_ilr_composition", glm_effect06(fit))
      if (name == "C") assert06(abs(coef(fit)[["residual"]] -
        context$primary$beta[context$primary$axis == a]) < 1e-10, "Primary model reproduction failed")
    }
  }
  all <- do.call(rbind, effects)
  five <- all[all$sensitivity_id %in% variants06, ]
  five <- five[order(match(five$axis, axes06), match(five$sensitivity_id, variants06)), ]
  write06(five, "results/stage06_robustness/five_primary_sensitivities.tsv")
  write06(all[all$sensitivity_id == "arm_specific", ],
    "results/stage06_robustness/arm_specific_residual_sensitivity.tsv")
  p <- context$primary[match(five$axis, context$primary$axis), ]
  five$same_direction <- sign(five$beta) == sign(p$beta)
  five$within_primary_95ci <- five$beta >= log(p$ci95_lower) & five$beta <= log(p$ci95_upper)
  five$stable <- five$same_direction & five$within_primary_95ci
  five$stable_count <- ave(as.integer(five$stable), five$axis, FUN = sum)
  five$sensitivity_4_of_5_gate <- five$stable_count >= 4
  write06(five, "results/stage06_robustness/five_sensitivity_stability.tsv")
  write06(do.call(rbind, comparisons), "results/stage06_robustness/composition_model_sensitivity.tsv")
  saveRDS(fits, path("models/stage06/outcome_sensitivities.rds"))
  diagnostics <- do.call(rbind, lapply(names(fits), function(name) {
    f <- fits[[name]]
    separation <- glm(formula(f), data = f$model, family = binomial("logit"),
      method = detectseparation::detect_separation)
    data.frame(model = name, n = nobs(f), rank = f$rank, parameters = length(coef(f)),
      converged = f$converged, finite = all(is.finite(coef(f))) && all(is.finite(vcov(f))),
      separation = any(is.infinite(coef(separation))))
  }))
  write06(diagnostics, "results/stage06_robustness/sensitivity_model_diagnostics.tsv")
  assert06(all(diagnostics$converged & diagnostics$finite & !diagnostics$separation),
    "Sensitivity model diagnostics failed", "MODIFY_REVIEW")
}
