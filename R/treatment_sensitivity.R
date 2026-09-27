arm_inner06 <- function(context) {
  output <- list()
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  for (a in c("CT", "CT/A")) for (k in 1:10) {
    x <- context$outer[context$outer$arm == a & context$outer$outer_fold != k, ]
    x <- x[order(x$patient_id, method = "radix"), ]
    set.seed(319641 + 100 * match(a, c("CT", "CT/A")) + k)
    output[[paste(a, k)]] <- data.frame(outer_fold = k, patient_id = x$patient_id,
      arm = a, inner_fold = sample(rep(1:5, length.out = nrow(x))))
  }
  out <- do.call(rbind, output)
  write06(out, "metadata/sensitivity/arm_specific_inner_folds_v1.tsv")
  out
}
treatment_models06 <- function(context) {
  x <- context$x
  result <- list()
  for (a in axes06) {
    f <- outcome_frame06(x, a)
    fit <- fit_glm06(f, pcr ~ residual * arm + baseline + cancer + immune + stromal)
    b <- coef(fit); v <- vcov(fit)
    result[[a]] <- data.frame(axis = a, residual_main_beta = b[["residual"]],
      interaction_beta = b[["residual:arm"]], interaction_se = sqrt(v["residual:arm", "residual:arm"]),
      interaction_ci95_lower = exp(b[["residual:arm"]] - qnorm(.975) * sqrt(v["residual:arm", "residual:arm"])),
      interaction_ci95_upper = exp(b[["residual:arm"]] + qnorm(.975) * sqrt(v["residual:arm", "residual:arm"])),
      interaction_p = 2 * pnorm(abs(b[["residual:arm"]] / sqrt(v["residual:arm", "residual:arm"])), lower.tail = FALSE),
      residual_OR_CT = exp(b[["residual"]]), residual_OR_CTA = exp(b[["residual"]] + b[["residual:arm"]]))
  }
  write06(do.call(rbind, result), "results/stage06_robustness/outcome_arm_interaction.tsv")
}
