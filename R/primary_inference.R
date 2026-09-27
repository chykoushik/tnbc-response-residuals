effect05 <- function(fit) {
  beta <- unname(coef(fit)[["residual_z_axis"]])
  se <- sqrt(vcov(fit)["residual_z_axis", "residual_z_axis"])
  z <- beta / se
  data.frame(beta = beta, standard_error = se, odds_ratio = exp(beta),
    ci95_lower = exp(beta - qnorm(0.975) * se), ci95_upper = exp(beta + qnorm(0.975) * se),
    wald_z = z, p_value = 2 * pnorm(abs(z), lower.tail = FALSE))
}
inference05 <- function(fits, spec) {
  primary <- comparison <- list()
  critical <- qnorm(1 - (spec$familywise_alpha / length(spec$axes)) / 2)
  for (axis in spec$axes) {
    adjusted <- effect05(fits[[paste(axis, "composition_adjusted_model", sep = "__")]])
    base <- effect05(fits[[paste(axis, "base_model", sep = "__")]])
    direction <- spec$directions[[axis]]
    row <- data.frame(axis = axis, prespecified_direction = direction, n = 150L, pcr_n = 63L, rd_n = 87L,
      adjusted, p_bonferroni = min(1, adjusted$p_value * length(spec$axes)),
      fwer_ci_lower = exp(adjusted$beta - critical * adjusted$standard_error),
      fwer_ci_upper = exp(adjusted$beta + critical * adjusted$standard_error))
    row$direction_matches <- if (direction == "positive") row$odds_ratio > 1 else row$odds_ratio < 1
    row$effect_size_gate <- if (direction == "positive") row$odds_ratio >= spec$effect_size_threshold_positive else
      row$odds_ratio <= spec$effect_size_threshold_negative
    row$multiplicity_gate <- row$p_bonferroni < spec$familywise_alpha &&
      (row$fwer_ci_lower > 1 || row$fwer_ci_upper < 1)
    primary[[axis]] <- row
    attenuation <- if (abs(base$beta) <= spec$attenuation_zero_tolerance) NA_real_ else
      1 - abs(adjusted$beta) / abs(base$beta)
    comparison[[axis]] <- data.frame(axis = axis, model = c("base_model", "composition_adjusted_model"),
      rbind(base, adjusted), beta_change = adjusted$beta - base$beta, relative_attenuation = attenuation)
  }
  primary <- do.call(rbind, primary)
  comparison <- do.call(rbind, comparison)
  list(primary = primary, comparison = comparison)
}
save_inference05 <- function(inference) {
  p <- inference$primary
  write_tsv(p, path("results/stage05_primary/primary_residual_associations.tsv"))
  write_tsv(inference$comparison, path("results/stage05_primary/composition_adjustment_comparison.tsv"))
  table <- data.frame(axis = p$axis, prespecified_direction = p$prespecified_direction,
    OR_per_SD = p$odds_ratio, CI95 = sprintf("%.4f to %.4f", p$ci95_lower, p$ci95_upper),
    Bonferroni_p = p$p_bonferroni, FWER_CI = sprintf("%.4f to %.4f", p$fwer_ci_lower, p$fwer_ci_upper),
    direction_matches = p$direction_matches, effect_size_gate = p$effect_size_gate,
    multiplicity_gate = p$multiplicity_gate)
  write_tsv(table, path("tables/stage05_primary_residual_associations.tsv"))
}
