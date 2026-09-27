influence06 <- function(context) {
  rows <- summaries <- list()
  for (a in axes06) {
    frame <- outcome_frame06(context$x, a)
    full <- fit_glm06(frame)
    beta <- unname(coef(full)[["residual"]])
    assert06(abs(beta - context$primary$beta[context$primary$axis == a]) < 1e-10, "Primary influence fit differs")
    loo <- vapply(seq_len(nrow(frame)), function(i) unname(coef(fit_glm06(frame[-i, ]))[["residual"]]), numeric(1))
    relative <- if (abs(beta) < 1e-8) rep(NA_real_, nrow(frame)) else abs(loo - beta) / abs(beta)
    z <- data.frame(axis = a, omitted_patient_id = context$x$patient_id, beta_full = beta,
      beta_leave_one_out = loo, absolute_beta_change = abs(loo - beta), relative_effect_change = relative,
      dfbeta_residual = dfbeta(full)[, "residual"], cooks_distance = cooks.distance(full), hat_value = hatvalues(full))
    rows[[a]] <- z
    summaries[[a]] <- data.frame(axis = a, refits = nrow(frame), beta_full = beta,
      max_relative_influence = if (all(is.na(relative))) NA_real_ else max(relative),
      influence_gate = if (all(is.na(relative))) FALSE else max(relative) <= .25)
  }
  write06(do.call(rbind, rows), "results/stage06_robustness/influence_analysis.tsv")
  write06(do.call(rbind, summaries), "results/stage06_robustness/influence_summary.tsv")
}
