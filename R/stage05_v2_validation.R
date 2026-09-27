validate_results_v2 <- function() {
  verify_v2()
  s <- read_spec_v2()
  x <- read_numeric05(analysis_v2, c("patient_id", "arm"))
  fits <- readRDS(path("results/stage05_primary/primary_models_v2.rds"))
  checks <- list()
  for (a in s$axes) for (name in c("base_model", s$primary_model_name)) {
    fit <- fits[[paste(a, name, sep = "__")]]
    terms <- if (name == "base_model") s$base_terms else s$primary_terms
    frame <- frame_v2(x, a, s)
    expected <- model.frame(reformulate(terms, response = "pcr"), data = frame)
    if (!isTRUE(all.equal(fit$model, expected, check.attributes = FALSE)))
      stop05("Saved model frame differs")
    if (fit$family$family != "binomial" || fit$family$link != "logit") stop05("Saved model family differs")
    m <- model.matrix(fit)
    probabilities <- plogis(as.numeric(m %*% coef(fit)))
    score_error <- max(abs(crossprod(m, x$pcr - probabilities)))
    independent_se <- sqrt(diag(solve(crossprod(m * sqrt(probabilities * (1 - probabilities))))))
    se_error <- max(abs(independent_se - sqrt(diag(vcov(fit)))))
    checks[[paste(a, name)]] <- data.frame(axis = a, model = name,
      frame_verified = TRUE, likelihood_score_max = score_error,
      information_se_max_error = se_error, passed = score_error < 1e-4 && se_error < 1e-4)
  }
  checks <- do.call(rbind, checks)
  write_tsv(checks, path("results/stage05_primary/model_reproduction_check_v2.tsv"))
  if (any(!checks$passed)) stop05("Model likelihood verification failed")
  inferred <- infer_v2(fits, s)
  for (item in list(list(file = "primary_residual_associations_v2.tsv", expected = inferred$primary),
                    list(file = "composition_adjustment_comparison_v2.tsv", expected = inferred$comparison))) {
    saved <- read_tsv(path("results/stage05_primary", item$file))
    if (!identical(names(saved), names(item$expected)) || nrow(saved) != nrow(item$expected))
      stop05("Inference output schema differs")
    for (f in names(saved)) {
      correct <- if (is.numeric(item$expected[[f]])) isTRUE(all.equal(as.numeric(saved[[f]]),
        item$expected[[f]], tolerance = 1e-11, check.attributes = FALSE)) else
        identical(saved[[f]], as.character(item$expected[[f]]))
      if (!correct) stop05("Inference reproduction differs")
    }
  }
  p <- inferred$primary
  if (!isTRUE(all.equal(p$p_bonferroni, p.adjust(p$p_value, method = "bonferroni", n = 6), tolerance = 1e-12)))
    stop05("Multiplicity validation failed")
  critical <- qnorm(1 - 0.05 / 12)
  if (any(abs(log(p$fwer_ci_lower) - (p$beta - critical * p$standard_error)) > 1e-12) ||
      any(abs(log(p$fwer_ci_upper) - (p$beta + critical * p$standard_error)) > 1e-12))
    stop05("Family-wise intervals differ")
  invisible(inferred)
}
