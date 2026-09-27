fit_logistic05 <- function(data, spec) {
  models <- list()
  diagnostics <- list()
  for (axis in spec$axes) {
    frame <- data.frame(pcr = data$pcr, residual_z_axis = data[[paste0("residual_z_", axis)]],
      arm = unname(unlist(spec$arm_encoding)[data$arm]),
      baseline_z_axis = data[[paste0("baseline_z_", axis)]],
      delta_cancer_epithelial_10 = data$delta_cancer_epithelial_10,
      delta_immune_fraction_10 = data$delta_immune_fraction_10)
    for (model_name in c("base_model", "composition_adjusted_model")) {
      terms <- if (model_name == "base_model") spec$base_terms else spec$primary_terms
      formula <- reformulate(terms, response = "pcr")
      warnings <- character()
      fit <- withCallingHandlers(glm(formula, data = frame, family = binomial(link = "logit"),
        na.action = na.fail), warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      se <- sqrt(diag(vcov(fit)))
      separation_warning <- any(grepl("separation|did not converge|numerically 0 or 1", warnings, ignore.case = TRUE))
      unstable <- isTRUE(fit$boundary) || any(!is.finite(coef(fit))) ||
        any(!is.finite(se)) || any(se <= 0)
      row <- data.frame(axis = axis, model = model_name, converged = fit$converged,
        iterations = fit$iter, finite_coefficients = all(is.finite(coef(fit))),
        finite_standard_errors = all(is.finite(se)), rank = fit$rank,
        number_parameters = ncol(model.matrix(fit)), number_observations = nobs(fit),
        number_events = sum(frame$pcr), minimum_fitted_probability = min(fitted(fit)),
        maximum_fitted_probability = max(fitted(fit)), warnings = paste(unique(warnings), collapse = ";"),
        separation_warning = separation_warning, numerical_instability = unstable)
      row$passed <- row$converged && row$finite_coefficients && row$finite_standard_errors &&
        row$rank == row$number_parameters && row$number_observations == 150 &&
        !row$separation_warning && !row$numerical_instability
      key <- paste(axis, model_name, sep = "__")
      models[[key]] <- fit
      diagnostics[[key]] <- row
    }
  }
  diagnostics <- do.call(rbind, diagnostics)
  write_tsv(diagnostics, path("results/stage05_primary/model_diagnostics.tsv"))
  saveRDS(models, path("results/stage05_primary/primary_models_v1.rds"))
  if (any(!diagnostics$passed)) stop05("Logistic model diagnostics require review", "MODIFY_REVIEW")
  list(models = models, diagnostics = diagnostics)
}
