spec_v2 <- "metadata/outcome/stage05_primary_spec_v2.yml"
amend_v2 <- "metadata/outcome/stage05_methodology_amendment_v2.yml"
analysis_v2 <- "data/processed/outcome/neotrip_primary_analysis_v2.tsv"
verify_v2 <- function() {
  for (f in c("metadata/composition/stage01_input_checksums.tsv",
    "results/stage02_composition/stage02_output_checksums.tsv",
    "results/stage03_pathways/stage03_output_checksums.tsv",
    "results/stage04_residuals/stage04_output_checksums.tsv",
    "results/stage05_primary/stage05_output_checksums.tsv")) verify_table05(f)
  previous <- jsonlite::read_json(path("results/stage05_primary/stage05_status.json"), simplifyVector = TRUE)
  if (!identical(previous$decision, "MODIFY_REVIEW") || previous$primary_models != 0L)
    stop05("Prior model-free amendment prerequisite failed")
  state <- jsonlite::read_json(path("metadata/outcome/stage05_execution_v1.json"), simplifyVector = TRUE)
  if (!identical(state$models_attempted, FALSE)) stop05("Prior outcome models were attempted")
  status <- jsonlite::read_json(path("results/stage04_residuals/stage04_status.json"), simplifyVector = TRUE)
  if (status$decision != "PROCEED" || status$patients != 150L || status$pathways != 6L ||
      status$crossfitted_predictions != 900L || status$patient_separation_failures != 0L ||
      status$leakage_failures != 0L || status$residual_qc_failures != 0L ||
      status$pathways_beating_arm_only != 6L || !isTRUE(status$all_residuals_frozen)) stop05("Stage 4 prerequisite failed")
  freeze <- yaml::read_yaml(path(freeze_file05))
  if (!identical(freeze$outcome_information_used, FALSE)) stop05("Residual construction leakage")
  report <- read_tsv(path("results/stage05_primary/residual_freeze_validation.tsv"))
  report$current_sha256 <- vapply(path(report$file), sha256, "")
  report$matched <- report$current_sha256 == report$expected_sha256
  write_tsv(report, path("results/stage05_primary/residual_freeze_validation_v2.tsv"))
  if (any(!report$matched) || sha256(path(freeze_file05)) != status$residual_freeze_sha256)
    stop05("Frozen residual checksum mismatch")
  invisible(TRUE)
}
freeze_yaml_v2 <- function(value, file) {
  checksum <- sub("[.]yml$", ".sha256", file)
  if (file.exists(path(file)) || file.exists(path(checksum))) {
    if (!file.exists(path(file)) || !file.exists(path(checksum)) ||
        !identical(readLines(path(checksum)), sha256(path(file)))) stop05("Frozen v2 file mismatch")
    old <- yaml::read_yaml(path(file))
    canonical <- yaml::yaml.load(yaml::as.yaml(value))
    fields <- setdiff(names(canonical), "freeze_timestamp")
    if (!isTRUE(all.equal(old[fields], canonical[fields], check.attributes = FALSE)))
      stop05("Frozen v2 specification differs")
    return(old)
  }
  yaml::write_yaml(value, path(file))
  writeLines(sha256(path(file)), path(checksum))
  yaml::read_yaml(path(file))
}
amendment_v2 <- function() {
  list(amendment_version = "v2", freeze_timestamp = utc_now(), previous_status = "MODIFY_REVIEW",
    previous_models_fitted = 0L, reason = "four_component_composition",
    previous_assumption = "cancer + immune + stromal sums to one",
    verified_definition = "cancer + normal_epithelial + immune + stromal sums to one",
    primary_reference_component = "normal_epithelial",
    primary_coordinates = c("delta_cancer_epithelial", "delta_immune_fraction", "delta_stromal_fraction"),
    outcome_results_seen_before_amendment = FALSE, residual_construction_changed = FALSE,
    pathway_panel_changed = FALSE, outcome_definition_changed = FALSE,
    basis = "User-authorized correction of verified composition definition before association modelling")
}
specification_v2 <- function() {
  s <- read_spec05()
  s$version <- "v2"
  s$freeze_timestamp <- utc_now()
  s$composition_definition <- list(components = c("cancer_epithelial", "normal_epithelial", "immune", "stromal"))
  s$composition_coordinates <- c("delta_cancer_epithelial", "delta_immune_fraction", "delta_stromal_fraction")
  s$omitted_composition_component <- NULL
  s$omitted_reference_component <- "delta_normal_epithelial"
  s$composition_scale <- "10_percentage_point_absolute_change"
  s$primary_terms <- c(s$primary_terms, "delta_stromal_fraction_10")
  s$primary_model_name <- "composition_adjusted_model_v2"
  s$sample_sum_tolerance <- 1e-8
  s$separation_diagnostic <- "detectseparation::detect_separation; diagnostic only, no estimator substitution"
  s$separation_package_version <- as.character(packageVersion("detectseparation"))
  s$separation_documentation <- "https://ikosmidis.r-universe.dev/detectseparation/doc/manual.html"
  s$methodology_amendment_sha256 <- sha256(path(amend_v2))
  s
}
read_spec_v2 <- function() {
  for (file in c(amend_v2, spec_v2)) {
    sidecar <- sub("[.]yml$", ".sha256", file)
    if (!file.exists(path(file)) || !file.exists(path(sidecar)) ||
        !identical(readLines(path(sidecar)), sha256(path(file)))) stop05("V2 freeze missing or changed")
  }
  s <- yaml::read_yaml(path(spec_v2))
  expected <- specification_v2()
  fields <- setdiff(names(expected), "freeze_timestamp")
  if (!isTRUE(all.equal(s[fields], expected[fields], check.attributes = FALSE)))
    stop05("V2 specification differs")
  s
}
prepare_v2 <- function(s) {
  x <- read_numeric05(analysis_file05, c("patient_id", "arm"))
  wide <- read_numeric05(wide_file05, c("patient_id", "arm"))
  manifest <- read_manifest05()
  i <- match(x$patient_id, manifest$patient_id)
  labels <- c(pCR = 1L, RD = 0L)
  mapped <- unname(labels[manifest$pcr])
  identity <- sum(is.na(i)) + sum(duplicated(x$patient_id)) +
    as.integer(!identical(x$patient_id, wide$patient_id)) +
    sum(manifest$arm[i] != x$arm, na.rm = TRUE) +
    sum(x$pcr != mapped[i], na.rm = TRUE) +
    sum(manifest$pcr_num != as.character(mapped), na.rm = TRUE)
  if (anyNA(mapped) || any(!manifest$pcr_num %in% c("0", "1")) || identity != 0L)
    stop05("Outcome identity or coding mismatch")
  for (a in s$axes) for (field in paste0(c("residual_", "residual_z_"), a)) {
    if (!isTRUE(all.equal(x[[field]], wide[[field]], tolerance = 1e-12))) stop05("Residual values changed")
  }
  full <- read_numeric05("data/processed/composition/neotrip_instaprism_full.tsv.gz",
    c("sample_id", "gsm", "patient_id", "timepoint"))
  pairs <- read_tsv(path("data/processed/composition/neotrip_pair_composition_change.tsv"))
  sample <- full[, c("sample_id", "gsm", "patient_id", "timepoint")]
  sample$cancer_epithelial <- full$Cancer_Epithelial
  sample$normal_epithelial <- full$Normal_Epithelial
  sample$immune_fraction <- rowSums(full[, c("B-cells", "T-cells", "Myeloid")])
  sample$stromal_fraction <- rowSums(full[, c("Endothelial", "CAFs", "PVL")])
  component_names <- c("cancer_epithelial", "normal_epithelial", "immune_fraction", "stromal_fraction")
  sample$composition_sum <- rowSums(sample[, component_names])
  sample$absolute_error <- abs(sample$composition_sum - 1)
  sample$tolerance <- s$sample_sum_tolerance
  sample$passed <- sample$absolute_error <= sample$tolerance &
    apply(sample[, component_names], 1L, function(z) all(is.finite(z) & z >= 0 & z <= 1))
  write_tsv(sample, path("results/stage05_primary/composition_sample_sum_check_v2.tsv"))
  if (nrow(sample) != 401L || anyDuplicated(sample$sample_id) || any(!sample$passed))
    stop05("Sample composition reconstruction failed", "MODIFY_REVIEW")
  j <- match(x$patient_id, pairs$patient_id)
  b <- match(pairs$baseline_sample_id[j], sample$sample_id)
  d <- match(pairs$d1c2_sample_id[j], sample$sample_id)
  if (anyNA(c(j, b, d)) || anyDuplicated(pairs$patient_id) ||
      any(sample$patient_id[b] != x$patient_id | sample$patient_id[d] != x$patient_id) ||
      any(sample$timepoint[b] != "Baseline" | sample$timepoint[d] != "D1C2"))
    stop05("Composition pair identity failed")
  for (a in component_names) {
    delta <- sample[[a]][d] - sample[[a]][b]
    field <- paste0("delta_", a)
    if (field %in% names(x) && !isTRUE(all.equal(x[[field]], delta, tolerance = 1e-10)))
      stop05("Composition delta reconstruction differs", "MODIFY_REVIEW")
    if (!field %in% names(x)) x[[field]] <- delta
    x[[paste0(field, "_10")]] <- x[[field]] / 0.10
  }
  fields <- paste0("delta_", component_names)
  constraint <- data.frame(patient_id = x$patient_id, x[, fields], delta_sum = rowSums(x[, fields]))
  constraint$absolute_error <- abs(constraint$delta_sum)
  constraint$tolerance <- s$composition_constraint_tolerance
  constraint$passed <- constraint$absolute_error <= constraint$tolerance
  write_tsv(constraint, path("results/stage05_primary/composition_constraint_check_v2.tsv"))
  if (any(!constraint$passed)) stop05("Four-component paired constraint failed", "MODIFY_REVIEW")
  missing <- do.call(rbind, lapply(names(x), function(field) {
    bad <- if (is.numeric(x[[field]])) !is.finite(x[[field]]) else missing_value(x[[field]])
    data.frame(variable = field, patient_count = sum(bad), patient_ids = paste(x$patient_id[bad], collapse = ";"))
  }))
  write_tsv(missing, path("results/stage05_primary/missing_data_check_v2.tsv"))
  counts <- c(patients = nrow(x), pCR = sum(x$pcr == 1L), RD = sum(x$pcr == 0L),
    CT = sum(x$arm == "CT"), CTA = sum(x$arm == "CT/A"), identity_failures = identity,
    missing_values = sum(missing$patient_count))
  target <- c(150L, 63L, 87L, 79L, 71L, 0L, 0L)
  qc <- data.frame(check = names(counts), observed = as.integer(counts), expected = target,
    passed = counts == target)
  write_tsv(qc, path("results/stage05_primary/outcome_qc_v2.tsv"))
  write_tsv(x, path(analysis_v2))
  if (any(!qc$passed)) stop05("Outcome QC failed",
    if (sum(missing$patient_count) > 0) "MODIFY_REVIEW" else "STOP_REVIEW")
  x
}
frame_v2 <- function(x, a, s) {
  data.frame(pcr = x$pcr, residual_z_axis = x[[paste0("residual_z_", a)]],
    arm = unname(unlist(s$arm_encoding)[x$arm]), baseline_z_axis = x[[paste0("baseline_z_", a)]],
    delta_cancer_epithelial_10 = x$delta_cancer_epithelial_10,
    delta_immune_fraction_10 = x$delta_immune_fraction_10,
    delta_stromal_fraction_10 = x$delta_stromal_fraction_10)
}
design_check_v2 <- function(m, axis, model) {
  finite <- all(is.finite(m))
  rank <- if (finite) qr(m)$rank else NA_integer_
  zero <- colnames(m)[vapply(seq_len(ncol(m)), function(j) length(unique(m[, j])) <= 1L, logical(1))]
  duplicates <- sum(duplicated(t(m)))
  data.frame(axis = axis, model = model, n_rows = nrow(m), n_columns = ncol(m), matrix_rank = rank,
    zero_variance_columns = paste(zero, collapse = ";"),
    zero_variance_predictors = sum(zero != "(Intercept)"), duplicate_columns = duplicates,
    nonfinite_values = sum(!is.finite(m)),
    passed = finite && nrow(m) == 150L && rank == ncol(m) &&
      sum(zero != "(Intercept)") == 0L && duplicates == 0L)
}
fit_v2 <- function(x, s) {
  designs <- list()
  for (a in s$axes) for (name in c("base_model", s$primary_model_name)) {
    terms <- if (name == "base_model") s$base_terms else s$primary_terms
    frame <- frame_v2(x, a, s)
    m <- model.matrix(reformulate(terms), data = frame)
    designs[[paste(a, name)]] <- design_check_v2(m, a, name)
  }
  design <- do.call(rbind, designs)
  write_tsv(design, path("results/stage05_primary/design_matrix_check_v2.tsv"))
  if (any(!design$passed)) stop05("Design matrix validation failed", "MODIFY_REVIEW")
  models <- diagnostics <- list()
  for (a in s$axes) for (name in c("base_model", s$primary_model_name)) {
    terms <- if (name == "base_model") s$base_terms else s$primary_terms
    formula <- reformulate(terms, response = "pcr")
    frame <- frame_v2(x, a, s)
    warnings <- character()
    detection <- withCallingHandlers(glm(formula, data = frame, family = binomial("logit"),
      method = detectseparation::detect_separation, na.action = na.fail),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    separation <- any(is.infinite(coef(detection)))
    fit <- withCallingHandlers(glm(formula, data = frame, family = binomial("logit"), na.action = na.fail),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    se <- sqrt(diag(vcov(fit)))
    m <- model.matrix(fit)
    weighted <- m * sqrt(fit$weights)
    condition <- if (all(is.finite(weighted))) rcond(crossprod(weighted)) else 0
    unstable <- isTRUE(fit$boundary) || !is.finite(condition) || condition <= .Machine$double.eps ||
      any(grepl("did not converge|numerically 0 or 1", warnings, ignore.case = TRUE))
    row <- data.frame(axis = a, model = name, converged = fit$converged, iterations = fit$iter,
      finite_coefficients = all(is.finite(coef(fit))), finite_standard_errors = all(is.finite(se)),
      rank = fit$rank, number_parameters = ncol(m), number_observations = nobs(fit),
      pcr_n = sum(frame$pcr), rd_n = sum(frame$pcr == 0),
      minimum_fitted_probability = min(fitted(fit)), maximum_fitted_probability = max(fitted(fit)),
      warnings = paste(unique(warnings), collapse = ";"), separation = separation,
      numerical_instability = unstable, information_rcond = condition)
    row$passed <- row$converged && row$finite_coefficients && row$finite_standard_errors &&
      row$rank == row$number_parameters && !separation && !unstable
    key <- paste(a, name, sep = "__")
    models[[key]] <- fit
    diagnostics[[key]] <- row
  }
  diagnostics <- do.call(rbind, diagnostics)
  write_tsv(diagnostics, path("results/stage05_primary/model_diagnostics_v2.tsv"))
  saveRDS(models, path("results/stage05_primary/primary_models_v2.rds"))
  if (any(!diagnostics$passed)) stop05("Logistic model diagnostics failed", "MODIFY_REVIEW")
  models
}
infer_v2 <- function(models, s) {
  adapted <- models
  for (a in s$axes) adapted[[paste(a, "composition_adjusted_model", sep = "__")]] <-
    models[[paste(a, s$primary_model_name, sep = "__")]]
  out <- inference05(adapted, s)
  out$comparison$model[out$comparison$model == "composition_adjusted_model"] <- s$primary_model_name
  out
}
save_inference_v2 <- function(out) {
  p <- out$primary
  write_tsv(p, path("results/stage05_primary/primary_residual_associations_v2.tsv"))
  write_tsv(out$comparison, path("results/stage05_primary/composition_adjustment_comparison_v2.tsv"))
  write_tsv(data.frame(axis = p$axis, prespecified_direction = p$prespecified_direction,
    OR_per_SD = p$odds_ratio, CI95 = sprintf("%.4f to %.4f", p$ci95_lower, p$ci95_upper),
    Bonferroni_p = p$p_bonferroni, FWER_CI = sprintf("%.4f to %.4f", p$fwer_ci_lower, p$fwer_ci_upper),
    direction_matches = p$direction_matches, effect_size_gate = p$effect_size_gate,
    multiplicity_gate = p$multiplicity_gate), path("tables/stage05_primary_residual_associations_v2.tsv"))
}
figures_v2 <- function(out, s) {
  for (comparison in c(FALSE, TRUE)) {
    name <- if (comparison) "composition_adjustment_comparison_v2" else "primary_residual_forest_v2"
    z <- if (comparison) out$comparison else out$primary
    draw <- function() {
      limits <- range(c(z$ci95_lower, z$ci95_upper, 1))
      par(mar = c(5, 11, 3, 2), oma = c(if (comparison) 2 else 0, 0, 0, 0))
      plot(NA_real_, NA_real_, xlim = limits, ylim = c(0.5, 6.5), log = "x", yaxt = "n",
        xlab = "Odds ratio per 1 SD residual (ordinary 95% CI)", ylab = "",
        main = if (comparison) "Composition adjustment comparison" else "Composition-adjusted residual associations")
      axis(2, at = 6:1, labels = gsub("_", " ", s$axes), las = 1)
      abline(v = 1, lty = 2, col = "grey60")
      y <- 7 - match(z$axis, s$axes)
      if (comparison) y <- y + ifelse(z$model == "base_model", 0.12, -0.12)
      colour <- if (comparison) ifelse(z$model == "base_model", "#777777", "#21618C") else "#21618C"
      segments(z$ci95_lower, y, z$ci95_upper, y, col = colour, lwd = 1.5)
      points(z$odds_ratio, y, pch = if (comparison) ifelse(z$model == "base_model", 1, 16) else 16, col = colour)
      if (comparison) mtext("Open grey: base model     Filled blue: composition-adjusted v2",
        side = 1, outer = TRUE, line = 0.6, cex = 0.85)
    }
    pdf(path("figures/stage05_primary", paste0(name, ".pdf")), width = 10, height = 6)
    tryCatch(draw(), finally = dev.off())
    png(path("figures/stage05_primary", paste0(name, ".png")), width = 2000, height = 1200, res = 200)
    tryCatch(draw(), finally = dev.off())
  }
}
