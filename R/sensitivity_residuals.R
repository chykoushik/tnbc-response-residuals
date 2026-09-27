ridge_target06 <- function(x, axis, context, label, predictors = base_predictors06,
                           rule = "lambda.1se", interaction = FALSE, arm_specific = FALSE, standardize_residual = TRUE) {
  m <- x[, predictors, drop = FALSE]
  if ("arm" %in% names(m)) m$arm <- as.integer(m$arm == "CT/A")
  if (interaction) m$arm_baseline_interaction <- as.integer(x$arm == "CT/A") * x[[paste0("baseline_", axis)]]
  m <- as.matrix(m)
  assert06(all(is.finite(m)) && !length(forbidden06(colnames(m))), "Invalid blinded ridge design")
  predicted <- rep(NA_real_, nrow(x)); seen <- integer(nrow(x))
  audit <- models <- list()
  groups <- if (arm_specific) c("CT", "CT/A") else "all"
  for (group in groups) for (k in 1:10) {
    use <- if (group == "all") rep(TRUE, nrow(x)) else x$arm == group
    train <- which(context$outer$outer_fold != k & use)
    test <- which(context$outer$outer_fold == k & use)
    inner <- if (arm_specific) {
      table <- context$arm_inner
      table <- table[table$arm == group & table$outer_fold == k, ]
      table[match(x$patient_id[train], table$patient_id), ]
    } else {
      table <- context$inner[context$inner$outer_fold == k, ]
      table[match(x$patient_id[train], table$patient_id), ]
    }
    assert06(!anyNA(inner$patient_id) && identical(inner$patient_id, x$patient_id[train]) &&
      !length(intersect(x$patient_id[test], inner$patient_id)) && !length(intersect(train, test)),
      "Cross-fitting patient separation failed")
    warnings <- character()
    fit <- withCallingHandlers(glmnet::cv.glmnet(m[train, , drop = FALSE], x[[paste0("delta_", axis)]][train],
      family = "gaussian", alpha = 0, intercept = TRUE, standardize = TRUE, type.measure = "mse",
      nlambda = 100, lambda.min.ratio = 1e-4, nfolds = if (arm_specific) 5 else 10,
      foldid = inner$inner_fold, grouped = TRUE, alignment = "lambda", parallel = FALSE),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    assert06(!length(warnings) && fit$glmnet.fit$jerr == 0L, "Ridge fitting warning", "MODIFY_REVIEW")
    expected <- as.numeric(predict(fit, newx = m[test, , drop = FALSE], s = rule))
    selected <- fit[[rule]]
    co <- as.matrix(coef(fit, s = rule))
    assert06(max(abs(expected - cbind(1, m[test, , drop = FALSE]) %*% co)) < 1e-10, "Stored coefficient prediction differs")
    predicted[test] <- expected; seen[test] <- seen[test] + 1L
    key <- paste(group, k, sep = "_")
    models[[key]] <- list(train_patient_ids = x$patient_id[train], test_patient_ids = x$patient_id[test],
      inner_folds = inner, predictors = colnames(m), axis = axis, selected_rule = rule,
      selected_lambda = selected, fit = fit,
      training_data_sha256 = digest::digest(list(m[train, , drop = FALSE], x[[paste0("delta_", axis)]][train]), algo = "sha256"))
    audit[[key]] <- data.frame(variant = label, axis = axis, arm_group = group, outer_fold = k,
      n_train = length(train), n_test = length(test), lambda_min = fit$lambda.min, lambda_1se = fit$lambda.1se,
      selected_lambda = selected, selected_rule = rule, separation_violations = 0L, warnings = length(warnings))
  }
  residual <- x[[paste0("delta_", axis)]] - predicted
  assert06(all(seen == 1L) && all(is.finite(residual)) && sd(residual) > 0, "Incomplete OOF residuals")
  out <- data.frame(patient_id = x$patient_id, arm = x$arm, outer_fold = context$outer$outer_fold,
    axis = axis, baseline_score = x[[paste0("baseline_", axis)]],
    observed_delta = x[[paste0("delta_", axis)]], expected_delta = predicted, residual = residual)
  if (standardize_residual) out$residual_z <- as.numeric(scale(residual))
  assert06(!length(forbidden06(names(out))), "Forbidden residual fields")
  file <- path("models/stage06", label, paste0(axis, ".rds"))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(models, file)
  list(residuals = out, audit = do.call(rbind, audit),
    scaling = data.frame(variant = label, axis = axis, mean = mean(residual), sd = sd(residual), n = length(residual)))
}
construct_variant06 <- function(x, context, label, ...) {
  results <- lapply(axes06, function(a) ridge_target06(x, a, context, label, ...))
  output <- do.call(rbind, lapply(results, "[[", "residuals"))
  write06(output, paste0("data/processed/sensitivity/residuals_", tolower(label), ".tsv"))
  write06(do.call(rbind, lapply(results, "[[", "audit")),
    paste0("results/stage06_robustness/tuning_", tolower(label), ".tsv"))
  write06(do.call(rbind, lapply(results, "[[", "scaling")),
    paste0("metadata/sensitivity/scaling_", tolower(label), ".tsv"))
  output
}
expression06 <- function() {
  f <- "data/raw/geo/GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz"
  x <- read06(f, blinded = TRUE)
  m <- as.matrix(x[, -(1:3)])
  rownames(m) <- x$HGNC_approved_symbol
  assert06(nrow(m) == 23890 && ncol(m) == 401 && all(is.finite(m)) &&
    !anyDuplicated(rownames(m)) && !anyNA(rownames(m)), "Expression matrix invalid")
  m
}
pair_scores06 <- function(scores, context) {
  x <- context$x
  b <- match(context$pairs$baseline_sample_id, scores$sample_id)
  d <- match(context$pairs$d1c2_sample_id, scores$sample_id)
  assert06(!anyNA(c(b, d)) && all(scores$patient_id[b] == x$patient_id) &&
    all(scores$patient_id[d] == x$patient_id), "Sensitivity score pair mismatch")
  for (a in axes06) {
    x[[paste0("baseline_", a)]] <- scores[[a]][b]
    x[[paste0("delta_", a)]] <- scores[[a]][d] - scores[[a]][b]
  }
  x
}
