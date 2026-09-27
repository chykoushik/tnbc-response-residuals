permutation06 <- function(context) {
  s <- read_spec06()
  x <- context$x; nperm <- s$permutations
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(s$permutation_seed)
  indices <- replicate(nperm, sample.int(nrow(x)))
  saveRDS(indices, path("metadata/sensitivity/permutation_indices_v1.rds"))
  dir.create(path("results/stage06_robustness/permutation_checkpoints"), recursive = TRUE, showWarnings = FALSE)
  rows <- summaries <- list()
  for (a in axes06) {
    frame <- outcome_frame06(x, a)
    m <- model.matrix(pcr ~ residual + arm + baseline + cancer + immune + stromal, frame)
    source_hash <- digest::digest(list(m, x$pcr, indices), algo = "sha256")
    file <- path("results/stage06_robustness/permutation_checkpoints", paste0(a, ".rds"))
    if (file.exists(file)) {
      saved <- readRDS(file)
      assert06(identical(saved$source_hash, source_hash), "Permutation checkpoint changed")
      z <- saved$results
    } else {
      beta <- rep(NA_real_, nperm); passed <- logical(nperm)
      for (i in seq_len(nperm)) {
        warnings <- character()
        fit <- withCallingHandlers(glm.fit(m, x$pcr[indices[, i]], family = binomial("logit")),
          warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
        beta[i] <- fit$coefficients[["residual"]]
        passed[i] <- isTRUE(fit$converged) && fit$rank == ncol(m) && all(is.finite(fit$coefficients)) && !length(warnings)
      }
      sign <- if (context$primary$prespecified_direction[context$primary$axis == a] == "positive") 1 else -1
      z <- data.frame(axis = a, permutation = seq_len(nperm), beta = beta,
        directional_statistic = sign * beta, passed = passed)
      saveRDS(list(source_hash = source_hash, results = z), file)
    }
    assert06(all(z$passed) && nrow(z) == nperm, "Permutation model failed", "MODIFY_REVIEW")
    direction <- if (context$primary$prespecified_direction[context$primary$axis == a] == "positive") 1 else -1
    observed <- direction * context$primary$beta[context$primary$axis == a]
    threshold <- unname(quantile(z$directional_statistic, .95, type = 7))
    summaries[[a]] <- data.frame(axis = a, permutations = nperm, observed_T = observed,
      permutation_95th_percentile = threshold,
      permutation_p = (1 + sum(z$directional_statistic >= observed)) / (nperm + 1),
      permutation_gate = observed > threshold, technical_failures = sum(!z$passed))
    rows[[a]] <- z
    cat(a, "permutations complete\n"); flush.console()
  }
  write06(do.call(rbind, rows), "results/stage06_robustness/permutation_nulls.tsv.gz")
  write06(do.call(rbind, summaries), "results/stage06_robustness/permutation_summary.tsv")
}
