random_sets06 <- function(m) {
  file <- "metadata/sensitivity/random_matched_gene_sets_v1.tsv.gz"
  means <- rowMeans(m)
  sets <- read06("metadata/pathways/confirmatory_gene_sets_v1.tsv", blinded = TRUE)
  eligible <- !rownames(m) %in% sets$gene_symbol
  breaks <- c(-Inf, quantile(means[eligible], seq(.05, .95, .05), names = FALSE, type = 7), Inf)
  assert06(!anyDuplicated(breaks), "Mean-expression quantile bins ambiguous", "MODIFY_REVIEW")
  bins <- cut(means, breaks = breaks, labels = FALSE, include.lowest = TRUE, right = TRUE)
  names(bins) <- rownames(m)
  write06(data.frame(gene_symbol = rownames(m), mean_expression = means, expression_bin = bins,
    eligible = eligible), "metadata/sensitivity/random_gene_expression_bins_v1.tsv.gz")
  write06(data.frame(boundary = seq_along(breaks), value = breaks),
    "metadata/sensitivity/random_gene_bin_boundaries_v1.tsv")
  if (file.exists(path(file))) {
    assert06(identical(readLines(path(paste0(file, ".sha256"))), sha256(path(file))), "Null memberships changed")
    return(read06(file, blinded = TRUE))
  }
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(319643)
  output <- list()
  for (a in axes06) {
    target <- intersect(sets$gene_symbol[sets$axis == a], rownames(m))
    counts <- table(factor(bins[target], levels = 1:20))
    for (j in 1:100) {
      chosen <- unlist(lapply(1:20, function(bin) {
        pool <- rownames(m)[eligible & bins == bin]
        assert06(length(pool) >= counts[bin], "Insufficient matched-null pool", "MODIFY_REVIEW")
        if (!counts[bin]) character() else sample(pool, counts[bin], replace = FALSE)
      }), use.names = FALSE)
      assert06(length(chosen) == length(target) && !anyDuplicated(chosen), "Null gene count differs")
      output[[paste(a, j)]] <- data.frame(axis = a, null_set = j, gene_symbol = chosen, expression_bin = bins[chosen])
    }
  }
  out <- do.call(rbind, output)
  write06(out, file)
  writeLines(sha256(path(file)), path(paste0(file, ".sha256")))
  out
}
random_residuals06 <- function(m, context) {
  sets <- random_sets06(m)
  ranks <- singscore::rankGenes(m, tiesMethod = "min")
  b <- match(context$pairs$baseline_sample_id, colnames(m))
  d <- match(context$pairs$d1c2_sample_id, colnames(m))
  dir.create(path("data/processed/sensitivity/null_checkpoints"), recursive = TRUE, showWarnings = FALSE)
  all <- audits <- list()
  for (a in axes06) for (j in 1:100) {
    key <- sprintf("%s_%03d", a, j)
    file <- paste0("data/processed/sensitivity/null_checkpoints/", key, ".rds")
    if (file.exists(path(file))) {
      z <- readRDS(path(file))
      assert06(identical(z$membership_sha256, sha256(path("metadata/sensitivity/random_matched_gene_sets_v1.tsv.gz"))),
        "Null checkpoint provenance differs")
    } else {
      genes <- sets$gene_symbol[sets$axis == a & sets$null_set == j]
      score <- singscore::simpleScore(ranks, upSet = genes, centerScore = TRUE, knownDirection = TRUE)$TotalScore
      x <- context$x
      x[[paste0("baseline_", a)]] <- score[b]
      x[[paste0("delta_", a)]] <- score[d] - score[b]
      z <- ridge_target06(x, a, context, paste0("random_", key), standardize_residual = FALSE)
      z$residuals$null_set <- j
      z$membership_sha256 <- sha256(path("metadata/sensitivity/random_matched_gene_sets_v1.tsv.gz"))
      saveRDS(z, path(file))
    }
    all[[key]] <- z$residuals
    audits[[key]] <- z$audit
    if (j %% 10 == 0L) {
      cat(a, j, "\n")
      flush.console()
    }
  }
  raw <- do.call(rbind, all)
  write06(raw, "data/processed/sensitivity/random_null_residuals_raw_v1.tsv.gz")
  freeze_yaml06(list(timestamp = utc_now(), clinical_labels_accessed = FALSE,
    files = hashes06(c("metadata/sensitivity/random_matched_gene_sets_v1.tsv.gz",
      "data/processed/sensitivity/random_null_residuals_raw_v1.tsv.gz"))),
    "metadata/sensitivity/random_null_raw_freeze_v1.yml")
  groups <- interaction(raw$axis, raw$null_set, drop = TRUE)
  raw$residual_z <- ave(raw$residual, groups, FUN = function(z) as.numeric(scale(z)))
  raw$baseline_z <- ave(raw$baseline_score, groups, FUN = function(z) as.numeric(scale(z)))
  assert06(all(is.finite(raw$residual_z)) && nrow(raw) == 90000, "Null residual coverage failed")
  write06(raw, "data/processed/sensitivity/random_null_residuals_v1.tsv.gz")
  write06(do.call(rbind, audits), "results/stage06_robustness/random_null_crossfit_audit.tsv.gz")
  invisible(NULL)
}
random_outcomes06 <- function(context) {
  null <- read06("data/processed/sensitivity/random_null_residuals_v1.tsv.gz", blinded = TRUE)
  rows <- list()
  for (a in axes06) for (j in 1:100) {
    z <- null[null$axis == a & null$null_set == j, ]
    z <- z[match(context$x$patient_id, z$patient_id), ]
    assert06(identical(z$patient_id, context$x$patient_id), "Random null patient mismatch")
    frame <- outcome_frame06(context$x, a)
    frame$residual <- z$residual_z; frame$baseline <- z$baseline_z
    effect <- glm_effect06(fit_glm06(frame))
    rows[[paste(a, j)]] <- data.frame(axis = a, null_set = j, effect)
  }
  out <- do.call(rbind, rows)
  write06(out, "results/stage06_robustness/random_gene_null_results.tsv.gz")
  summary <- do.call(rbind, lapply(axes06, function(a) {
    values <- abs(out$beta[out$axis == a]); observed <- abs(context$primary$beta[context$primary$axis == a])
    data.frame(axis = a, null_sets = length(values), primary_absolute_beta = observed,
      null_absolute_beta_95th = unname(quantile(values, .95, type = 7)),
      primary_percentile = 100 * mean(values <= observed),
      null_absolute_beta_at_least_primary = sum(values >= observed))
  }))
  write06(summary, "results/stage06_robustness/random_gene_null_summary.tsv")
}
