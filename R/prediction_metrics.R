metrics07 <- function(scores, x) {
  assert07(nrow(scores) == 7500 && !anyDuplicated(scores[, c("repeat", "patient_id", "model")]) &&
    all(scores$predicted_probability > 0 & scores$predicted_probability < 1), "OOF coverage invalid")
  rows <- list()
  for (r in 1:10) for (model in names(features07)) {
    z <- scores[scores[["repeat"]] == r & scores$model == model, ]
    assert07(nrow(z) == 150 && setequal(z$patient_id, x$patient_id) &&
      all(z$pcr == x$pcr[match(z$patient_id, x$patient_id)]), "OOF patient identity failed")
    y <- z$pcr; p <- z$predicted_probability; n1 <- sum(y); n0 <- sum(y == 0)
    auc <- (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
    row <- data.frame(repetition = r, model = model, roc_auc = auc,
      brier_score = mean((y - p)^2), log_loss = -mean(y * log(p) + (1 - y) * log1p(-p)),
      as.list(calibration07(y, p)))
    names(row)[1] <- "repeat"; rows[[paste(r, model)]] <- row
  }
  metrics <- do.call(rbind, rows); summary <- list()
  for (model in names(features07)) for (metric in names(metrics)[-(1:2)]) {
    v <- metrics[metrics$model == model, metric]
    summary[[paste(model, metric)]] <- data.frame(model = model, metric = metric, mean = mean(v),
      standard_deviation = sd(v), median = median(v), minimum = min(v), maximum = max(v),
      percentile_2_5 = unname(quantile(v, .025, type = 7)),
      percentile_97_5 = unname(quantile(v, .975, type = 7)), range_label = "empirical repeat range")
  }
  differences <- list()
  for (pair in list(c("P1", "P0"), c("P2", "P1"), c("P3", "P1"), c("P4", "P2"), c("P4", "P3"))) {
    later <- metrics[metrics$model == pair[1], ]; earlier <- metrics[metrics$model == pair[2], ]
    assert07(identical(later[["repeat"]], earlier[["repeat"]]), "Metric pairing invalid")
    row <- data.frame(comparison = paste(pair, collapse = " - "), repetition = later[["repeat"]],
      auc_difference = later$roc_auc - earlier$roc_auc, brier_difference = later$brier_score - earlier$brier_score,
      logloss_difference = later$log_loss - earlier$log_loss)
    names(row)[2] <- "repeat"; differences[[pair[1:2] |> paste(collapse = "_")]] <- row
  }
  differences <- do.call(rbind, differences); paired <- list()
  for (comparison in unique(differences$comparison)) for (metric in names(differences)[3:5]) {
    v <- differences[differences$comparison == comparison, metric]
    paired[[paste(comparison, metric)]] <- data.frame(comparison = comparison, metric = metric,
      mean_difference = mean(v), median_difference = median(v),
      repeats_favoring_later = sum(if (metric == "auc_difference") v > 0 else v < 0),
      minimum = min(v), maximum = max(v))
  }
  paired <- do.call(rbind, paired)
  descriptive <- paired[paired$comparison %in% c("P3 - P1", "P4 - P2"), ]
  descriptive$mean_improves <- ifelse(descriptive$metric == "auc_difference",
    descriptive$mean_difference > 0, descriptive$mean_difference < 0)
  descriptive$interpretation <- "Secondary descriptive result; no significance threshold or biological success criterion"
  write07(metrics, "results/stage07_prediction/repeat_metrics.tsv")
  write07(do.call(rbind, summary), "results/stage07_prediction/prediction_summary.tsv")
  write07(differences, "results/stage07_prediction/paired_metric_differences.tsv")
  write07(paired, "results/stage07_prediction/paired_metric_summary.tsv")
  write07(descriptive, "tables/stage07_prediction/residual_predictive_description.tsv")
  figures07(scores, metrics, differences)
  invisible(metrics)
}
