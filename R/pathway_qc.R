score_qc <- function(scores) {
  do.call(rbind, lapply(names(axes), function(a) {
    x <- as.numeric(scores[[a]])
    q <- data.frame(axis = a, n = length(x), missing = sum(is.na(x)),
      nonfinite = sum(!is.finite(x)), minimum = min(x), maximum = max(x),
      mean = mean(x), sd = sd(x), median = median(x), iqr = IQR(x),
      unique_values = length(unique(x)))
    q$passed <- q$n == 401L & q$missing == 0L & q$nonfinite == 0L &
      q$unique_values > 1L & is.finite(q$sd) & q$sd > 0
    q
  }))
}
primary_files <- c("metadata/neotrip_analysis_manifest_blinded_v1.tsv",
  "data/processed/pathways/neotrip_confirmatory_scores_v1.tsv.gz",
  "data/processed/pathways/neotrip_confirmatory_pair_scores_v1.tsv",
  "data/processed/pathways/neotrip_model_inputs_blinded_v1.tsv")
leakage_qc <- function() {
  result <- do.call(rbind, lapply(primary_files, function(f) {
    present <- file.exists(path(f))
    fields <- if (present) forbidden_fields(read_tsv_header(path(f))) else character()
    data.frame(file = f, exists = present, forbidden_fields = length(fields),
      field_names = paste(fields, collapse = ";"), passed = present && !length(fields))
  }))
  write_tsv(result, path("results/stage03_pathways/leakage_check.tsv"))
  if (any(result$forbidden_fields > 0L)) stage03_stop("Outcome fields detected", "STOP_REVIEW")
  result
}
pathway_figures <- function(scores) {
  draw <- function() {
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    for (a in names(axes)) {
      groups <- factor(scores$timepoint, levels = c("Baseline", "D1C2"))
      boxplot(scores[[a]] ~ groups, main = gsub("_", " ", a), xlab = "",
        ylab = "Centered singscore", col = c("#AACDE5", "#EFCB95"), outline = TRUE)
      abline(h = 0, col = "grey60", lty = 3)
      mtext(paste("n =", sum(groups == "Baseline"), "/", sum(groups == "D1C2")), side = 1, line = 2.5, cex = 0.75)
    }
    mtext("Confirmatory pathway scores: technical distributions", outer = TRUE, cex = 1.1)
  }
  pdf(path("figures/stage03_pathways/score_distributions_v1.pdf"), width = 12, height = 7)
  tryCatch(draw(), finally = dev.off())
  png(path("figures/stage03_pathways/score_distributions_v1.png"), width = 2400, height = 1400, res = 200)
  tryCatch(draw(), finally = dev.off())
}
