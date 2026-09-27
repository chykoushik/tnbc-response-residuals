eligibility06 <- function() {
  context <- analysis_context06()
  p <- context$primary
  stable <- read06("results/stage06_robustness/five_sensitivity_stability.tsv")
  permutation <- read06("results/stage06_robustness/permutation_summary.tsv")
  influence <- read06("results/stage06_robustness/influence_summary.tsv")
  output <- do.call(rbind, lapply(axes06, function(a) {
    primary <- p[p$axis == a, ]
    s <- stable[stable$axis == a, ]
    perm <- permutation[permutation$axis == a, ]
    inf <- influence[influence$axis == a, ]
    assert06(nrow(s) == 5 && setequal(s$sensitivity_id, variants06) && nrow(perm) == 1 && nrow(inf) == 1,
      "Eligibility inputs incomplete")
    row <- data.frame(axis = a, prespecified_direction = primary$prespecified_direction,
      primary_OR = primary$odds_ratio, primary_95CI = sprintf("%.6f to %.6f", primary$ci95_lower, primary$ci95_upper),
      primary_bonferroni_p = primary$p_bonferroni,
      primary_direction_gate = primary$direction_matches,
      primary_effect_size_gate = primary$effect_size_gate,
      primary_multiplicity_gate = primary$multiplicity_gate,
      sensitivity_stable_count = sum(s$stable), sensitivity_4_of_5_gate = sum(s$stable) >= 4,
      permutation_p = perm$permutation_p, permutation_95th_percentile = perm$permutation_95th_percentile,
      permutation_gate = perm$permutation_gate, max_relative_influence = inf$max_relative_influence,
      influence_gate = inf$influence_gate)
    row$full_transportability_eligible <- all(unlist(row[c("primary_direction_gate", "primary_effect_size_gate",
      "primary_multiplicity_gate", "sensitivity_4_of_5_gate", "permutation_gate", "influence_gate")]))
    row$classification <- if (row$full_transportability_eligible) "FULL_ELIGIBLE" else
      if (row$primary_direction_gate && row$sensitivity_4_of_5_gate) "HYPOTHESIS_GENERATING" else "NOT_ROBUST"
    row
  }))
  assert06(!any(output$full_transportability_eligible & !output$primary_multiplicity_gate),
    "Sensitivity cannot replace primary multiplicity")
  write06(output, "results/stage06_robustness/stage06_eligibility.tsv")
  write06(output, "tables/stage06_robustness/stage06_eligibility.tsv")
  eligible <- output$axis[output$full_transportability_eligible]
  jsonlite::write_json(list(eligible_pathways = eligible, eligible_count = length(eligible),
    decision = if (length(eligible)) "ELIGIBLE_PENDING_SEPARATE_STAGE" else "SKIP",
    reason = if (length(eligible)) "Frozen full eligibility criteria met; no external data opened in Stage 6." else
      "No confirmatory residual pathway satisfies the frozen full eligibility rule."),
    path("results/stage06_robustness/external_transportability_gate.json"), auto_unbox = TRUE, pretty = TRUE)
  output
}
figures06 <- function() {
  five <- read06("results/stage06_robustness/five_primary_sensitivities.tsv")
  primary <- read06("results/stage05_primary/primary_residual_associations_v2.tsv")
  permutation <- read06("results/stage06_robustness/permutation_nulls.tsv.gz")
  perm_summary <- read06("results/stage06_robustness/permutation_summary.tsv")
  influence <- read06("results/stage06_robustness/influence_summary.tsv")
  random <- read06("results/stage06_robustness/random_gene_null_results.tsv.gz")
  output <- function(name, draw) {
    pdf(path("figures/stage06_robustness", paste0(name, ".pdf")), width = 13, height = 8)
    tryCatch(draw(), finally = dev.off())
    png(path("figures/stage06_robustness", paste0(name, ".png")), width = 2600, height = 1600, res = 200)
    tryCatch(draw(), finally = dev.off())
  }
  output("five_sensitivity_comparison", function() {
    par(mfrow = c(2, 3), mar = c(4, 5, 3, 1), oma = c(0, 0, 2, 0))
    for (a in axes06) {
      z <- five[five$axis == a, ]; z <- z[match(variants06, z$sensitivity_id), ]
      p <- primary[primary$axis == a, ]
      plot(NA_real_, NA_real_, xlim = range(c(z$ci95_lower, z$ci95_upper, 1, p$odds_ratio)),
        ylim = c(.5, 5.5), log = "x", yaxt = "n", main = gsub("_", " ", a),
        xlab = "OR per 1 SD residual (95% CI)", ylab = "")
      axis(2, at = 5:1, labels = c("S1", "S2", "S3", "S4", "S5"), las = 1)
      segments(z$ci95_lower, 5:1, z$ci95_upper, 5:1, col = "#276A98")
      points(z$odds_ratio, 5:1, pch = 16, col = "#276A98")
      abline(v = 1, lty = 3, col = "grey50"); abline(v = p$odds_ratio, lty = 2, col = "#A34A28")
    }
    mtext("Five frozen sensitivities; dashed orange line: primary OR", outer = TRUE, cex = 1.1)
  })
  output("permutation_distributions", function() {
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    for (a in axes06) {
      z <- permutation$directional_statistic[permutation$axis == a]
      p <- perm_summary[perm_summary$axis == a, ]
      hist(z, breaks = 35, xlim = range(c(z, p$observed_T)), main = gsub("_", " ", a),
        xlab = "Prespecified directional beta", col = "#AACDE5", border = "white")
      abline(v = p$observed_T, col = "#A34A28", lwd = 2)
      abline(v = p$permutation_95th_percentile, col = "grey35", lty = 2)
    }
    mtext("Permutation nulls: orange observed statistic; dashed 95th percentile", outer = TRUE, cex = 1.1)
  })
  output("maximum_patient_influence", function() {
    par(mar = c(5, 11, 3, 1))
    z <- influence[match(axes06, influence$axis), ]
    plot(z$max_relative_influence, 6:1, yaxt = "n", pch = 16, col = "#276A98",
      xlim = range(c(0, .25, z$max_relative_influence), na.rm = TRUE), ylim = c(.5, 6.5),
      xlab = "Maximum relative beta change", ylab = "", main = "Leave-one-out influence")
    axis(2, at = 6:1, labels = gsub("_", " ", axes06), las = 1)
    abline(v = .25, lty = 2, col = "#A34A28")
  })
  output("matched_random_gene_distributions", function() {
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    for (a in axes06) {
      z <- abs(random$beta[random$axis == a]); b <- abs(primary$beta[primary$axis == a])
      hist(z, breaks = "FD", xlim = range(c(z, b)), main = gsub("_", " ", a),
        xlab = "Absolute residual beta", col = "#AACDE5", border = "white")
      abline(v = b, col = "#A34A28", lwd = 2)
    }
    mtext("Expression-matched random-gene controls; orange: primary absolute beta", outer = TRUE, cex = 1.1)
  })
}
