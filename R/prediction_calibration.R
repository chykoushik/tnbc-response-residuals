calibration07 <- function(y, p) {
  eta <- qlogis(pmin(1 - 1e-6, pmax(1e-6, p)))
  warnings <- character()
  fits <- withCallingHandlers(list(slope = glm(y ~ eta, family = binomial()),
    intercept = glm(y ~ 1 + offset(eta), family = binomial())),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  result <- c(calibration_intercept = unname(coef(fits$intercept)[1]),
    calibration_slope = unname(coef(fits$slope)[2]))
  assert07(!length(warnings) && all(vapply(fits, function(f) f$converged, logical(1))) &&
    all(is.finite(result)), "Calibration estimation failed", "MODIFY_REVIEW")
  result
}
figures07 <- function(scores, metrics, differences) {
  colors <- c("#444444", "#0072B2", "#009E73", "#D55E00", "#CC79A7")
  render <- function(name, draw) {
    for (ext in c("pdf", "png")) {
      file <- path("figures/stage07_prediction", paste0(name, ".", ext))
      if (ext == "pdf") pdf(file, width = 8, height = 6) else png(file, width = 1600, height = 1200, res = 200)
      tryCatch(draw(), finally = dev.off())
    }
  }
  means <- aggregate(predicted_probability ~ patient_id + pcr + model, scores, mean)
  bins <- do.call(rbind, lapply(names(features07), function(model) {
    z <- means[means$model == model, ]; z <- z[order(z$predicted_probability, z$patient_id), ]
    z$bin <- rep(1:5, each = 30)
    b <- aggregate(cbind(predicted_probability, pcr) ~ bin, z, mean)
    b$model <- model; b$n <- 30L; b
  }))
  write07(means, "tables/stage07_prediction/calibration_patient_means.tsv")
  write07(bins, "tables/stage07_prediction/calibration_display_bins.tsv")
  render("calibration_plot", function() {
    par(mar = c(5, 5, 4, 2))
    plot(0:1, 0:1, type = "n", xlab = "Mean OOF predicted probability", ylab = "Observed pCR proportion",
      main = "Calibration: patient-level mean predictions")
    abline(0, 1, lty = 2, col = "grey60")
    for (i in 1:5) { b <- bins[bins$model == names(features07)[i], ]; lines(b$predicted_probability, b$pcr, type = "b", col = colors[i], pch = 15 + i) }
    legend("topleft", names(features07), col = colors, lty = 1, pch = 16:20, bty = "n")
    mtext("Display only: 10-repeat mean per patient; five equal-count bins per model", side = 3, cex = .75)
  })
  render("auc_by_repeat", function() {
    boxplot(roc_auc ~ model, metrics, col = adjustcolor(colors, .3), ylim = range(c(.5, metrics$roc_auc)),
      xlab = "Prediction model", ylab = "Repeat-level ROC-AUC", main = "OOF ROC-AUC across 10 repeats", outline = FALSE)
    for (i in 1:5) points(i + seq(-.12, .12, length.out = 10), metrics$roc_auc[metrics$model == names(features07)[i]], pch = 16, col = colors[i])
    abline(h = .5, lty = 2, col = "grey60")
  })
  render("residual_increment", function() {
    par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
    for (metric in c("auc_difference", "brier_difference")) {
      d <- differences[differences$comparison %in% c("P3 - P1", "P4 - P2"), ]
      limits <- range(c(0, d[[metric]])); if (diff(limits) == 0) limits <- c(-.001, .001)
      plot(c(1, 2), limits, type = "n", xaxt = "n", xlim = c(.7, 2.3), xlab = "Residual addition",
        ylab = if (metric == "auc_difference") "ROC-AUC difference" else "Brier score difference",
        main = if (metric == "auc_difference") "Positive favors residuals" else "Negative favors residuals")
      axis(1, 1:2, c("P3 - P1", "P4 - P2")); abline(h = 0, lty = 2, col = "grey50")
      for (r in 1:10) {
        z <- d[d[["repeat"]] == r, ]; z <- z[match(c("P3 - P1", "P4 - P2"), z$comparison), ]
        lines(1:2, z[[metric]], col = adjustcolor("#0072B2", .45), type = "b", pch = 16)
      }
    }
  })
}
