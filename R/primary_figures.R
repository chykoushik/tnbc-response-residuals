figures05 <- function(inference, spec) {
  device <- function(name, draw) {
    pdf(path("figures/stage05_primary", paste0(name, ".pdf")), width = 9, height = 5.5)
    tryCatch(draw(), finally = dev.off())
    png(path("figures/stage05_primary", paste0(name, ".png")), width = 1800, height = 1100, res = 200)
    tryCatch(draw(), finally = dev.off())
  }
  draw <- function(x, comparison = FALSE) {
    limits <- range(c(x$ci95_lower, x$ci95_upper, 1))
    if (any(!is.finite(limits)) || any(limits <= 0)) stop05("Forest interval not finite", "MODIFY_REVIEW")
    par(mar = c(5, 11, 3, 2))
    plot(NA, xlim = limits, ylim = c(0.5, 6.5), log = "x", yaxt = "n",
      xlab = "Odds ratio per 1 SD residual (ordinary 95% CI)", ylab = "",
      main = if (comparison) "Composition adjustment comparison" else "Composition-adjusted residual associations")
    axis(2, at = 6:1, labels = gsub("_", " ", spec$axes), las = 1)
    abline(v = 1, lty = 2, col = "grey60")
    y <- 7 - match(x$axis, spec$axes)
    if (comparison) y <- y + ifelse(x$model == "base_model", 0.12, -0.12)
    colour <- if (comparison) ifelse(x$model == "base_model", "#777777", "#21618C") else "#21618C"
    segments(x$ci95_lower, y, x$ci95_upper, y, col = colour)
    points(x$odds_ratio, y, pch = if (comparison) ifelse(x$model == "base_model", 1, 16) else 16, col = colour)
    if (comparison) legend("bottomright", c("Base model", "Composition-adjusted model"),
      col = c("#777777", "#21618C"), pch = c(1, 16), bty = "n", cex = 0.8)
  }
  device("primary_residual_forest", function() draw(inference$primary))
  device("composition_adjustment_comparison", function() draw(inference$comparison, TRUE))
}
