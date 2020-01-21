#' p-values Q-Q Plot
#'
#' \code{qq_plot} plots a Q-Q plot of the p-values computed.
#'
#' @param x an output from \code{outlier} containing the p-values of interest.
#' @param K an integer specifying which principal component to display when \code{method="componentwise"}.
#'
#' @keywords internal
#'
#' @import ggplot2
#'
#' @export
#'
qq_plot = function(x, K = 1) {
  
  if (attr(x, "method") == "componentwise") {
    sorted.pval <- sort(x$pvalues[x$maf >= attr(x, "min.maf"), K])
  } else {
    sorted.pval <- sort(x$pvalues[x$maf >= attr(x, "min.maf")])
  }
  expected.p <- stats::ppoints(length(sorted.pval))
  qplot(-log10(expected.p), -log10(sorted.pval), col = "red", 
        xlab = "Expected -log10(p-values)", 
        ylab = "Observed -log10(p-values)") + 
    geom_abline() +
    ggtitle("Q-Q plot") + 
    guides(colour = FALSE)
}
