#' p-values Q-Q Plot
#'
#' \code{pvalqq.plotting} plots a Q-Q plot of the p-values computed.
#'
#' @param x an output from \code{outlier} containing the p-values of interest.
#' @param K an integer specifying which principal component to display when \code{method="componentwise"}.
#' @param threshold a real number between \code{0} et \code{1}.
#'
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom ggplot2 qplot geom_vline guides geom_abline ggtitle aes
#'
#' @export
#'
pvalqq.plotting = function(x, K, threshold){
  if (attr(x, "method") == "componentwise"){
    sorted.pval <- sort(x$pvalues[x$maf >= attr(x, "min.maf"), K])
  } else {
    sorted.pval <- sort(x$pvalues[x$maf >= attr(x, "min.maf")])
  }
  p <- length(sorted.pval)
  expected.p <- 1:p / p
  p0 <- ggplot2::qplot(-log10(expected.p), -log10(sorted.pval), col = "red", xlab = "Expected -log10(p-values)", ylab = "Observed -log10(p-values)") + 
    ggplot2::geom_abline()
  if (!missing(threshold)){
    q <- floor(threshold * p)
    pval.thresh <- expected.p[q]
    p0 <- p0 + ggplot2::geom_vline(aes(xintercept = -log10(pval.thresh)), colour = "blue")
  }
  p0 <- p0 + ggplot2::ggtitle("Q-Q plot") + ggplot2::guides(colour = FALSE)
  print(p0)
}
