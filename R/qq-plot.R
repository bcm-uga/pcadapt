#' p-values Q-Q Plot
#'
#' \code{qq_plot} plots a Q-Q plot of the p-values computed.
#'
#' @param x an output from \code{outlier} containing the p-values of interest.
#' @param K an integer specifying which principal component to display when \code{method="componentwise"}.
#'
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom ggplot2 qplot guides geom_abline ggtitle aes
#'
#' @export
#'
qq_plot = function(x, K = 1) {
  
  if (attr(x, "method") == "componentwise") {
    sorted.pval <- sort(x$pvalues[x$maf >= attr(x, "min.maf"), K])
  } else {
    sorted.pval <- sort(x$pvalues[x$maf >= attr(x, "min.maf")])
  }
  expected.p <- 1:length(sorted.pval) / length(sorted.pval)
  ggplot2::qplot(-log10(expected.p), 
                 -log10(sorted.pval), 
                 col = "red", 
                 xlab = "Expected -log10(p-values)", 
                 ylab = "Observed -log10(p-values)") + 
    ggplot2::geom_abline() +
    ggplot2::ggtitle("Q-Q plot") + 
    ggplot2::guides(colour = FALSE)
  
}
