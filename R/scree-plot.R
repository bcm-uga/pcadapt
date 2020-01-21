#' Principal Components Analysis Scree Plot
#'
#' \code{scree_plot} plots the scee plot associated with the principal components 
#' analysis performed on the dataset. NB : \code{pcadapt} has to be run on the
#' dataset in order to get an output readable by \code{plot.screePlot}
#'
#' @param x an output from \code{pcadapt} containing the singular values.
#' @param K an integer specifying the number of components to plot.
#'
#' @keywords internal
#' 
#' @import ggplot2
#' 
#' @export
#'
scree_plot = function(x, K = NULL) {
  
  if (is.null(K)) K <- attr(x, "K")
  
  if (K < 2) {
    warning("The scree plot is not available for K=1.")
  } else {
    p0 <- qplot(x = 1:K, y = (x$singular.values[1:K]) ^ 2, col = "red", 
                xlab = "PC", ylab = "Proportion of explained variance") + 
      geom_line() + guides(colour = FALSE) +
      ggtitle(paste("Scree Plot - K =", K))
    print(p0)
  }
}