#' Principal Components Analysis Scree Plot
#'
#' \code{scree_plot} plots the scee plot associated with the principal components analysis performed on the dataset.
#' NB : \code{pcadapt} has to be run on the dataset in order to get an output readable by \code{plot.screePlot}
#'
#' @param x an output from \code{pcadapt} containing the singular values.
#' @param K an integer specifying the number of components to take into account in the scree plot.
#'
#' @examples
#' ## see ?fastpcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom ggplot2 qplot geom_line guides ggtitle
#' @importFrom plotly plot_ly layout
#' 
#' @export
#'
scree_plot = function(x, K) {
  
  if (is.null(K)) {
    m <- attr(x, "K")
  } else {
    m <- K
  }
  
  if (m < 2) {
    warning("the scree plot is not available.")
  } else {
    p0 <- ggplot2::qplot(x = 1:m, 
                         y = (x$singular.values[1:m]) ^ 2, 
                         col = "red", 
                         xlab = "PC", 
                         ylab = "Proportion of explained variance") + 
      ggplot2::geom_line() + ggplot2::guides(colour = FALSE) +
      ggplot2::ggtitle(paste("Scree Plot - K =", m))
    print(p0)
  }
}