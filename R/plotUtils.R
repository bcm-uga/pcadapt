#' pcadapt visualization tool
#'
#' \code{plot.pcadapt} is a method designed for objects of class \code{pcadapt}.
#' It provides a plotting utile for quick visualization of \code{pcadapt} 
#' objects. Different options are currently available : \code{"screeplot"}, 
#' \code{"scores"}, \code{"stat.distribution"}, \code{"manhattan"} and 
#' \code{"qqplot"}. \code{"screeplot"} shows the decay of the genotype matrix 
#' singular values and provides a figure to help with the choice of \code{K}.
#' \code{"scores"} plots the projection of the individuals onto the first two 
#' principal components. \code{"stat.distribution"} displays the histogram of 
#' the selected test statistics, as well as the estimated distribution for the 
#' neutral SNPs. \code{"manhattan"} draws the Manhattan plot of the p-values 
#' associated with the statistic of interest. \code{"qqplot"} draws a Q-Q plot 
#' of the p-values associated with the statistic of interest.
#'
#' @param x an object of class "pcadapt" generated with \code{pcadapt}.
#' @param ... \dots
#' @param option a character string specifying the figures to be displayed. If 
#' \code{NULL} (the default), all three plots are printed.
#' @param i an integer indicating onto which principal component the individuals
#' are projected when the "scores" option is chosen.
#' Default value is set to \code{1}.
#' @param j an integer indicating onto which principal component the individuals 
#' are projected when the "scores" option is chosen.
#' Default value is set to \code{2}.
#' @param pop a list of integers or strings specifying which subpopulation the 
#' individuals belong to.
#' @param col a list of colors to be used in the score plot.
#' @param chr.info a list containing the chromosome information for each marker.
#' @param snp.info a list containing the names of all genetic markers present in
#' the input.
#' @param plt.pkg a character string specifying the package to be used to 
#' display the graphical outputs. Use \code{"plotly"} for interactive plots, or 
#' \code{"ggplot"} for static plots.
#' @param K an integer specifying the principal component of interest. \code{K} 
#' has to be specified only when using the \code{"componentwise"} method.
#' 
#' @examples
#' ## see ?pcadapt for examples
#'
#' @method plot pcadapt
#'
#' @export
plot.pcadapt = function(x, 
                        ..., 
                        option = "manhattan", 
                        i = 1, 
                        j = 2, 
                        pop, 
                        col,
                        chr.info = NULL,
                        snp.info = NULL,
                        plt.pkg = "ggplot",
                        K = NULL) {
  
  if (!(option %in% c("screeplot", 
                      "scores", 
                      "manhattan", 
                      "qqplot", 
                      "stat.distribution"))){
    warning(paste("Plotting option", 
                  option, 
                  "not valid, options currently available are: 'screeplot', 'scores', 'manhattan', 'qqplot', 'stat.distribution'."))
  } else {
    if (option == "screeplot"){
      scree_plot(x, K)
    } else if (option == "scores"){
      
      if (missing(pop)){
        score_plot(x, i, j, plt.pkg = plt.pkg)
      } else {
        if (missing(col)){
          score_plot(x, i, j, pop, plt.pkg = plt.pkg)
        } else {
          score_plot(x, i, j, pop, col, plt.pkg = plt.pkg)
        }
      }
      
    } else if (option == "stat.distribution"){
      if ((attr(x, "method") %in% c("mahalanobis", "communality")) == FALSE){
        if (is.null(K)){
          warning("K has to be specified.")
        } else {
          hist_plot(x, K)
        }
      } else {
        hist_plot(x, 1)
      }
    } else if (option == "manhattan"){
      if ((attr(x,"method") %in% c("mahalanobis", "communality")) == FALSE){
        if (is.null(K)){
          warning("K has to be specified.")
        } else {
          manhattan_plot(x,
                         chr.info = chr.info, 
                         snp.info = snp.info,
                         plt.pkg = plt.pkg, 
                         K = K)
        }
      } else {
        manhattan_plot(x, 
                       chr.info = chr.info, 
                       snp.info = snp.info,
                       plt.pkg = plt.pkg)
      }
    } else if (option == "qqplot"){
      if ((attr(x, "method") %in% c("mahalanobis", "communality")) == FALSE){
        if (is.null(K)){
          warning("K has to be specified")
        } else{
          qq_plot(x, K)
        }
      } else {
        qq_plot(x, K = 1)
      }
    }
  }
}