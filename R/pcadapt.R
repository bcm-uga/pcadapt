#' Principal Component Analysis for outlier detection
#'
#' \code{pcadapt} performs principal component analysis and computes p-values to
#' test for outliers. The test for outliers is based on the correlations between
#' genetic variation and the first \code{K} principal components. \code{pcadapt}
#' also handles Pool-seq data for which the statistical analysis is performed on
#' the genetic markers frequencies. Returns an object of class \code{pcadapt}.
#'
#' @details First, a principal component analysis is performed on the scaled and 
#' centered genotype data. To account for missing data, the correlation matrix 
#' between individuals is computed using only the markers available for each
#' pair of individuals. Depending on the specified \code{method}, different test 
#' statistics can be used.
#'
#' \code{mahalanobis} (default): the robust Mahalanobis distance is computed for 
#' each genetic marker using a robust estimate of both mean and covariance 
#' matrix between the \code{K} vectors of z-scores.
#'
#' \code{communality}: the communality statistic measures the proportion of 
#' variance explained by the first \code{K} PCs.
#'
#' \code{componentwise}: returns a matrix of z-scores.
#'
#' To compute p-values, test statistics (\code{stat}) are divided by a genomic 
#' inflation factor (\code{gif}) when \code{method="mahalanobis"}. When 
#' \code{method="communality"}, the test statistic is first multiplied by 
#' \code{K} and divided by the percentage of variance explained by the first 
#' \code{K} PCs before accounting for genomic inflation factor. When using 
#' \code{method="mahalanobis"} or \code{"communality"}, the scaled statistics 
#' (\code{chi2_stat}) should follow a chi-squared distribution with \code{K} 
#' degrees of freedom. When using \code{method="componentwise"}, the z-scores 
#' should follow a chi-squared distribution with \code{1} degree of freedom. For 
#' Pool-seq data, \code{pcadapt} provides p-values based on the Mahalanobis 
#' distance for each SNP.
#'
#' @param input a genotype matrix or a character string specifying the name of 
#' the file to be processed with \code{pcadapt}.
#' @param K an integer specifying the number of principal components to retain.
#' @param method a character string specifying the method to be used to compute
#' the p-values. Three statistics are currently available, \code{"mahalanobis"},
#' \code{"communality"} and \code{"componentwise"}.
#' @param min.maf a value between \code{0} and \code{0.45} specifying the 
#' threshold of minor allele frequencies above which p-values are computed.
#' 
#' @return The returned value \code{x} is an object of class \code{pcadapt}.
#' 
#' @useDynLib pcadapt
#'
#' @export
#'
pcadapt = function(input, 
                   K = 2, 
                   method = "mahalanobis", 
                   data.type = "genotype",
                   min.maf = 0.05) {
  
  #############################################
  ########## test arguments and init ##########
  #############################################
  
  if (missing(input)) {
    appDir = system.file("shiny-examples/app-pcadapt", package = "pcadapt")
    if (appDir == "") {
      stop("Could not find Shiny app in pcadapt.", call. = FALSE)
    }
    shiny::runApp(appDir, display.mode = "normal")  
  } else {
    
    if (!(class(K) %in% c("numeric", "integer")) || K <= 0){
      stop("K has to be a positive integer.")
    }
    
    if (!(method %in% c("mahalanobis", "communality", "componentwise"))) {
      warning("Unknown method. Default method will be used.")
      method <- "mahalanobis"
    }
    
    if (class(min.maf) != "numeric" || min.maf < 0 || min.maf > 0.45) {
      warning("min.maf has to be a real number between 0 and 0.45. Default 
                value will be used.")
      min.maf <- 0.05
    }
    
    if (is.character(input) && !file.exists(input)) {
      stop(paste("File", input, "does not exist."))
    } 
    
    obj.pca <- iram(input, K = K, min.maf = min.maf)
    res <- get_statistics(as.matrix(obj.pca$zscores), 
                          method = method, 
                          values = obj.pca$d)
    output <- list(scores = obj.pca$u,
                   singular.values = sqrt(obj.pca$d * nrow(obj.pca$v) / (nrow(obj.pca$u) - 1)),
                   loadings = obj.pca$v,
                   zscores = obj.pca$zscores,
                   chi2.stat = res$chi2.stat,
                   gif = res$gif,
                   pvalues = res$pvalues)
    class(output) <- "pcadapt"
    attr(output, "K") <- K
    attr(output, "method") <- method
    attr(output, "min.maf") <- min.maf
    return(output)
  } 
}

#' Shiny app
#'
#' \code{pcadapt} comes with a Shiny interface.
#'
#' @export
run.pcadapt <- function() {
  appDir <- system.file("shiny-examples", "app-pcadapt", package = "pcadapt")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `pcadapt`.", 
         call. = FALSE)
  }
  shiny::runApp(appDir, launch.browser = TRUE)
}
