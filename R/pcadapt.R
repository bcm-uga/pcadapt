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
#' @param LD.clumping a logical value indicating whether LD clumping should be
#' performed.
#' @param pca.only a logical value indicating whether multiple linear regression
#' should be performed.
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
                   min.maf = 0.05, 
                   LD.clumping = FALSE,
                   pca.only = FALSE) {
  
  #############################################
  ########## test arguments and init ##########
  #############################################
  
  if (!missing(input)) {
    
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
    
    obj.pca <- iram(input, K = K, min.maf = min.maf, LD.clumping = LD.clumping)
    res <- get_statistics(as.matrix(obj.pca$zscores), 
                          method = method, 
                          values = obj.pca$d,
                          pass = obj.pca$pass)
    output <- list(scores = obj.pca$u,
                   singular.values = sqrt(obj.pca$d * nrow(obj.pca$v) / (nrow(obj.pca$u) - 1)),
                   loadings = obj.pca$v,
                   zscores = obj.pca$zscores,
                   maf = obj.pca$maf,
                   chi2.stat = res$chi2.stat,
                   gif = res$gif,
                   pvalues = res$pvalues,
                   pass = obj.pca$pass)
    class(output) <- "pcadapt"
    attr(output, "K") <- K
    attr(output, "method") <- method
    attr(output, "min.maf") <- min.maf
    return(output)
  } else {
    appDir = system.file("shiny-examples/app-pcadapt", package = "pcadapt")
    if (appDir == "") {
      stop("Could not find Shiny app in pcadapt.", call. = FALSE)
    }
    shiny::runApp(appDir, display.mode = "normal")  
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

#' pcadapt statistics
#'
#' \code{get_statistics} returns chi-squared distributed statistics. 
#'
#' @param zscores a numeric matrix containing the z-scores.
#' @param method a character string specifying the method to be used to compute
#' the p-values. Three statistics are currently available, \code{"mahalanobis"},
#' \code{"communality"}, and \code{"componentwise"}.
#' @param values a numeric vector containing the singular values.
#' 
#' @return The returned value is a list containing the test statistics and the associated p-values.
#' 
#' @importFrom stats median na.omit pchisq qchisq
#' @importFrom MASS cov.rob
#' @importFrom utils head
#' 
#' @export
#'
get_statistics = function(zscores, 
                          method = c("mahalanobis", 
                                     "communality", 
                                     "componentwise"),
                          values, 
                          pass = rep(TRUE, nrow(zscores))) {
  nSNP <- nrow(zscores)
  K <- ncol(zscores)
  if (method == "mahalanobis") {
    res <- rep(NA, nSNP)
    if (K == 1) {
      one.d.cov <- as.vector(MASS::cov.rob(zscores[pass])) 
      res <- (zscores - one.d.cov$center)^2 / one.d.cov$cov[1]
    } else if (K > 1) {
      ogk <- covRob_cpp(zscores[pass, ])
      res[pass] <- ogk$dist
    }
    gif <- median(res, na.rm = TRUE) / qchisq(0.5, df = K)
    res.gif <- res / gif
    pval <- as.numeric(pchisq(res.gif, df = K, lower.tail = FALSE))
  } else if (method == "communality") {
    res <- sapply(1:nSNP, FUN = function(h) {sum(zscores[h, ]^2 * values^2 / nSNP)})
    c <- sum(values^2) / K
    gif <- median(res * nSNP / c, na.rm = TRUE) / qchisq(0.5, df = K)
    res.gif <- res * nSNP / (c * gif)
    pval <- pchisq(res.gif, df = K, lower.tail = FALSE)
  } else if (method == "componentwise"){
    res <- apply(zscores, MARGIN = 2, FUN = function(h) {h^2})
    gif <- sapply(1:K, FUN = function(h) {median(zscores[, h]^2, na.rm = TRUE) / qchisq(0.5, df = 1)})
    res.gif = res / gif
    pval <- NULL
    for (k in 1:K){
      pval <- cbind(pval, pchisq(res.gif[, k], df = 1, lower.tail = FALSE))
    }
  }
  return(list(stat = res, gif = gif, chi2.stat = res.gif, pvalues = pval))
}

