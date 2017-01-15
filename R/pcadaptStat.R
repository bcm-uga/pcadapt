#' pcadapt statistics
#'
#' \code{cmpt.stat} returns chi-squared distributed statistics. 
#'
#' @param zscores a real-valued matrix obtained from the multiple linear regression of the genotype on the scores.
#' @param s.v a vector containing the K larger singular values.
#' @param K an integer specifying the number of principal components to retain.
#' @param method a character string specifying the method to be used to compute
#' the p-values. Three statistics are currently available, \code{"mahalanobis"},
#' \code{"communality"}, and \code{"componentwise"}.
#' @param nSNP an integer specifying the number of genetic markers present in the data.
#' 
#' @return The returned value is a list containing the test statistics and the associated p-values.
#' 
#' @export
#'
cmpt.stat = function(zscores, s.v, K, method, nSNP, maf, min.maf){
  zsc <- zscores
  zsc[maf < min.maf] <- NA
  if (method == "mahalanobis"){
    xstat <- array(NA, nSNP)
    not.nan <- which(!is.na(apply(abs(zsc), 1, sum)))
    if (K > 1){
      xstat[not.nan] <- as.vector(robust::covRob(zsc, na.action = na.omit, estim = "pairwiseGK")$dist)
    } else if (K == 1){
      one.d.cov <- as.vector(MASS::cov.rob(zsc[not.nan, 1]))
      xstat <- (zsc[, 1] - one.d.cov$center)^2 / one.d.cov$cov[1]
    }
    gif <- median(xstat, na.rm = TRUE) / qchisq(0.5, df = K)
    ystat <- xstat / gif
    pval <- as.numeric(pchisq(ystat, df = K, lower.tail = FALSE))
  } else if (method == "communality"){
    xstat <- sapply(1:nSNP, FUN = function(h){sum(zsc[h, 1:K]^2 * s.v[1:K]^2 / (nSNP))})
    c <- sum(s.v[1:K]^2) / K
    gif <- median(xstat * nSNP / c, na.rm = TRUE) / qchisq(0.5, df = K)
    ystat <- xstat * nSNP / (c * gif)
    pval <- as.numeric(pchisq(ystat, df = K, lower.tail = FALSE))
  } else if (method == "componentwise"){
    xstat <- apply(zsc, MARGIN = 2, FUN = function(h){h^2})
    gif <- sapply(1:K, FUN = function(h){median(zsc[, h]^2, na.rm = TRUE) / qchisq(0.5, df = 1)})
    ystat = xstat / gif
    pval <- NULL
    for (k in 1:K){
      pval <- cbind(pval, pchisq(ystat[, k], df = 1, lower.tail = FALSE))
    }
  }
  return(list(stat = xstat, gif = gif, chi2.stat = ystat, pvalues = pval))
}
