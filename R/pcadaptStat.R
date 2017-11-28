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
                          values) {
  nSNP <- nrow(zscores)
  K <- ncol(zscores)
  if (method == "mahalanobis") {
    res <- rep(NA, nSNP)
    not.NA <- which(!is.na(zscores[, 1]))
    if (K == 1) {
      one.d.cov <- as.vector(MASS::cov.rob(zsc[not.NA])) 
      res <- (zscores - one.d.cov$center)^2 / one.d.cov$cov[1]
    } else if (K > 1) {
      ogk <- covRob_cpp(zscores[not.NA, ])
      res[not.NA] <- ogk$dist
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
