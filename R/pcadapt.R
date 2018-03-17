################################################################################

#' @useDynLib pcadapt, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

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
#' variance explained by the first \code{K} PCs. Deprecated in version 4.0.0.
#'
#' \code{componentwise}: returns a matrix of z-scores.
#'
#' To compute p-values, test statistics (\code{stat}) are divided by a genomic 
#' inflation factor (\code{gif}) when \code{method="mahalanobis"}. When using 
#' \code{method="mahalanobis"}, the scaled statistics 
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
#' the p-values. Two statistics are currently available, \code{"mahalanobis"},
#' and \code{"componentwise"}.
#' @param min.maf a value between \code{0} and \code{0.45} specifying the 
#' threshold of minor allele frequencies above which p-values are computed.
#' @param LD.clumping Default is \code{NULL} and doesn't use any SNP thinning.
#'   If you want to use SNP thinning, provide a named list with parameters 
#'   \code{size} and \code{thr} which corresponds respectively to the window 
#'   radius and the squared correlation threshold. A good default value would 
#'   be \code{list(size = 200, thr = 0.1)}.
#' @param pca.only a logical value indicating whether PCA results should be 
#'   returned (before computing any statistic).
#' @param ploidy Number of trials, parameter of the binomial distribution. 
#'   Default is 2, which corresponds to diploidy, such as for the human genome.
#' 
#' @return The returned value \code{x} is an object of class \code{pcadapt}.
#' 
#' @useDynLib pcadapt
#' 
#' @name pcadapt
#'
#' @export
#'
pcadapt <- function(input, 
                    K = 2, 
                    method = "mahalanobis", 
                    min.maf = 0.05, 
                    ploidy = 2,
                    LD.clumping = NULL,
                    pca.only = FALSE) {
  
  if (missing(input)) {
    appDir <- system.file("shiny-examples/app-pcadapt", package = "pcadapt")
    if (appDir == "") {
      stop("Could not find Shiny app in pcadapt.", call. = FALSE)
    }
    shiny::runApp(appDir, display.mode = "normal")  
  } else {
    UseMethod("pcadapt")
  }
}

#' @rdname pcadapt
#' @export
pcadapt.pcadapt_matrix <- function(input, 
                                   K = 2, 
                                   method = c("mahalanobis", "componentwise"), 
                                   min.maf = 0.05, 
                                   ploidy = 2,
                                   LD.clumping = NULL,
                                   pca.only = FALSE) {
  
  pcadapt0(input, K, match.arg(method), min.maf, ploidy, LD.clumping, pca.only)
}

#' @rdname pcadapt
#' @export
pcadapt.pcadapt_bed <- function(input, 
                                K = 2, 
                                method = c("mahalanobis", "componentwise"), 
                                min.maf = 0.05, 
                                ploidy = 2,
                                LD.clumping = NULL,
                                pca.only = FALSE) {
  
  # File mapping
  n <- attr(input, "n")
  p <- attr(input, "p")
  xptr <- structure(bedXPtr(input, n, p), n = n, p = p, class = "xptr_bed")
  
  pcadapt0(xptr, K, match.arg(method), min.maf, ploidy, LD.clumping, pca.only)
}

#' @rdname pcadapt
#' @export
pcadapt.pcadapt_pool <- function(input, 
                                 K = (nrow(input) - 1),
                                 method = "mahalanobis",
                                 min.maf = 0.05,
                                 ploidy = NULL,
                                 LD.clumping = NULL,
                                 pca.only = FALSE) {
  
  w <- matrix(NA, nrow = ncol(input), ncol = K)
  
  tmat <- scale(input, center = TRUE, scale = FALSE) 
  tmat[is.na(tmat)] <- 0 # mean imputation
  
  mean_freq <- attr(tmat, "scaled:center")
  mean_freq <- pmin(mean_freq, 1 - mean_freq)

  pass <- mean_freq > min.maf
  
  if (nrow(input) == 2) {
    obj.pca <- list()
    obj.pca$u <- matrix(0, nrow = 1, ncol = 2)
    obj.pca$v <- tmat[1, pass, drop = FALSE]
    obj.pca$d <- 1
  } else {
    obj.pca <- RSpectra::svds(tmat[, pass, drop = FALSE], k = K)
  }
  
  w[pass, ] <- obj.pca$v 
  res <- get_statistics(w, 
                        method = method, 
                        pass = pass)
  
  structure(
    list(
      scores = obj.pca$u,
      singular.values = sqrt(obj.pca$d * nrow(w) / (nrow(obj.pca$u) - 1)),
      loadings = w,
      zscores = w,
      af = attr(tmat, "scaled:center"),
      maf = mean_freq,
      chi2.stat = res$chi2.stat,
      stat = res$stat,
      gif = res$gif,
      pvalues = res$pvalues,
      pass = pass
    ),
    K = K, method = method, min.maf = min.maf, class = "pcadapt"
  )
}

################################################################################

#' pcadapt statistics
#'
#' \code{get_statistics} returns chi-squared distributed statistics. 
#'
#' @param zscores a numeric matrix containing the z-scores.
#' @param method a character string specifying the method to be used to compute
#' the p-values. Two statistics are currently available, \code{"mahalanobis"},
#' and \code{"componentwise"}.
#' @param pass a boolean vector.
#' 
#' @return The returned value is a list containing the test statistics and the 
#' associated p-values.
#' 
#' @importFrom stats median na.omit pchisq qchisq
#'
get_statistics = function(zscores, method, pass) {
  
  nSNP <- nrow(zscores)
  K <- ncol(zscores)
  if (method == "mahalanobis") {
    res <- rep(NA, nSNP)
    if (K == 1) {
      res[pass] <- (zscores[pass] - median(zscores[pass]))^2 
    } else if (K > 1) {
      # covRob_cpp(zscores[pass, ])
      res[pass] <- robust::covRob(zscores, na.action = na.omit, 
                                  estim = "pairwiseGK")$dist
    }
    gif <- median(res, na.rm = TRUE) / qchisq(0.5, df = K)
    res.gif <- res / gif
    pval <- as.numeric(pchisq(res.gif, df = K, lower.tail = FALSE))
  } else if (method == "communality") {
    # res <- sapply(1:nSNP, FUN = function(h) {sum(zscores[h, ]^2 * values^2 / nSNP)})
    # c <- sum(values^2) / K
    # gif <- median(res * nSNP / c, na.rm = TRUE) / qchisq(0.5, df = K)
    # res.gif <- res * nSNP / (c * gif)
    # pval <- pchisq(res.gif, df = K, lower.tail = FALSE)
  } else if (method == "componentwise") {
    res <- apply(zscores, MARGIN = 2, FUN = function(h) {h^2})
    gif <- sapply(1:K, FUN = function(h) {
      median(zscores[, h]^2, na.rm = TRUE) / qchisq(0.5, df = 1)
    })
    res.gif = res / gif
    pval <- NULL
    for (k in 1:K) {
      pval <- cbind(pval, pchisq(res.gif[, k], df = 1, lower.tail = FALSE))
    }
  }
  return(list(stat = res, gif = gif, chi2.stat = res.gif, pvalues = pval))
}

################################################################################

pcadapt0 <- function(input, K, method, min.maf, ploidy, LD.clumping, pca.only) {
  
  # Test arguments and init
  if (!(class(K) %in% c("numeric", "integer")) || K <= 0)
    stop("K has to be a positive integer.")
  
  if (class(min.maf) != "numeric" || min.maf < 0 || min.maf > 0.45) 
    stop("min.maf has to be a real number between 0 and 0.45.")
  
  # Compute PCs and z-scores    
  obj.pca <- iram_and_reg(input, K = K, 
                          min.maf = min.maf, 
                          ploidy = ploidy,
                          LD.clumping = LD.clumping)
  if (pca.only) return(obj.pca)
  
  res <- get_statistics(obj.pca$zscores, 
                        method = method, 
                        pass = obj.pca$pass)
  
  structure(
    list(
      scores = obj.pca$u,
      singular.values = sqrt(obj.pca$d * nrow(obj.pca$v) / (nrow(obj.pca$u) - 1)),
      loadings = obj.pca$v,
      zscores = obj.pca$zscores,
      af = obj.pca$af,
      maf = pmin(obj.pca$af, 1 - obj.pca$af),
      chi2.stat = res$chi2.stat,
      stat=res$stat,
      gif = res$gif,
      pvalues = res$pvalues,
      pass = obj.pca$pass
    ),
    K = K, method = method, min.maf = min.maf, class = "pcadapt"
  )
}

################################################################################

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

################################################################################