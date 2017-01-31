#' pcadapt objects
#'
#' \code{create.pcadapt} computes the numerical quantities needed to compute the test
#' statistics, and stores them in an object of class \code{pcadapt}.
#'
#' @param input a genotype matrix or a character string specifying the name of the file to be processed with \code{pcadapt}.
#' @param K an integer specifying the number of principal components to retain.
#' @param method a character string specifying the method to be used to compute
#' the p-values. Three statistics are currently available, \code{"mahalanobis"},
#' \code{"communality"}, and \code{"componentwise"}.
#' @param min.maf a value between \code{0} and \code{0.45} specifying the threshold
#' of minor allele frequencies above which p-values are computed.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' @param type an integer specifying the input type.
#' 
#' @return The returned value is a list containing all the numerical quantities needed to create an object of class \code{pcadapt}.
#' 
#' @importFrom RSpectra eigs_sym
#' 
#' @export
#'
create.pcadapt = function(input, K, method, min.maf, ploidy, type){
  if (type == 0){
    unused_matrix <- matrix(0, nrow = 2, ncol = 2)
    aux <- cmpt_cov_cpp(filename = input, xmatrix = unused_matrix, min_maf = min.maf, ploidy = ploidy, type = 0)
    nIND <- aux$nIND
    nSNP <- aux$nSNP
    xcov <- aux$xcov
    xsvd <- RSpectra::eigs_sym(xcov, k = K)
    lr <- lrfunc_cpp(filename = input,
                     xmatrix = unused_matrix,
                     scores = xsvd$vectors, 
                     nIND = nIND,
                     nSNP = nSNP, 
                     K = K, 
                     ploidy = ploidy,
                     min_maf = min.maf, 
                     sigma = xsvd$values,
                     type = 0)
  } else if (type == 1){
    unused_string <- "null" 
    aux <- cmpt_cov_cpp(filename = unused_string, xmatrix = as.matrix(input), min_maf = min.maf, ploidy = ploidy, type = 1)
    nIND <- aux$nIND
    nSNP <- aux$nSNP
    xcov <- aux$xcov
    xsvd <- RSpectra::eigs_sym(xcov, k = K)
    lr <- lrfunc_cpp(filename = unused_string,
                     xmatrix = as.matrix(input),
                     scores = xsvd$vectors, 
                     nIND = nIND,
                     nSNP = nSNP, 
                     K = K, 
                     ploidy = ploidy,
                     min_maf = min.maf, 
                     sigma = xsvd$values,
                     type = 1)
  }
  obj.stat <- cmpt.stat(x = lr$zscores, 
                        sqrt(xsvd$values / (nIND - 1)),
                        K = K,
                        method = method,
                        nSNP = nSNP,
                        maf = lr$maf,
                        min.maf = min.maf)
  return(list(scores = xsvd$vectors, 
              singular.values = sqrt(xsvd$values / (nIND - 1)), 
              zscores = lr$zscores, 
              loadings = lr$loadings,
              maf = lr$maf, 
              missing = lr$missing,
              stat = obj.stat$stat,
              gif = obj.stat$gif,
              chi2.stat = obj.stat$chi2.stat,
              pvalues = obj.stat$pvalues))
}
