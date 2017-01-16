#' pcadapt objects
#'
#' \code{create.pcadapt.file} computes the numerical quantities needed to compute the test
#' statistics, and stores them in an object of class \code{pcadapt}.
#'
#' @param input a character string specifying the name of the file to be processed with \code{pcadapt}.
#' @param K an integer specifying the number of principal components to retain.
#' @param method a character string specifying the method to be used to compute
#' the p-values. Three statistics are currently available, \code{"mahalanobis"},
#' \code{"communality"}, and \code{"componentwise"}.
#' @param min.maf a value between \code{0} and \code{0.45} specifying the threshold
#' of minor allele frequencies above which p-values are computed.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' 
#' @return The returned value is a list containing all the numerical quantities needed to create an object of class \code{pcadapt}.
#' 
#' @importFrom RSpectra eigs_sym
#' 
#' @export
#'
create.pcadapt.file = function(input, K, method, min.maf, ploidy){
  aux <- cmpt_cov_file(path = input, min_maf = min.maf, ploidy = ploidy)  
  nIND <- aux$nIND
  nSNP <- aux$nSNP
  xcov <- aux$xcov
  xsvd <- RSpectra::eigs_sym(xcov, k = K)
  lr <- lrfunc_file(filename = input, 
                    scores = xsvd$vectors, 
                    nIND = nIND,
                    nSNP = nSNP, 
                    K = K, 
                    ploidy = ploidy,
                    min_maf = min.maf)
  obj.stat <- cmpt.stat(x = lr$zscores, 
                        s.v = sqrt(xsvd$values / (nIND - 1)),
                        K = K,
                        method = method,
                        nSNP = nSNP,
                        maf = lr$maf,
                        min.maf = min.maf)
  return(list(scores = xsvd$vectors, 
              singular.values = sqrt(xsvd$values / (nIND - 1)), 
              zscores = lr$zscores, 
              maf = lr$maf, 
              missing = lr$missing,
              stat = obj.stat$stat,
              gif = obj.stat$gif,
              chi2.stat = obj.stat$chi2.stat,
              pvalues = obj.stat$pvalues))
}

#' pcadapt objects
#'
#' \code{create.pcadapt.matrix} computes the numerical quantities needed to compute the test
#' statistics, and stores them in an object of class \code{pcadapt}.
#'
#' @param input a genotype matrix.
#' @param K an integer specifying the number of principal components to retain.
#' @param method a character string specifying the method to be used to compute
#' the p-values. Three statistics are currently available, \code{"mahalanobis"},
#' \code{"communality"}, and \code{"componentwise"}.
#' @param min.maf a value between \code{0} and \code{0.45} specifying the threshold
#' of minor allele frequencies above which p-values are computed.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' 
#' @return The returned value is a list containing all the numerical quantities needed to create an object of class \code{pcadapt}.
#' 
#' @importFrom RSpectra eigs_sym
#' 
#' @export
#'
create.pcadapt.matrix = function(input, K, method, min.maf, ploidy){
  aux <- cmpt_cov_matrix(as.matrix(input), min_maf = min.maf, ploidy = ploidy)
  nIND <- aux$nIND
  nSNP <- aux$nSNP
  xcov <- aux$xcov
  xsvd <- RSpectra::eigs_sym(xcov, k = K)
  lr <- lrfunc_matrix(Geno = input, 
                      scores = xsvd$vectors,
                      nIND = nIND,
                      nSNP = nSNP,
                      K = K,
                      ploidy = ploidy,
                      min_maf = min.maf)
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
              maf = lr$maf, 
              missing = lr$missing,
              stat = obj.stat$stat,
              gif = obj.stat$gif,
              chi2.stat = obj.stat$chi2.stat,
              pvalues = obj.stat$pvalues))
}





