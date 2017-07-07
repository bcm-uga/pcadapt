#' Singular Value Decomposition
#'
#' \code{svd.pcadapt} computes the numerical quantities needed to compute the 
#' test statistics, and stores them in an object of class \code{pcadapt}.
#'
#' @param input a genotype matrix or a character string specifying the name of 
#' the file to be processed with \code{pcadapt}.
#' @param K an integer specifying the number of principal components to retain.
#' @param min.maf a value between \code{0} and \code{0.45} specifying the 
#' threshold of minor allele frequencies above which p-values are computed.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' @param type an integer specifying the input type: \code{0} for genotype 
#' matrices or \code{1} for genotype files.
#' 
#' @return The returned value is a list containing the Singular Value 
#' Decomposition of the centered and scaled genotype matrix.
#' 
#' @importFrom RSpectra eigs_sym
#' 
#' @export
#'
svd.pcadapt = function(input, K, min.maf = 0.05, ploidy = 2, type){
  if (type == 0){
    unused_matrix <- matrix(0, nrow = 2, ncol = 2)
    aux <- cmpt_cov_cpp(filename = input, xmatrix = unused_matrix, 
                        min_maf = min.maf, ploidy = ploidy, type = 0)
    nIND <- aux$nIND
    nSNP <- aux$nSNP
    xcov <- aux$xcov
    xsvd <- RSpectra::eigs_sym(xcov, k = K)
    ld <- cmpt_loadings(filename = input,
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
    aux <- cmpt_cov_cpp(filename = unused_string, xmatrix = as.matrix(input), 
                        min_maf = min.maf, ploidy = ploidy, type = 1)
    nIND <- aux$nIND
    nSNP <- aux$nSNP
    xcov <- aux$xcov
    xsvd <- RSpectra::eigs_sym(xcov, k = K)
    ld <- cmpt_loadings(filename = unused_string,
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
  return(list(u = xsvd$vectors, 
              d = sqrt(xsvd$values), 
              v = ld / sqrt(nSNP)))
}
