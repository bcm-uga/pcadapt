#' Principal Component Analysis based on the correlation matrix
#'
#' \code{corpca} is an auxiliary function that performs principal components analysis on a dataset. It returns an object \code{x}
#' which contains the loadings, the scores and the singular values of the \code{K} first principal components.
#' It handles missing values in a dataset and actually computes the eigen elements of the \code{n x n}
#' covariance matrix, where \code{n} is the number of individuals.
#'
#' @param data a data matrix or a data frame.
#' @param K an integer specifying the number of principal components that are retained.
#'
#' @importFrom stats cov
#'
#' @examples
#' x <- NULL
#'
#' @keywords internal
#'
#' @export
corpca = function(data,K){
  n <- dim(data)[1]
  p <- dim(data)[2]
  cat(paste0("Number of SNPs: ",p,"\n")) 
  cat(paste0("Number of populations: ",n,"\n")) 
  data_aux <- scale(data,scale=FALSE)*sqrt(p/(n-1))
  covmat <- cov(t(data_aux),use="pairwise.complete.obs")
  res <- NULL
  aux <- eigen(covmat,symmetric=TRUE)
  sdev <- aux$values[1:K]
  res$scores <- aux$vectors[,1:K]
  aux_ldgs <- t(aux$vectors)%*%data_aux
  res$loadings <- array(0,dim=c(p,K))
  res$loadings[,] <- t((1/(sqrt(sdev)))*aux_ldgs[1:K,])
  res$singular.values <- sqrt(abs(sdev))
  return(res)
}