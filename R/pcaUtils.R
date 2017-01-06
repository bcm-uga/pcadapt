#' Compute covariance matrix for large genotype data
#'
#' @param input.filename a character string specifying the name of the file to be processed with \code{pcadapt}.
#' @return The returned value is a square symmetric matrix.
#'
#' @useDynLib pcadapt compute_covariance
#'
#' @export
#'
compute.covariance <- function(input.filename){
  input.size <- getsize(input.filename)
  nIND <- input.size[1]
  res <- .C("compute_covariance",as.character(input.filename),0.05,0,"tmp.pcadapt",PACKAGE = "pcadapt",result=as.double(array(0,dim=nIND*nIND)));
  return(matrix(res$result,nrow = nIND,ncol=nIND));
}

#' Compute covariance matrix for loaded genotype data
#'
#' @param x a genotype matrix.
#' @param ploidy an integer.
#' @return The returned value is a square symmetric matrix.
#'
#' @export
#'
cov.std <- function(x,ploidy=2){
  n <- ncol(x)
  p <- nrow(x)
  af <- apply(x, FUN = function(x){mean(x, na.rm = TRUE)}, MARGIN = 1) / ploidy
  dts <- scale(t(x), scale = sqrt(ploidy * (af * (1 - af))))*sqrt(p / (n-1))
  res <- cov(t(dts), use = "pairwise.complete.obs") * (n - 1)
  return(res)
}

#' Return the number of rows and the the number of columns of large genotype matrices
#'
#' @param input.filename a character string specifying the name of the file to be processed with \code{pcadapt}.
#' @return The returned value is a list of two integers.
#'
#' @useDynLib pcadapt get_size
#'
#' @export
#'
getsize <- function(input.filename){
  res <- .C("get_size",as.character(input.filename),size=as.integer(numeric(2)));
  return(res$size);
}

#' SVD for large genotype matrices
#'
#' @param input.filename a character string specifying the name of the file to be processed with \code{pcadapt}.
#' @return The returned value is an object containing eigenvectors and eigenvalues.
#'
#' @importFrom RSpectra eigs_sym 
#'
#' @export
#'
geteigen <- function(input.filename,K){
  cov.mat <- compute.covariance(input.filename)
  res <- RSpectra::eigs_sym(cov.mat,k=K)   
  return(list(scores=res$vectors,singular.values=sqrt(res$values/(ncol(cov.mat)-1))))
}

