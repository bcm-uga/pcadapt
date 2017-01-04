#' Regression 
#'
#' @param mat a matrix.
#' @param filename a character string.
#' @param nind an integer.
#' @param nloci an integer.
#' @param K an integer.
#' @param ploidy an integer.
#' @param min.maf a real number.
#' 
#' @useDynLib pcadapt lrfunc_
#'
#' @export
#'
lrfunc <- function(matrix, filename, nind, nloci, K, ploidy, min.maf){
  .C("lrfunc_", as.double(matrix), as.character(filename), as.integer(nind), as.integer(nloci), as.integer(K), as.integer(ploidy), as.double(min.maf));
}