#' Redundant Analysis for outlier detection
#'
#' \code{rdadapt}
#'
#' @param input a genotype matrix or a character string specifying the name of 
#' the file to be processed with \code{pcadapt}.
#' @param env a matrix of environmental or phenotypic variables.
#' 
#' @return The returned value \code{x} is an object of class \code{pcadapt}.
#' 
#' @useDynLib pcadapt
#'
#' @export
#'
rdadapt = function(input, env) {
  nSNP <- nrow(input)
  maf <- cmpt_minor_af(input, 2)
  normalized.input <- scale(t(input), 
                            center = TRUE, 
                            scale = sqrt(2 * maf * (1 - maf)))
  normalized.env <- scale(env, center = TRUE, scale = TRUE)
  svd.X <- svd(normalized.env)
  Yhat <- get_fitted_matrix(normalized.input, svd.X$u)
  obj.rda <- RSpectra::svds(Yhat, k = ncol(env))
  obj.rda$res <- normalized.input - Yhat
  return(obj.rda)
}