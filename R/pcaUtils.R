#' Covariance for loaded genotype data 
#'
#' @param input a genotype matrix.
#' @param K an integer specifying the number of principal components to retain.
#' 
#' @return The returned value is the covariance matrix.
#'
#' @export
#'
cmpt.cov.matrix <- function(input, ploidy = 2){
  n <- ncol(input)
  p <- nrow(input)
  af <- apply(input, FUN = function(x){mean(x, na.rm = TRUE)}, MARGIN = 1) / ploidy
  dts <- scale(t(input), scale = sqrt(ploidy * (af * (1 - af)))) * sqrt(p / (n - 1))
  res <- tAA_cpp(t(dts), nrow = ncol(dts), ncol = nrow(dts))
  return(res)
}



