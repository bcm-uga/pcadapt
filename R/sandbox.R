#' Residuals-based statistics
#'
#' \code{residuals_to_stat}
#'
#' @param geno a scaled genotype matrix.
#' @param obj.svd an object with \code{u}, \code{d} and \code{v}.
#' @param K an integer.
#' @param pop a list of integers or strings specifying which subpopulation the 
#' individuals belong to.
#' @param admixed a string specifying the label of the hybrid population.
#' @param window.size an integer.
#' 
#' @return The returned value is a numeric vector.
#' 
#' @importFrom RSpectra eigs_sym
#' @importFrom RcppRoll roll_mean
#' @importFrom Matrix Diagonal
#' 
#' @export
#'
residuals_to_stat <- function(geno, obj.svd, K = 1, pop, admixed, window.size = 100){
  D <- Matrix::Diagonal(obj.svd$d)
  res <- geno - (cbind(obj.svd$u[, K] * obj.svd$d[ K])) %*% (t(obj.svd$v[, K]))
  mean.stat <- apply(res[pop == admixed, ],
                     MARGIN = 2,
                     FUN = function(X){mean(X, na.rm = TRUE)}
  )
  stat <- (mean.stat - median(mean.stat, na.rm = TRUE)) / mad(mean.stat, na.rm = TRUE)
  smooth.stat <- RcppRoll::roll_mean(stat, n = window.size, by = 1, align = "center")
  final.stat <- c(rep(smooth.stat[1], window.size / 2),
                  smooth.stat,
                  rep(tail(smooth.stat, n = 1), window.size / 2 - 1))  
  return(final.stat)
}
