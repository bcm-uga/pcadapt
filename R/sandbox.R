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

#' LD clumping
#'
#' \code{clumping} is adapted from the snp_clumping function implemented in 
#' the bigsnpr package developed by Florian Prive.
#'
#' @param G a genotype matrix. 
#' @param chr.info a vector containing the chromosome information for each 
#' marker.
#' @param size an integer.
#' @param thr a numerical value.
#' 
#' @return The returned value is a logical vector.
#' 
#' @importFrom stats var
#'
#' @export
#'
clumping3 = function(G, chr.info, size = 100, thr = 0.2) {
  n <- ncol(G) # Number of individuals
  p <- nrow(G) # Number of SNPs
  S <- cmpt_minor_af(G, 2)
  sumX <- apply(G, MARGIN = 1, FUN = function(h) {sum(h, na.rm = TRUE)})
  denoX <- apply(G, MARGIN = 1, FUN = function(h) {var(h, na.rm = TRUE)})
  denoX <- (n - 1) * denoX # Does not account for missing values -> (n - na - 1)
  if (missing(chr.info)) {
    ord.chr <- order(S, decreasing = TRUE)
    remain <- rep(TRUE, length(ord.chr))
    ind.keep <- clumping_cpp(G,
                             ord.chr,
                             remain,
                             sumX,
                             denoX,
                             size, 
                             thr)
  } else {
    ind.chrs <- split(seq_along(chr.info), chr.info)
    ind.keep <- NULL
    for (n.chr in 1:length(ind.chrs)) {
      ord.chr <- order(S[ind.chrs[[n.chr]]], decreasing = TRUE)
      remain <- rep(TRUE, length(ord.chr))
      aux <- clumping_cpp(G[ind.chrs[[n.chr]], ],
                          ord.chr,
                          remain,
                          sumX[ind.chrs[[n.chr]]],
                          denoX[ind.chrs[[n.chr]]],
                          size, 
                          thr)
      ind.keep <- c(ind.keep, aux)
    }
  }
  return(ind.keep)
}

clumping2 = function(input, size = 100, thr = 0.2) {
  
  lookup_byte <- getCode()
  
  if (class(input) == "character") {
    path_to_bed <- normalizePath(input)
    p <- nrow(data.table::fread(sub("\\.bed$", ".bim", path_to_bed)))
    n <- nrow(data.table::fread(sub("\\.bed$", ".fam", path_to_bed)))
    
    ### File mapping
    xptr <- bedXPtr(path_to_bed, n, p)
    
  } else if (class(input) == "matrix") {
    # an input matrix has nIND rows and nSNP columns
    xptr <- input
    n <- nrow(xptr)
    p <- ncol(xptr)
  }
  
  tmp <- pcadapt:::af(xptr, rbind(rep(0, p), 1, 2, 3), lookup_byte)
  S <- pmin(tmp, 1 - tmp)
  sumX <- get_sumX(xptr, rbind(rep(0, p), 1, 2, 3), lookup_byte)
  denoX <- get_denoX(xptr, rbind(rep(0, p), 1, 2, 3), lookup_byte, tmp)
  
  ord.chr <- order(S, decreasing = TRUE)
  remain <- rep(TRUE, length(ord.chr))
  ind.keep <- clumping(xptr,
                       rbind(rep(0, p), 1, 2, 3),
                       lookup_byte,
                       ord.chr,
                       remain,
                       sumX,
                       denoX,
                       size, 
                       thr)
  
  return(ind.keep)
}

