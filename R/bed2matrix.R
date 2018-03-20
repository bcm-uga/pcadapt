#' Convert a bed to a matrix
#'
#' @param bedfile Path to a bed file. 
#' @param n Number of samples. Default reads it from coresponding fam file.
#' @param p Number of SNPs. Default reads it from coresponding bim file.
#'
#' @return An integer matrix.
#' @export
#'
#' @examples
#' bedfile <- system.file("extdata", "geno3pops.bed", package = "pcadapt")
#' mat <- bed2matrix(bedfile)
#' dim(mat)
#' table(mat)
bed2matrix <- function(bedfile, n = NULL, p = NULL) {
  
  if (is.null(n) || is.null(p)) {
    bed <- read.pcadapt(bedfile, type = "bed")
    n <- attr(bed, "n")
    p <- attr(bed, "p")
  } 
  
  xptr <- bedXPtr(bedfile, n, p)
  mat <- bed2mat(xptr)
  mat[is.na(mat)] <- NA
  mat
}
