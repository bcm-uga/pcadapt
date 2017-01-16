#' Convert genotypes to pooled samples
#'
#' \code{get.pool.matrix} creates a pooled-sequenced data out of a genotype matrix,
#' given the labels of each individuals.
#' 
#' @param data a matrix with n rows and p columns where n is the number of individuals and p is the number of markers. 
#' @param pop a list of integers or strings specifying which subpopulation the individuals belong to.
#' @param ploidy an integer specifying the ploidy of the individuals.
#'
#' @export
#' 
get.pool.matrix = function(data, pop, ploidy = 2){
  nSNP <- ncol(data)
  pop.names <- unique(pop)
  freq <- array(0, dim = c(length(pop.names), nSNP))
  for (k in 1:length(pop.names)){
    geno.k <- data[pop == pop.names[k], ]
    geno.k[geno.k == 9] <- NA
    freq[k, ] <- apply(geno.k, MARGIN = 2, FUN = function(h){mean(h, na.rm = TRUE)}) / ploidy
  }
  return(freq)
}