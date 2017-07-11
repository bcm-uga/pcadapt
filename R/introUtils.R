#' Convert string labels to integer labels
#'
#' \code{assign.int.labels} returns a vector of integers.
#'
#' @param pop a vector of integers or strings specifying which subpopulation the
#' individuals belong to.
#' 
#' @return The returned value is a vector of integers.
#' 
#' @export
#'
assign.int.labels = function(pop){
  res <- vector(length = length(pop), mode = "numeric")
  tmp <- get.pop.names(pop)
  for (i in 1:length(tmp)){
    res[pop == tmp[i]] <- i
  }
  return(res)
}

#' Introgression
#'
#' \code{scan.intro} computes statistics to detect excesses of local ancestry 
#' based on a PCA approach.
#'
#' @param input a genotype matrix or a character string specifying the name of 
#' the file to be imputed.
#' @param K a vector of integers specifying the components along which local 
#' ancestries may vary.
#' @param pop a vector of integers or strings specifying which subpopulation the
#' individuals belong to.
#' @param ancestral a string specifying the label of the donor population.
#' @param admixed a string specifying the label of the hybrid population.
#' @param min.maf a value between \code{0} and \code{0.45} specifying the 
#' threshold of minor allele frequencies above which p-values are computed.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' @param window.size an integer specifying the window size.
#' @param impute a logical value indicating whether the input data has to 
#' imputed. 
#' @param chr.info a list containing the chromosome information for each marker.
#' @param map a numeric vector containing the genetic positions.
#' 
#' @return The returned value is a list containing the test statistics.
#' 
#' @importFrom data.table fread
#' @importFrom MASS cov.rob
#' @importFrom stats approx median mad
#' @importFrom utils flush.console
#' 
#' @export
#'
scan.intro = function(input, 
                      K = 2, 
                      pop, 
                      ancestral,
                      admixed,
                      min.maf = 0.05, 
                      ploidy = 2, 
                      window.size = 1000, 
                      impute = FALSE, 
                      chr.info,
                      map){
  if (impute){
    geno <- (impute.pcadapt(input = input, pop = pop))$x
  } else if (!impute){
    if (is.character(input)){
      geno <- as.matrix(data.table::fread(input, data.table = FALSE))
    } else if (class(input) %in% c("array", "matrix", "data.frame")){
      geno <- as.matrix(input)
    } else {
      stop("Wrong argument.")
    }
  }
  nSNP <- nrow(geno)
  
  if (length(K) == 1 && K == 1){
    k = 2
  } else {
    k = max(K)
  }
  
  maf <- cmpt_minor_af(xmatrix = geno, ploidy = ploidy)
  
  if (missing(map)){
    gmap <- (1:nSNP)[maf >= min.maf]
  } else {
    gmap <- map[maf >= min.maf]
  }
  
  geno <- geno[maf >= min.maf, ]
  filtered.maf <- maf[maf >= min.maf]
  
  sd <- sqrt(ploidy * filtered.maf * (1 - filtered.maf))
  cat("Scaling the genotype matrix...")
  scaled.geno <- scale(t(geno), center = TRUE, scale = sd) 
  cat("DONE\n")
  cat("Performing PCA...\n")
  obj.svd <- svd.pcadapt(input = geno, K = k, min.maf = min.maf, 
                         ploidy = ploidy, type = 1)
  cat("Computing the statistics...\n")
  
  stat <- slidingWindows(sgeno = as.matrix(scaled.geno),
                         d = as.vector(obj.svd$d),
                         v = as.matrix(obj.svd$v),
                         pop = pop,
                         popUnique = unique(pop),
                         admixed = admixed,
                         window_size = window.size,
                         gmap)
  
  stat.med <- apply(stat, MARGIN = 2, FUN = function(x){median(x, na.rm = TRUE)})
  obj.stat <- matrix(NA, nrow = nSNP, ncol = ncol(stat))
  obj.stat[maf >= min.maf, ] <- stat
  
  obj.stat[1, ] <- stat.med
  obj.stat[nSNP, ] <- stat.med
  
  for (k in 1:ncol(stat)){
    subset <- which(!is.na(obj.stat[, k]))
    y.int <- approx(subset, obj.stat[subset, k], 1:nSNP)    
    median.y <- stats::median(y.int$y, na.rm = TRUE)
    mad.y <- stats::mad(y.int$y, na.rm = TRUE)
    obj.stat[, k] <- (y.int$y - median.y) / mad.y
  }
  flush.console()
  cat("DONE\n")
  class(obj.stat) <- "pcadapt"
  attr(obj.stat, "method") <- "introgression"
  attr(obj.stat, "min.maf") <- min.maf
  attr(obj.stat, "window.size") <- window.size
  return(obj.stat)
} 