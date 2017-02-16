#' Genotype matrix imputation
#'
#' \code{impute.pcadapt} imputes values based on medians.
#'
#' @param input a genotype matrix or a character string specifying the name of the file to be imputed.
#' @param pop a vector of integers or strings specifying which subpopulation the individuals belong to.
#' @param skip.return a logical value specifying whether the list of markers to be skipped should be returned or not. 
#' 
#' @return The returned value is a list containing the test statistics and the associated p-values.
#' 
#' @export
#'
impute.pcadapt = function(input, pop, skip.return = FALSE){
  if (is.character(input)){
    dt <- as.matrix(data.table::fread(input))
  } else if (class(input) %in% c("array", "matrix", "data.frame")){
    dt <- as.matrix(input)
  } else {
    stop("Wrong argument.")
  }
  
  if (missing(pop)){
    y <- impute_geno(dt)   
  } else if (!missing(pop)){
    pop.names <- pcadapt::get.pop.names(pop)
    y <- impute_geno_pop(dt, pop, pop.names)   
  }
  if (skip.return == FALSE){
    return(list(x = y$x[y$skip == 0, ]))  
  } else if (skip.return == TRUE){
    return(list(x = y$x[y$skip == 0, ], ix = which(y$skip == 1))) 
  }
}

scan.intro = function(input, 
                      K = 2, 
                      pop, 
                      min.maf = 0.05, 
                      ploidy = 2, 
                      window.size = 1000, 
                      ancstrl.1,
                      ancstrl.2,
                      admxd,
                      impute = FALSE){
  if (impute){
    geno <- (impute.pcadapt(input = input, pop = pop))$x
  } else if (!impute){
    if (is.character(input)){
      geno <- as.matrix(data.table::fread(input))
    } else if (class(input) %in% c("array", "matrix", "data.frame")){
      geno <- as.matrix(input)
    } else {
      stop("Wrong argument.")
    }
  }
  
  if (length(K) == 1 && K == 1){
    k = 2
  } else {
    k = max(K)
  }
  
  axis.vector <- vector(length = k, mode = "numeric")
  axis.vector[K] <- 1
  
  maf <- cmpt_minor_af(xmatrix = geno, ploidy = ploidy)
  geno <- geno[maf >= min.maf, ]
  maf <- maf[maf >= min.maf]
  sd <- sqrt(ploidy * maf * (1 - maf))
  cat("Scaling the genotype matrix...")
  scaled.geno <- scale(t(geno), center = TRUE, scale = sd) 
  cat("DONE\n")
  cat("Performing PCA...\n")
  obj.svd <- svd.pcadapt(input = geno, K = k, min.maf = min.maf, ploidy = ploidy, type = 1)
  cat("DONE\n")
  cat("Computing the statistics...")
  
  stat <- cmpt_all_stat(geno = scaled.geno, 
                        V = obj.svd$v, 
                        sigma = obj.svd$d, 
                        window_size = window.size,  
                        direction = 0, 
                        lab = pop, 
                        ancstrl1 = ancstrl.1,
                        ancstrl2 = ancstrl.2,
                        adm = admxd, 
                        axis = axis.vector)
  stat.sd <- sd(stat)
  pval <- pnorm(stat / stat.sd, lower.tail = FALSE)
  flush.console()
  cat("DONE\n")
  return(-log10(pval))
} 
