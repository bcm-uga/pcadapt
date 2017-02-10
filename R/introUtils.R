impute.pcadapt = function(input, lab, skip.return = FALSE){
  if (is.character(input)){
    dt <- as.matrix(data.table::fread(input))
  } else if (class(input) %in% c("array", "matrix", "data.frame")){
    dt <- as.matrix(input)
  } else {
    stop("Wrong argument.")
  }
  pop <- pcadapt::get.pop.names(lab)
  if (missing(lab)){
    y <- impute_geno(dt)   
  } else if (!missing(lab)){
    y <- impute_geno_pop(dt, lab, pop)   
  }
  if (skip.return == FALSE){
    return(list(x = y$x[y$skip == 0, ]))  
  } else if (skip.return == TRUE){
    return(list(x = y$x[y$skip == 0, ]), ix = which(y$skip == 1)) 
  }
}

scan.intro = function(input, 
                      K = 2, 
                      lab, 
                      min.maf = 0.05, 
                      ploidy = 2, 
                      window.size = 1000, 
                      direction = 1, 
                      ancstrl.1,
                      ancstrl.2,
                      admxd,
                      impute = FALSE){
  if (impute){
    geno <- (impute.pcadapt(input = input, lab = lab))$x
  } else if (!impute){
    if (is.character(input)){
      geno <- as.matrix(data.table::fread(input))
    } else if (class(input) %in% c("array", "matrix", "data.frame")){
      geno <- as.matrix(input)
    } else {
      stop("Wrong argument.")
    }
  }
  maf <- cmpt_minor_af(xmatrix = geno, ploidy = ploidy)
  geno <- geno[maf >= min.maf, ]
  maf <- maf[maf >= min.maf]
  sd <- sqrt(ploidy * maf * (1 - maf))
  cat("Scaling the genotype matrix...")
  scaled.geno <- scale(t(geno), center = TRUE, scale = sd) 
  cat("DONE\n")
  cat("Performing PCA...\n")
  obj.svd <- svd.pcadapt(input = geno, K = K, min.maf = min.maf, ploidy = ploidy, type = 1)
  flush.console()
  cat("DONE\n")
  cat("Computing the statistics...")
  stat <- cmpt_all_stat(geno = scaled.geno, 
                        V = obj.svd$v, 
                        sigma = obj.svd$d, 
                        window_size = window.size,  
                        direction = 0, 
                        lab = lab, 
                        ancstrl1 = ancstrl.1,
                        ancstrl2 = ancstrl.2,
                        adm = admxd, 
                        axis = 0)
  stat.sd <- sd(stat)
  pval <- pnorm(stat / stat.sd, lower.tail = FALSE)
  flush.console()
  cat("DONE\n")
  return(-log10(pval))
} 
