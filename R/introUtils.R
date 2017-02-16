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


#' Convert string labels to integer labels
#'
#' \code{assign.int.labels} returns a vector of integers.
#'
#' @param pop a vector of integers or strings specifying which subpopulation the individuals belong to.
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
#' \code{scan.intro} computes statistics to detect excesses of local ancestry based on a PCA approach.
#'
#' @param input a genotype matrix or a character string specifying the name of the file to be imputed.
#' @param K a vector of integers specifying the components along which local ancestries may vary.
#' @param pop a vector of integers or strings specifying which subpopulation the individuals belong to.
#' @param ancstrl.1 a string specifying the label of the ancestral population genetically closer to the hybrid population.
#' @param ancstrl.2 a string specifying the label of the ancestral population genetically further from the hybrid population.
#' @param admxd a string specifying the label of thehybrid population.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' @param window.size an integer specifying the window size.

#' @param impute a logical value indicating whether the input data has to imputed. 
#' 
#' @return The returned value is a list containing the test statistics and the associated p-values.
#' 
#' @export
#'
scan.intro = function(input, 
                      K = 2, 
                      pop, 
                      ancstrl.1,
                      ancstrl.2,
                      admxd,
                      min.maf = 0.05, 
                      ploidy = 2, 
                      window.size = 1000, 
  
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
  
  pop.names <- get.pop.names(pop)
  pop.int <- assign.int.labels(pop)
  ancstrl.int.1 <- which(pop.names == ancstrl.1)
  ancstrl.int.2 <- which(pop.names == ancstrl.2)
  admxd.int <- which(pop.names == admxd)
  
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
  cat("Computing the statistics...")
  
  stat <- cmpt_all_stat(geno = scaled.geno, 
                        V = obj.svd$v, 
                        sigma = obj.svd$d, 
                        window_size = window.size,  
                        direction = 0, 
                        lab = pop.int, 
                        ancstrl1 = ancstrl.int.1,
                        ancstrl2 = ancstrl.int.2,
                        adm = admxd.int, 
                        axis = axis.vector)
  stat.sd <- sd(stat)
  pval <- pnorm(stat / stat.sd, lower.tail = FALSE)
  flush.console()
  cat("DONE\n")
  return(-log10(pval))
} 

draw.local.pca = function(geno, V, sigma, uglob, beg, end, pop, i = 1, j = 2, ancstrl1, ancstrl2, adm){
  uloc <- cmpt_local_pca(geno, V, sigma = sigma, beg = beg, end = end)
  dloc <- vector(length = ncol(V), mode = "numeric")
  dglob <- vector(length = ncol(V), mode = "numeric")
  s <- vector(length = ncol(V), mode = "numeric")
  R <- matrix(0, nrow = ncol(V), ncol = ncol(V))
  cent <- cmpt_centroids(uglob, pop, ancstrl1, ancstrl2)
  cent.loc <- cmpt_centroids(uloc, pop, ancstrl1, ancstrl2)
  axis <- cent$m2 - cent$m2
  shape.1 <- as.matrix(t(cbind(cent$m1, cent$m2)))
  
  shape.2 <- as.matrix(t(cbind(cent.loc$m1, cent.loc$m2)))
  
  # R <- pca_rotation(shape.1, shape.2)
  # R <- matrix(0, nrow = ncol(V), ncol = ncol(V))
  # diag(R) <- 1
  
  cmpt_transformation(uloc, uglob, lab, ancstrl1, ancstrl2, s, dloc, dglob, R);
  usc <- rescale_local_pca(uloc, s, dloc, dglob, R);
  
  xmin <- min(min(uglob[, i]), min(usc[, i]))
  xmax <- max(max(uglob[, i]), max(usc[, i]))
  ymin <- min(min(uglob[, j]), min(usc[, j]))
  ymax <- max(max(uglob[, j]), max(usc[, j]))
  
  plot(usc[, i], usc[, j], col = as.factor(pop), pch = 19, cex = 0.5,  xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  points(uglob[, i], uglob[, j], col = as.factor(pop), cex = 1)
  arrows(cent$m1[1], cent$m1[2], cent$m2[1], cent$m2[2])
  arrows(uglob[pop == adm, i], uglob[pop == adm, j], usc[pop == adm, i], usc[pop == adm, j])
}