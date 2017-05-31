#' Genotype matrix imputation
#'
#' \code{impute.pcadapt} imputes values based on medians.
#'
#' @param input a genotype matrix or a character string specifying the name of 
#' the file to be imputed.
#' @param pop a vector of integers or strings specifying which subpopulation the 
#' individuals belong to.
#' @param skip.return a logical value specifying whether the list of markers to 
#' be skipped should be returned or not. 
#' 
#' @return The returned value is a list containing the test statistics and the 
#' associated p-values.
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
    pop.int <- assign.int.labels(pop)
    pop.names <- pcadapt::get.pop.names(pop.int)
    y <- impute_geno_pop(dt, pop.int, pop.names)   
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
#' @param ancstrl.1 a string specifying the label of the ancestral population 
#' genetically closer to the hybrid population.
#' @param ancstrl.2 a string specifying the label of the ancestral population 
#' genetically further from the hybrid population.
#' @param admxd a string specifying the label of the hybrid population.
#' @param min.maf a value between \code{0} and \code{0.45} specifying the 
#' threshold of minor allele frequencies above which p-values are computed.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' @param window.size an integer specifying the window size.
#' @param impute a logical value indicating whether the input data has to 
#' imputed. 
#' @param chr.info a list containing the chromosome information for each marker.
#' @param map a numeric vector containing the genetic positions.
#' @param side a character string specifying whether the window should be aligned on 
#' the left, middle or right.
#' 
#' @return The returned value is a list containing the test statistics and the 
#' associated p-values.
#' 
#' @importFrom data.table fread
#' @importFrom MASS cov.rob
#' @importFrom stats approx
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
                      impute = FALSE, 
                      chr.info,
                      map,
                      side = "middle"){
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
  nSNP <- nrow(geno)
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
  
  xaxis <- (1:nSNP)[maf >= min.maf]
  
  if (missing(map)){
    gmap <- (1:nSNP)[maf >= min.maf]
  } else {
    gmap <- map[maf >= min.maf]
  }
  
  geno <- geno[maf >= min.maf, ]
  maf <- maf[maf >= min.maf]
  
  sd <- sqrt(ploidy * maf * (1 - maf))
  cat("Scaling the genotype matrix...")
  scaled.geno <- scale(t(geno), center = TRUE, scale = sd) 
  cat("DONE\n")
  cat("Performing PCA...\n")
  obj.svd <- svd.pcadapt(input = geno, K = k, min.maf = min.maf, 
                         ploidy = ploidy, type = 1)
  cat("Computing the statistics...\n")
  
  if (missing(map)){
    with.map <- 0  
  } else {
    with.map <- 1
  }
  
  if (side == "left"){
    side.int <- -1
  } else if (side == "middle"){
    side.int <- 0
  } else if (side == "right"){
    side.int <- 1
  }
  
  if (missing(chr.info)){
    s_1 <- cmpt_stat_introgr(geno = as.matrix(scaled.geno), 
                          V = as.matrix(obj.svd$v), 
                          sigma = as.vector(obj.svd$d), 
                          window_size = as.integer(window.size),  
                          direction = as.integer(0), 
                          lab = as.vector(pop.int), 
                          ancstrl1 = as.integer(ancstrl.int.1),
                          ancstrl2 = as.integer(ancstrl.int.2),
                          adm = as.integer(admxd.int), 
                          axis = as.vector(axis.vector),
                          map = gmap,
                          with_map = with.map,
                          side = side.int)  
    yint <- approx(gmap[!is.na(s_1)], s_1[!is.na(s_1)], 1:nSNP)  
    aux <- MASS::cov.rob(yint$y)
    obj.stat <- list()
    obj.stat[[1]] <- (yint$y - aux$center[1]) / sqrt(aux$cov[1, 1])
    obj.stat[[2]] <- gmap
  } else {
    chr <- chr.info[maf >= min.maf]
    chr.it <- unique(chr)
    stat <- vector(mode = "numeric", length = length(chr))
    obj.stat <- list()
    for (k in chr.it){
      cat("Analyzing chromosome ", k, "\n")
      chr_k <- (chr == k)
      gmap_k = gmap[chr_k]
      scaled.geno_k <- as.matrix(scaled.geno[, chr_k])
      v_k <- as.matrix(obj.svd$v[chr_k, ])
      s_k <- cmpt_stat_introgr(geno = scaled.geno_k, 
                            V = v_k, 
                            sigma = as.vector(obj.svd$d), 
                            window_size = as.integer(window.size),  
                            direction = as.integer(0), 
                            lab = as.vector(pop.int), 
                            ancstrl1 = as.integer(ancstrl.int.1),
                            ancstrl2 = as.integer(ancstrl.int.2),
                            adm = as.integer(admxd.int), 
                            axis = as.vector(axis.vector),
                            map = gmap_k)   
      aux <- MASS::cov.rob(s_k)
      obj.stat[[2 * k - 1]] <- (s_k - aux$center[1]) / sqrt(aux$cov[1, 1])
      obj.stat[[2 * k]] <- 1:(length(s_k) - window.size)
    }
  }
  flush.console()
  cat("DONE\n")
  class(obj.stat) <- "pcadapt"
  attr(obj.stat, "K") <- K
  attr(obj.stat, "method") <- "introgression"
  attr(obj.stat, "min.maf") <- min.maf
  attr(obj.stat, "ancstrl.1") <- ancstrl.1
  attr(obj.stat, "ancstrl.2") <- ancstrl.2
  attr(obj.stat, "window.size") <- window.size
  return(obj.stat)
} 

#' Display local PCA
#'
#' \code{draw.pca} displays both local and global scores.
#'
#' @param geno a genotype matrix.
#' @param V a loading matrix.
#' @param sigma a vector of singular values.
#' @param uglob a matrix of global scores.
#' @param beg an integer specifying the first marker to be included.
#' @param end an integer specifying the first marker to be excluded.
#' @param pop a list of integers.
#' @param i an integer indicating onto which principal component the individuals
#' are projected when the "scores" option is chosen. Default value is set to 
#' \code{1}.
#' @param j an integer indicating onto which principal component the individuals
#' are projected when the "scores" option is chosen. Default value is set to 
#' \code{2}.
#' @param ancstrl1 an integer.
#' @param ancstrl2 an integer.
#' @param adm an integer.
#' 
#' @return The returned value is a list containing the test statistics and the 
#' associated p-values.
#' 
#' @importFrom graphics arrows plot points
#' @importFrom stats pnorm
#' @importFrom utils flush.console
#' 
#' @export
#'
draw.pca = function(geno, V, sigma, uglob, beg, end, pop, i = 1, j = 2, 
                    ancstrl1, ancstrl2, adm){
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
  
  #R <- pca_rotation(shape.1, shape.2)
  R <- matrix(0, nrow = ncol(V), ncol = ncol(V))
  diag(R) <- 1
  
  cmpt_transformation(uloc, uglob, pop, ancstrl1, ancstrl2, s, dloc, dglob, R);
  usc <- rescale_local_pca(uloc, s, dloc, dglob, R);
  
  xmin <- min(min(uglob[, i]), min(usc[, i]))
  xmax <- max(max(uglob[, i]), max(usc[, i]))
  ymin <- min(min(uglob[, j]), min(usc[, j]))
  ymax <- max(max(uglob[, j]), max(usc[, j]))
  
  plot(usc[, i], usc[, j], col = as.factor(pop), pch = 19, cex = 0.5,  
       xlim = c(xmin, xmax), 
       ylim = c(ymin, ymax))
  points(uglob[, i], uglob[, j], col = as.factor(pop), cex = 1)
  #arrows(cent$m1[1], cent$m1[2], cent$m2[1], cent$m2[2])
  #arrows(uglob[pop == adm, i], uglob[pop == adm, j], 
  #       usc[pop == adm, i], usc[pop == adm, j])
}