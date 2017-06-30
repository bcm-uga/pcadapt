#' @export
#' 
cmpt_stat_bary = function(scaled.geno, V, sigma, window_size, pop, adm, donor, anc, PCs){
  nSNP = ncol(scaled.geno)
  nIND = nrow(scaled.geno)
  map <- 1:nSNP
  hws = window_size / 2
  idx_old = get_window(hws, map, window_size, 0)    
  loop_beg = hws + 1
  loop_end = nSNP - hws
  uloc = cmpt_local_pca(scaled.geno, V, sigma, idx_old[1], idx_old[2]);
  n.pop <- length(unique(pop))
  stat <- matrix(NA, nrow = nSNP, ncol = n.pop - 1)
  K <- max(1, n.pop - 2)
  tmp <- matrix(0, nrow = nIND, ncol = K)
  for (i in loop_beg:loop_end){
    idx_new = get_window(i, map, window_size, 0);
    print(idx_new)
    updt_local_scores(uloc, scaled.geno, V, sigma, idx_old[1], idx_old[2], idx_new[1], idx_new[2]);
    tmp[, 1:K] <- uloc[, 1:K]
    stat[i, ] <- bary_to_anc(scores = tmp, pop = pop, adm = adm)
    idx_old <- idx_new;
  }
  return(stat)
}

#' @export
#'
scores_centroids = function(scores, pop){
  pop.names <- unique(pop)
  n.pop <- length(pop.names) 
  K <- ncol(scores)
  centroids <- matrix(0, nrow = n.pop, ncol = K)
  if (K > 1){
    for (i in 1:n.pop){
      if (sum(pop == pop.names[i]) > 1){
        centroids[i, ] <- apply(scores[pop == pop.names[i], ], 2, mean)  
      } else if (sum(pop == pop.names[i]) == 1){
        centroids[i, ] <- scores[pop == pop.names[i], ]
      }
    }
  } else if (K == 1){
    for (i in 1:n.pop){
      centroids[i, ] <- mean(scores[pop == pop.names[i], 1])
    }  
  }
  return(centroids)
}


#' @export
#'
centroids_to_simplex = function(centroids, pop, adm){
  pop.names <- unique(pop)
  K <- ncol(centroids)
  return(matrix(centroids[pop.names != adm, ], ncol = K))
}

#' @export
#' 
bary_to_anc = function(scores, pop, adm){
  cent <- scores_centroids(scores, pop)
  simp <- centroids_to_simplex(cent, pop, adm)
  bary.coord <- cart2bary(simp, as.matrix(scores[pop == adm, ]))  
  anc <- apply(bary.coord, MARGIN = 2, mean)
  return(anc)
}

#' @export
#' 
bary_to_coord = function(scores, pop, adm){
  cent <- scores_centroids(scores, pop)
  simp <- centroids_to_simplex(cent, pop, adm)
  bary.coord <- cart2bary(simp, as.matrix(scores[pop == adm, ]))  
  return(bary.coord)
}

#' @export
#' 
scalar_prod = function(scores, anc, donor, adm, PCs){
  cent <- scores_centroids(scores, pop)
  pop.name <- unique(pop)
  c1 <- cent[pop.name == anc, PCs]
  c2 <- cent[pop.name == donor, PCs]
  c.adm <- cent[pop.name == adm, PCs]
  res <- sum((c.adm - c1) * (c2 - c1)) / (sum((c.adm - c1) * (c.adm - c1)) * sum((c2 - c1) * (c2 - c1))) 
  return(res)
}

#' @export
#'
scan.intro.2 = function(input, 
                        K = 2, 
                        pop, 
                        admxd,
                        min.maf = 0.05, 
                        ploidy = 2, 
                        window.size = 1000
){
  
  if (is.character(input)){
    geno <- as.matrix(data.table::fread(input))
  } else if (class(input) %in% c("array", "matrix", "data.frame")){
    geno <- as.matrix(input)
  } else {
    stop("Wrong argument.")
  }
  
  nSNP <- nrow(geno)
  
  if (length(K) == 1 && K == 1){
    k = 2
  } else {
    k = max(K)
  }
  
  maf <- cmpt_minor_af(xmatrix = geno, ploidy = ploidy)
  xaxis <- (1:nSNP)[maf >= min.maf]
  gmap <- (1:nSNP)[maf >= min.maf]
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
  
  s_1 <- cmpt_stat_bary(scaled.geno = as.matrix(scaled.geno), 
                        V = as.matrix(obj.svd$v), 
                        sigma = as.vector(obj.svd$d), 
                        window_size = window.size,  
                        pop = pop, 
                        adm = admxd
  )    
  
  complete.stat <- rep(NA, length = nSNP)
  complete.stat[maf >= min.maf] <- s_1[ ,1]
  m <- median(complete.stat, na.rm = TRUE)
  s <- mad(complete.stat, na.rm = TRUE)
  complete.stat[1] <- m
  complete.stat[nSNP] <- m
  yint <- approx((1:nSNP)[!is.na(complete.stat)], complete.stat[!is.na(complete.stat)], 1:nSNP)  
  median.y <- stats::median(yint$y, na.rm = TRUE)
  mad.y <- stats::mad(yint$y, na.rm = TRUE)
  obj.stat <- list()
  obj.stat[[1]] <- (yint$y - median.y) / mad.y
  obj.stat$scores <- obj.svd$u
  # obj.stat[[2]] <- gmap
  
  flush.console()
  cat("DONE\n")
  
  return(obj.stat)
} 

#' @export
#' 
cmpt_stat_scalar = function(scaled.geno, V, sigma, window_size, pop, adm, donor, anc, PCs){
  nSNP = ncol(scaled.geno)
  nIND = nrow(scaled.geno)
  map <- 1:nSNP
  hws = window_size / 2
  idx_old = get_window(hws, map, window_size, 0)    
  loop_beg = hws + 1
  loop_end = nSNP - hws
  uloc = cmpt_local_pca(scaled.geno, V, sigma, idx_old[1], idx_old[2]);
  n.pop <- length(unique(pop))
  #stat <- matrix(NA, nrow = nSNP, ncol = n.pop - 1)
  stat <- rep(NA, length = nSNP)
  #K <- max(1, n.pop - 2)
  #tmp <- matrix(0, nrow = nIND, ncol = K)
  for (i in loop_beg:loop_end){
    idx_new = get_window(i, map, window_size, 0);
    print(idx_new)
    updt_local_scores(uloc, scaled.geno, V, sigma, idx_old[1], idx_old[2], idx_new[1], idx_new[2]);
    #tmp[, 1:K] <- uloc[, 1:K]
    #stat[i, ] <- bary_to_anc(scores = tmp, pop = pop, adm = adm)
    stat[i] <- scalar_prod(uloc, anc, donor, adm, PCs)
    idx_old <- idx_new;
  }
  return(stat)
}
