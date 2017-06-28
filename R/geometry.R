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
  #res <- sum((c.adm - c1) * (c2 - c1)) / (sum((c.adm - c1) * (c.adm - c1)) * sum((c2 - c1) * (c2 - c1))) 
  res <- sum((c.adm - c1) * (c2 - c1)) / sqrt((sum((c2 - c1) * (c2 - c1)))) 
  return(res)
}

