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
get.pool.matrix = function(data,pop,ploidy=2){
  nSNP <- ncol(data)
  pop.names <- unique(pop)
  freq <- array(0,dim=c(length(pop.names),nSNP))
  for (k in 1:length(pop.names)){
    geno.k <- data[pop==pop.names[k],]
    geno.k[geno.k==9] <- NA
    freq[k,] <- apply(geno.k,MARGIN=2,FUN=function(xx){mean(xx,na.rm=TRUE)})/ploidy
  }
  return(freq)
}

#' Sample genotype matrix from pooled samples
#'
#' \code{sample.geno} samples a genotype matrix from pooled samples.
#' 
#' @param pool.matrix a matrix with n rows and p columns where n is the number of pools and is the number of markers.
#' @param ploidy an integer specifying the ploidy.
#' @param cover.matrix a matrix with n rows and p columns where n is the number of pools and is the number of markers.
#' @param pop.sizes a list specifying the number of individuals for each pool.
#' @param method a character string indicating the method used for sampling.
#'
#' @examples
#' ## see also ?pcadapt for examples
#'
#' @importFrom stats rbinom rbeta
#'
#' @export
#'
sample.geno = function(pool.matrix=NULL,ploidy=2,cover.matrix=NULL,pop.sizes=NULL,method="per.pop"){
  nPOOL <- nrow(pool.matrix)  
  nSNP <- ncol(pool.matrix)
  if (missing(pop.sizes) || is.null(pop.sizes)){
    sample.size <- rep(1000,nPOOL)
  } else {
    sample.size <- pop.sizes
  }
  if (missing(cover.matrix) || is.null(cover.matrix)){
    #print("Coverage matrix missing. Drawing genotypes from a binomial distribution.")
    geno <- NULL
    for (k in 1:nPOOL){
      G <- rbinom(n=(sample.size[k]*nSNP),size=ploidy,prob=as.numeric(pool.matrix[k,]))
      G[which(is.na(G))] <- 9
      geno <- rbind(geno,t(matrix(G,nrow=nSNP,sample.size[k])))
    }
  } else {
    #print("Coverage matrix not missing. Drawing genotypes from a beta-binomial distribution.")
    geno <- NULL
    for (k in 1:nPOOL){
      cover.1 <- cover.matrix[k,]
      n.reads <- array(NA,dim=c(1,nSNP))
      nna <- which(pool.matrix[k,]>0)
      na <- which(pool.matrix[k,]==0)
      epsilon <- 0.0001
      n.reads[nna] <- cover.1[nna]/pool.matrix[k,nna]
      n.reads[na] <- cover.1[na]/epsilon
      cover.2 <- n.reads - cover.1
      if (method=="per.pop"){
        p <- matrix(rbeta(sample.size[k]*nSNP,cover.1+1,cover.2+1),nrow=sample.size[k],ncol=nSNP)
        p.aux <- matrix(p,nrow=1)
        G <- t(matrix(rbinom(sample.size[k]*nSNP,size=ploidy,prob=p.aux),ncol=sample.size[k],nrow=nSNP))
        geno <- rbind(geno,G)
      } else if (method=="per.ind"){
        for (ind in 1:sample.size[k]){
          p <- rbeta(nSNP,cover.1+1,cover.2+1)
          G <- matrix(rbinom(nSNP,size=ploidy,prob=p),nrow=1,ncol=nSNP)
          geno <- rbind(geno,G)
        }
      }
    }
  }
  return(geno)
}

#' Simulate frequency matrix from genotype and coverage matrices
#'
#' \code{cover.to.pool} creates a matrix of frequency estimates, given a genotype matrix and a coverage matrix.
#' 
#' @param data a matrix with n rows and p columns where n is the number of individuals and p is the number of markers.  
#' @param cover.matrix a matrix with n rows and p columns where n is the number of pools and is the number of markers.
#' @param pop a list of integers or strings specifying which subpopulation the individuals belong to.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' 
#' @importFrom stats runif
#' 
#' @export
#'
cover.to.pool = function(data,cover.matrix,pop,ploidy=2){
  nPOP <- nrow(cover.matrix)
  nSNP <- ncol(cover.matrix)
  pool.matrix <- array(0,dim=c(nPOP,nSNP))
  pop.lab <- unique(pop)
  for (n in 1:nPOP){
    idx <- which(pop==pop.lab[n])
    nIND.pop.n <- length(idx)
    for (p in 1:nSNP){
      c.np <- cover.matrix[n,p]
      if (c.np > 0){
        draw.np <- floor(runif(c.np,min = 1,max = nIND.pop.n) + 1)
        drawn.geno <- data[idx[draw.np],p]
        pool.matrix[n,p] <- sum(drawn.geno,na.rm = TRUE)/(ploidy*c.np)
      } else {
        pool.matrix[n,p] <- NA
      }
    }
  }
  return(pool.matrix)
}

