#' pcadapt objects
#'
#' \code{create.pcadapt} loads the numerical quantities needed to compute the test
#' statistics, and stores them in an object of class \code{pcadapt}.
#'
#' @param output.filename a character string indicating which outputs from PCAdapt
#' fast should be processed.
#' @param K an integer specifying the number of principal components to retain.
#' @param method a character string specifying the method to be used to compute
#' the p-values. Four statistics are currently available, \code{"mahalanobis"},
#' \code{"communality"}, and \code{"componentwise"}.
#' @param data.type a character string specifying the type of data being read, either
#' a \code{genotype} matrix (\code{data.type="genotype"}), or a matrix of allele
#' frequencies (\code{data.type="pool"}).
#' @param min.maf a value between \code{0} and \code{0.45} specifying the threshold
#' of minor allele frequencies above which p-values are computed.
#'
#' @importFrom robust covRob
#' @importFrom MASS cov.rob
#' @importFrom stats median pchisq na.omit qchisq
#' @importFrom utils read.table
#' 
#' @export
#'
create.pcadapt = function(output.filename,K,method,data.type,min.maf){

  # Load outputs from PCAdapt fast
  res <- NULL
  res$maf <- read.table(paste0(output.filename,".maf"))[,1]
  nSNP <- length(res$maf)
  res$loadings <- as.matrix(read.table(paste0(output.filename,".loadings")))*sqrt(nSNP)
  res$singular.values <- as.numeric(read.table(paste0(output.filename,".sigma")))
  #res$loadings <- as.matrix(read.table(paste0(output.filename,".loadings")))*nSNP*res$singular.values
  res$scores <- t(read.table(paste0(output.filename,".scores")))
  nIND <- nrow(res$scores)
  res$loadings[res$maf<min.maf,] <- NA
  res$zscores <- as.matrix(read.table(paste0(output.filename,".zscores")))
  Kbis <- ncol(res$zscores)
  finite.list <- which(!is.na(apply(abs(res$zscores),1,sum)))
  res$stat <- array(NA,dim=nSNP)
  
  if (method == "mahalanobis"){
    # Use covRob (robust) for K>1 and cov.rob (MASS) for K=1
    if (Kbis>1){
      if (K != Kbis){
        warning("Number of principal components inconsistent.")
      }
      res$stat[finite.list] <- as.vector(robust::covRob(res$zscores,na.action=na.omit,estim="pairwiseGK")$dist)
    } else if (Kbis==1){
      onedcov <- as.vector(MASS::cov.rob(res$zscores[finite.list,1]))
      res$stat <- (res$zscores[,1]-onedcov$center)^2/onedcov$cov[1]
    }
    df <- K
    res$gif <- median(res$stat,na.rm=TRUE)/qchisq(0.5,df=df)
    res$chi2.stat <- res$stat/res$gif
  } else if (method == "communality"){
    res$stat <- sapply(1:nSNP,FUN=function(xx){sum(res$loadings[xx,1:K]^2*res$singular.values[1:K]^2/(nSNP))})
    c <- sum(res$singular.values[1:K]^2)/K
    df <- K
    res$gif <- median(res$stat*nSNP/c,na.rm=TRUE)/qchisq(0.5,df=K)
    res$chi2.stat <- res$stat*nSNP/(c*res$gif)
  } else if (method == "componentwise"){
    res$stat <- apply(res$zscores,MARGIN=2,FUN=function(xx){xx^2})
    df <- 1
    res$gif <- sapply(1:K,FUN=function(xx){median(res$zscores[,xx]^2,na.rm=TRUE)/qchisq(0.5,df=df)})
    res$chi2.stat <- res$stat/res$gif
  }
  
  # Compute p-values
  res$pvalues <- compute.pval(res$chi2.stat,K,method)
  
  class(res) <- 'pcadapt'
  attr(res,"K") <- K
  attr(res,"method") <- method
  attr(res,"data.type") <- data.type
  attr(res,"min.maf") <- min.maf

  return(res)
}

#' Principal Components p-values
#'
#' \code{compute.pval} computes p-values for each genetic marker, depending on the
#' test statistics distribution.
#'
#' @param st a chi-square statistic.
#' @param K an integer specifying the number of principal components to retain.
#' @param method a character string specifying the method to be used to compute
#' the p-values.
#'
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom stats pchisq
#'
#' @export
compute.pval = function(st,K,method){
  if (method == "componentwise"){
    p <- NULL
    for (k in 1:K) {
      p <- cbind(p,pchisq(st[,k],df=1,lower.tail=FALSE))
    }
  } else {
    p <- as.numeric(pchisq(st,df=K,lower.tail=FALSE))
  }
  return(p)
}


