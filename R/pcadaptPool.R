#' Principal Component Analysis based on the correlation matrix
#'
#' \code{corpca} is an auxiliary function that performs principal components analysis on a dataset. It returns an object \code{x}
#' which contains the loadings, the scores and the singular values of the \code{K} first principal components.
#' It handles missing values in a dataset and actually computes the eigen elements of the \code{n x n}
#' covariance matrix, where \code{n} is the number of individuals.
#'
#' @param data a data matrix or a data frame.
#' @param K an integer specifying the number of principal components that are retained.
#'
#' @importFrom stats cov
#'
#' @examples
#' x <- NULL
#'
#' @keywords internal
#'
#' @export
#' 
corpca = function(data,K){
  n <- dim(data)[1]
  p <- dim(data)[2]
  cat(paste0("Number of SNPs: ",p,"\n")) 
  cat(paste0("Number of populations: ",n,"\n")) 
  data_aux <- scale(data,scale=FALSE)*sqrt(p/(n-1))
  covmat <- cov(t(data_aux),use="pairwise.complete.obs")
  res <- NULL
  aux <- eigen(covmat,symmetric=TRUE)
  sdev <- aux$values[1:K]
  print(K)
  res$scores <- aux$vectors[,1:K]
  aux_ldgs <- t(aux$vectors)%*%data_aux
  res$loadings <- array(0,dim=c(p,K))
  res$loadings[,] <- t((1/(sqrt(sdev)))*aux_ldgs[1:K,])
  res$singular.values <- sqrt(abs(sdev))
  return(res)
}

#' pcadapt objects
#'
#' \code{create.pcadapt.pool} creates an object of class \code{pcadapt} for Pool-Seq data.
#'
#' @param data a data matrix or a data frame containing the allele frequencies per SNP and per population.
#' @param K an integer specifying the number of principal components to retain.
#' @param min.maf a value between \code{0} and \code{0.45} specifying the threshold
#' of minor allele frequencies above which p-values are computed.
#' @param cover.matrix a matrix specifying the average coverage per genetic marker and per population.
#'
#' @importFrom robust covRob
#' @importFrom MASS cov.rob
#' @importFrom stats median pchisq na.omit qchisq
#' 
#' @keywords internal
#' 
#' @export
#'
create.pcadapt.pool = function(data,K,min.maf,cover.matrix=NULL){
  nSNP <- ncol(data)
  nPOP <- nrow(data)
  # New procedure
  z.matrix <- as.matrix(array(NA,dim=c(nPOP,nSNP)))
  for (n in 1:nPOP){
    n_i <- cover.matrix[n,]
    nnan <- which(!is.na(n_i) & (n_i != 0))
    n_i[!nnan] <- 1
    f_i <- as.numeric(data[n,])
    f_i[!nnan] <- NA
    se <- as.vector(array(NA,dim=length(n_i)))
    se[nnan] <- sqrt(abs(f_i[nnan]*(1-f_i[nnan]))/n_i[nnan])
    z.matrix[n,nnan] <- f_i[nnan]/se[nnan] 
  }
#   for (k in 1:nSNP){
#     for (n in 1:nPOP){
#       n_i <- cover.matrix[n,k]
#       f_i <- data[n,k]
#       se <- sqrt(abs(f_i*(1-f_i))/n_i)
#       if ((!is.na(se)) && (se > 0)){
#         z.matrix[n,k] <- f_i/se 
#       } else {
#         z.matrix[n,k] <- f_i/0.01
#       }
#     }
#   }
  # End new procedure
  res <- corpca(data=data,K=K)
  freq <- apply(data,2,FUN=function(x){mean(x,na.rm=TRUE)})
  res$maf <- as.vector(pmin(freq,1-freq))
  res$loadings[res$maf<min.maf] <- NA 
  z.matrix[,res$maf<min.maf] <- NA
  res$stat <- array(NA,dim=nSNP)
  finite.list <- which(!is.na(apply(abs(z.matrix),2,sum)))
#   if (K>1){
#     res$stat[finite.list] <- as.vector(robust::covRob(res$loadings,na.action=na.omit,estim="pairwiseGK")$dist)
#     res$gif <- median(res$stat,na.rm=TRUE)/qchisq(0.5,df=K)
#   } else {
#     onedcov <- as.vector(MASS::cov.rob(res$loadings[finite.list,1]))
#     res$gif <- onedcov$cov[1]
#     res$stat <- (res$zscores[,1]-onedcov$center)^2
#   }
#   res$chi2.stat <- res$stat/res$gif
  res$stat <- as.vector(robust::covRob(t(z.matrix),na.action=na.omit,estim = "pairwiseGK")$dist)
  res$gif <- median(res$stat,na.rm=TRUE)/qchisq(0.5,df=K)
  res$chi2.stat <- res$stat/res$gif
  res$z.mat <- apply(z.matrix,MARGIN=2,FUN=function(x){sum(x^2,na.rm = TRUE)})
  # Compute p-values
  res$pvalues <- compute.pval(res$chi2.stat,K,method="mahalanobis")
  class(res) <- 'pcadapt'
  attr(res,"K") <- K
  attr(res,"method") <- "mahalanobis"
  attr(res,"data.type") <- "pool"
  attr(res,"min.maf") <- min.maf
  return(res)
}