################################################################################

getCode <- function(NA.VAL = 3L) {
  geno.raw <- as.logical(rawToBits(as.raw(0:255)))
  s <- c(TRUE, FALSE)
  geno1 <- geno.raw[s]
  geno2 <- geno.raw[!s]
  geno <- geno1 + geno2
  geno[geno1 & !geno2] <- NA.VAL
  dim(geno) <- c(4, 256)
  return(geno)
}

################################################################################

dim.xptr_bed <- function(x) c(attr(x, "n"), attr(x, "p"))

iram_and_reg <- function(input, K, min.maf, ploidy, LD.clumping) {
  
  # Get dimensions
  n <- dim(input)[1]  ## can't use nrow()
  p <- dim(input)[2]  ## can't use ncol()
  
  # Get allele frequencies
  # Uses a non-scaled lookup table
  af <- get_af(input) * (2 / ploidy)
  
  # Create a logical vector to locate SNPs with MAF >= min.maf
  maf <- pmin(af, 1 - af)
  ind.pass.af <- which(maf >= min.maf)
  
  if (!is.null(LD.clumping)) {
    size <- LD.clumping$size
    thr  <- LD.clumping$thr
    if (is.null(size) || is.null(thr))
      stop("Incorrect parameter 'LD.clumping'.")
    
    # Create a logical vector to locate SNPs that have been clumped
    # Take also number of NAs into account?
    ord <- order(maf[ind.pass.af], decreasing = TRUE)
    pass <- clumping(input, ind.pass.af, 
                     ord, rep(TRUE, length(ord)), 
                     size, thr)
    ind.pass <- ind.pass.af[pass]
  } else {
    ind.pass <- ind.pass.af
  }
  p2 <- length(ind.pass)
  
  # Get number of non-missing values per row and per column
  nb_nona <- nb_nona(input, ind.pass)
  
  ### SVD using RSpectra
  obj.svd <- RSpectra::svds(
    A = function(x, args) {
      # When filtering, the actual number of SNPs that we have is actually
      # sum(pass) and not p anymore
      pMatVec4(input, ind.pass, af, ploidy, x) / nb_nona$p * p2
    }, 
    Atrans = function(x, args) {
      # NB: nb_nona$n depends on 'pass' as well
      cpMatVec4(input, ind.pass, af, ploidy, x) / nb_nona$n * n
    },
    k = K, 
    dim = c(n, p2),
    opts = list(tol = 1e-4, maxitr = 100)
  )
  
  # Multiple Linear Regression is performed also on SNPs that have been clumped,
  # that is why we recompute the lookup table
  Z <- matrix(NA_real_, p, K)
  auxreg <- multLinReg(input, ind.pass.af, af, ploidy, obj.svd$u)

  #New subsetting because of possible zscores equal to NaN when 
  #there are monomorphic sites
  theNaN<-rowSums(is.na(auxreg))!=0
  Z[ind.pass.af[!theNaN], ] <- auxreg[!theNaN]
  ind.pass.af <- ind.pass.af[!theNaN]
  obj.svd$zscores <- Z
  
  V <- matrix(NA_real_, p, K)
  V[ind.pass, ] <- obj.svd$v
  obj.svd$v <- V
  
  obj.svd$snps.included <- ind.pass
  # obj.svd$d <- obj.svd$d^2 / length(ind.pass)
  obj.svd$af <- af
  # obj.svd$nona1 <- nb_nona$p
  # obj.svd$nona2 <- nb_nona$n
  
  obj.svd$pass <- ind.pass.af
  
  obj.svd
}

################################################################################
