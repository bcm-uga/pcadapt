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

iram_and_reg <- function(input, 
                 K = 2, 
                 min.maf = 0.05, 
                 LD.clumping = FALSE, 
                 size = 100, 
                 thr = 0.2) {  # TODO: add exclude parameter
  
  if (class(input) == "character") {
    path_to_bed <- normalizePath(input)
    p <- nrow(data.table::fread(sub("\\.bed$", ".bim", path_to_bed)))
    n <- nrow(data.table::fread(sub("\\.bed$", ".fam", path_to_bed)))
    
    # File mapping
    xptr <- bedXPtr(path_to_bed, n, p)
    
  } else if (class(input) == "matrix") {
    # an input matrix has nIND rows and nSNP columns
    xptr <- input
    n <- nrow(xptr)
    p <- ncol(xptr)
  }
  
  # Get allele frequencies
  # Uses a non-scaled lookup table
  af <- get_af(xptr)
  
  # Create a logical vector to locate SNPs with mAF >= min.maf
  maf <- pmin(af, 1 - af)
  ind.pass.af <- which(maf >= min.maf)
  
  # Create a logical vector to locate SNPs that have been clumped
  if (LD.clumping) {
    # Take also number of NAs into account?
    ord <- order(maf[ind.pass.af], decreasing = TRUE)
    pass <- clumping(xptr, ind.pass.af, 
                     ord, rep(TRUE, length(ord)), 
                     size, thr)
    ind.pass <- ind.pass.af[pass]
  } else {
    ind.pass <- ind.pass.af
  }
  p2 <- length(ind.pass)
  
  # Get number of non-missing values per row and per column
  nb_nona <- nb_nona(xptr, ind.pass)
  
  ### SVD using RSpectra
  obj.svd <- RSpectra::svds(
    A = function(x, args) {
      # When filtering, the actual number of SNPs that we have is actually
      # sum(pass) and not p anymore
      pMatVec4(xptr, ind.pass, af, x) / nb_nona$p * p2
    }, 
    Atrans = function(x, args) {
      # NB: nb_nona$n depends on 'pass' as well
      cpMatVec4(xptr, ind.pass, af, x) / nb_nona$n * n
    },
    k = K, 
    dim = c(n, p2),
    opts = list(tol = 1e-4, maxitr = 100)
  )
  
  # Multiple Linear Regression is performed also on SNPs that have been clumped,
  # that is why we recompute the lookup table
  obj.svd$zscores <- multLinReg(xptr, ind.pass.af, af, obj.svd$u)
  
  V <- matrix(NA_real_, p, K)
  V[ind.pass, ] <- obj.svd$v
  obj.svd$v <- V
  
  obj.svd$snps.included <- ind.pass
  # obj.svd$d <- obj.svd$d^2 / length(ind.pass)
  obj.svd$af <- af
  # obj.svd$nona1 <- nb_nona$p
  # obj.svd$nona2 <- nb_nona$n
  
  obj.svd
}

################################################################################

