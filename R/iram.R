getCode = function(NA.VAL = 3L) {
  geno.raw <- as.logical(rawToBits(as.raw(0:255)))
  s <- c(TRUE, FALSE)
  geno1 <- geno.raw[s]
  geno2 <- geno.raw[!s]
  geno <- geno1 + geno2
  geno[geno1 & !geno2] <- NA.VAL
  dim(geno) <- c(4, 256)
  return(geno)
}

iram = function(input, 
                K = 2, 
                min.maf = 0.05, 
                LD.clumping = FALSE, 
                size = 100, 
                thr = 0.2) {
  
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
  
  lookup_byte <- getCode()
  lookup_geno <- rbind(rep(0, p), 1, 2, 3)
  
  # Get allele frequencies
  # Uses a non-scaled lookup table
  af <- get_af(xptr, lookup_geno, lookup_byte)
  
  # Create a logical vector to locate SNPs with mAF >= min.maf
  pass.af <- (pmin(af, 1 - af) >= min.maf)
  
  # Create a logical vector to locate SNPs that have been clumped
  if (LD.clumping) {
    pass <- clumping_r(xptr = xptr, 
                       n = n, 
                       p = p, 
                       lookup_geno = lookup_geno, 
                       lookup_byte = lookup_byte, 
                       ind_col = which(pass.af),
                       af = af[pass.af], 
                       size = size, 
                       thr = thr,
                       exclude = !pass.af)
  } else {
    pass <- pass.af
  }
  
  ind.pass <- which(pass)
  
  # Get number of non-missing values per row and per column
  # Uses a non-scaled lookup table
  nb_nona <- nb_nona(xptr, lookup_geno, lookup_byte, ind.pass)
  
  # Scaled lookup table 
  lookup_geno <- rbind(outer(0:2, af, function(g, p) {
    if (p > 0 || p < 1) {
      return((g - 2 * p) / sqrt(2 * p * (1 - p)))
    } else {
      return(0)
    }
  }), 0)
  
  # SNP filtering: we assign the value 0 to all SNPs that have a mAF lower
  # than min.maf or that have been clumped
  lookup_geno[, ind.pass] <- 0
  
  ### SVD using RSpectra
  obj.svd <- RSpectra::svds(
    A = function(x, args) {
      # When filtering, the actual number of SNPs that we have is actually
      # sum(pass) and not p anymore
      pMatVec4(xptr, x, lookup_geno, lookup_byte, ind.pass) / 
        nb_nona$p * length(ind.pass)
    }, 
    Atrans = function(x, args) {
      # NB: nb_nona$n depends on 'pass' as well
      cpMatVec4(xptr, x, lookup_geno, lookup_byte, ind.pass) / 
        nb_nona$n * n
    },
    k = K, 
    nv = 0,  
    dim = c(n, p),
    opts = list(tol = 1e-4, maxitr = 100)
  )
  
  # We calculate the loadings even for the SNPs that have been clumped.
  # The loadings are passed by reference and are computed in the multLinReg
  # function (not in svds)
  obj.svd$v <- matrix(0, nrow = p, ncol = K)
  
  # Lookup table
  lookup_geno <- rbind(outer(0:2, af, function(g, p) {
    if (p > 0 || p < 1) {
      return((g - 2 * p) / sqrt(2 * p * (1 - p)))
    } else {
      return(0)
    }
  }), 0)
  
  # Multiple Linear Regression is performed also on SNPs that have been clumped,
  # that is why we recompute the lookup table
  lookup_geno[, !pass.af] <- 0
  
  obj.svd$zscores <- multLinReg(xptr,
                                lookup_geno,
                                lookup_byte,
                                obj.svd$u,
                                obj.svd$d,
                                obj.svd$v)
  
  obj.svd$pass <- pass.af
  obj.svd$d <- obj.svd$d^2 / sum(pass)
  obj.svd$maf <- pmin(af, 1 - af)
  obj.svd$nona1 <- nb_nona$p
  obj.svd$nona2 <- nb_nona$n
  
  obj.svd
}


iram2 = function(X, k) {
  p <- apply(X, MARGIN = 1, FUN = function(h) {mean(h, na.rm = TRUE) / 2})
  A <- function(x, args) {
    return(prodtGx(X, x, p)) # Input vector of length p
  }
  Atrans <- function(x, args) {
    return(prodGx(G, x, p)) # Input vector of length n
  }
  res <- RSpectra::svds(A, k, nu = k, nv = 0, Atrans = Atrans,
                        opts = list(tol = 1e-4, maxitr = 100),
                        dim = c(ncol(X), nrow(X)))
  res$maf <- pmin(p, 1 - p)
  return(res)
}