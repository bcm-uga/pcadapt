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
    
    ### File mapping
    xptr <- bedXPtr(path_to_bed, n, p)
    
  } else if (class(input) == "matrix") {
    # an input matrix has nIND rows and nSNP columns
    xptr <- input
    n <- nrow(xptr)
    p <- ncol(xptr)
  }
  
  lookup_byte <- getCode()
  lookup_geno <- rbind(rep(0, p), 1, 2, 3)
  
  ### Get allele frequencies
  tmp <- get_af(xptr, lookup_geno, lookup_byte)
  
  ### Get number of non-missing values per row and per column
  nb_nona <- nb_nona(xptr, lookup_geno, lookup_byte)
  
  pass.af <- (pmin(tmp, 1 - tmp) >= min.maf)
  if (LD.clumping) {
    pass.LD <- clumping_r(xptr = xptr, 
                          n = n, 
                          p = p, 
                          lookup_geno = lookup_geno, 
                          lookup_byte = lookup_byte, 
                          af = tmp, 
                          size = size, 
                          thr = thr)
  } else {
    pass.LD <- rep(TRUE, length(pass.af))
  }
  pass <- (pass.LD & pass.af)
  
  ### Lookup table
  lookup_scale <- rbind(outer(0:2, tmp, function(g, p) {
    if (p > 0 || p < 1) {
      return((g - 2 * p) / sqrt(2 * p * (1 - p)))
    } else {
      return(0)
    }
  }), 0)
  
  lookup_scale[, !pass] <- 0
  
  ### SVD using RSpectra
  obj.svd <- RSpectra::svds(
    A = function(x, args) {
      pMatVec4(xptr, x, lookup_scale, lookup_byte) / nb_nona[[1]] * sum(pass)
    }, 
    Atrans = function(x, args) {
      cpMatVec4(xptr, x, lookup_scale, lookup_byte) / nb_nona[[2]] * n
    },
    k = K, 
    dim = c(n, p),
    opts = list(tol = 1e-4)
  )
  
  obj.svd$zscores <- multLinReg(xptr, 
                                lookup_scale, 
                                lookup_byte, 
                                obj.svd$u, 
                                obj.svd$d, 
                                obj.svd$v)
  obj.svd$pass <- pass
  obj.svd$d <- obj.svd$d^2 / p
  obj.svd$maf <- pmin(tmp, 1 - tmp)
  return(obj.svd)
}