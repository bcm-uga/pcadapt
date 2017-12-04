clumping_r = function(xptr, 
                      n, 
                      p, 
                      lookup_geno, 
                      lookup_byte, 
                      af, 
                      size = 100, 
                      thr = 0.2,
                      exclude = NULL) {
  
  maf   <- pmin(af, 1 - af)
  sumX  <- get_sumX(xptr, lookup_geno, lookup_byte)
  denoX <- get_denoX(xptr, lookup_geno, lookup_byte, af)
  
  ord <- order(maf, decreasing = TRUE)
  remain <- rep(TRUE, length(ord))
  remain[exclude] <- FALSE
  
  clumping(xptr, lookup_geno, lookup_byte,
           ord, remain, sumX, denoX, size, thr)
}

