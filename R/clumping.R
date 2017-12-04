################################################################################

clumping_r = function(xptr, n, p, lookup_geno, lookup_byte, ind_col, 
                      af, size = 100, thr = 0.2) {
  
  maf   <- pmin(af, 1 - af)
  ord <- order(maf, decreasing = TRUE)
  remain <- rep(TRUE, length(ord))
  
  tmpX <- get_sumX_denoX(xptr, lookup_geno, lookup_byte, ind_col, af)
  clumping(xptr, lookup_geno, lookup_byte, ind_col,
           ord, remain, tmpX$sumX, tmpX$denoX, size, thr)
}

################################################################################