################################################################################

clumping_r = function(xptr, lookup_geno, lookup_byte, ind_col, maf, size, thr) {
  
  ord <- order(maf, decreasing = TRUE)
  remain <- rep(TRUE, length(ord))
  clumping(xptr, lookup_geno, lookup_byte, ind_col, ord, remain, size, thr)
}

################################################################################