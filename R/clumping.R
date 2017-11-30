clumping_r = function(xptr, 
                      n, 
                      p, 
                      lookup_geno, 
                      lookup_byte, 
                      af, 
                      size = 100, 
                      thr = 0.2) {
  
  S <- pmin(af, 1 - af)
  sumX <- get_sumX(xptr, lookup_geno, lookup_byte)
  denoX <- get_denoX(xptr, lookup_geno, lookup_byte, af)
  
  ord.chr <- order(S, decreasing = TRUE)
  remain <- rep(TRUE, length(ord.chr))
  ind.keep <- clumping(xptr,
                       lookup_geno,
                       lookup_byte,
                       ord.chr,
                       remain,
                       sumX,
                       denoX,
                       size, 
                       thr)
  
  return(ind.keep)
  
}

