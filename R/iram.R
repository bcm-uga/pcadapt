#' @export
#'
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

#' @export
#'
iram = function(input) {
  
  lookup_byte <- getCode()
  
  if (class(input) == "character") {
    path_to_bed <- normalizePath(input)
    p <- nrow(data.table::fread(sub("\\.bed$", ".bim", path_to_bed)))
    n <- nrow(data.table::fread(sub("\\.bed$", ".fam", path_to_bed)))
    
    ### File mapping
    xptr <- bedXPtr(path_to_bed, n, p)
    
    ### Get number of missing values per row and per column
    nb_nona <- nb_nona(xptr, rbind(rep(0, p), 1, 2, 3), lookup_byte)
    
    ### Get allele frequencies
    tmp <- af(xptr, rbind(rep(0, p), 1, 2, 3), lookup_byte)
    
    ### Lookup table
    lookup_scale <- rbind(outer(0:2, tmp, function(g, p) {
      (g - 2 * p) / sqrt(2 * p * (1 - p))
    }), 0)
    
    obj.svd <- RSpectra::svds(
      A = function(x, args) {
        cat(".")
        pMatVec4(xptr, x, lookup_scale, lookup_byte) / nb_nona[[1]] * n
      }, 
      k = 5, 
      Atrans = function(x, args) {
        cpMatVec4(xptr, x, lookup_scale, lookup_byte) / nb_nona[[2]] * p
      },
      dim = c(n, p),
      opts = list(tol = 1e-4)
    )
    
    return(obj.svd)
    
  } else if (class(obj) == "matrix") {
    cat("Not implemented yet.")
  }
}