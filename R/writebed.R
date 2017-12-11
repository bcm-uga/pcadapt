################################################################################

getInverseCode <- function() {
    
  geno <- getCode()
  r <- raw(256); dim(r) <- rep(4, 4)
  for (i in 1:256) {
    ind <- geno[, i] + 1
    r[ind[1], ind[2], ind[3], ind[4]] <- as.raw(i - 1)
  }
  r
}

CODE_0123 <- mmapcharr:::CODE_012; CODE_0123[is.na(CODE_0123)] <- 3L

################################################################################

write.table2 <- function(x, file, sep = "\t") {
  utils::write.table(x = x, file = file, sep = sep, quote = FALSE,  
                     row.names = FALSE, col.names = FALSE)
}

################################################################################

write_fake_bim_fam <- function(n, m, bedfile) {
  
  # fam
  fam <- data.frame(0L, paste0("ind_", 1:n), 0L, 0L, 0L, -9L,
                    stringsAsFactors = FALSE)
  famfile <- sub("\\.bed$", ".fam", bedfile)
  write.table2(fam, famfile)
  
  # map
  map <- data.frame(1L, paste0("snp_", 1:m), 0L, 0L,
                    ifelse(cond <- (stats::runif(m) > 0.5), "A", "T"),
                    ifelse(!cond, "A", "T"),
                    stringsAsFactors = FALSE)
  bimfile <- sub("\\.bed$", ".bim", bedfile)
  write.table2(map, bimfile)
}

################################################################################

#' Write PLINK files
#'
#' Function to write bed/bim/fam files from a pcadapt or an lfmm file.
#' Files shouldn't already exists.
#'
#' @param x A [mmapchar][mmapchar-class] object associated 
#'   with a pcadapt or lfmm file.
#' @param is.pcadapt
#'
#' @return The input `bedfile` path.
#' 
#' @export
#' 
writeBed <- function(x, is.pcadapt) {
  
  # Get path to new bed file
  file <- x$backingfile
  bedfile <- paste0(file, ".bed")
  if (file.exists(bedfile)) stop("The bed file already exists!", call. = FALSE)
  
  # Write files
  ## Write bed file
  writebed(bedfile, x$copy(code = CODE_0123), getInverseCode(), is.pcadapt)
  ## Write other files
  write_fake_bim_fam(n = `if`(is.pcadapt, ncol(x), nrow(x)), 
                     m = `if`(is.pcadapt, nrow(x), ncol(x)), 
                     bedfile = bedfile)
  
  bedfile
}

################################################################################