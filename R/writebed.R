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

# write.table2 <- function(x, file, sep = "\t") {
#   utils::write.table(x = x, file = file, sep = sep, quote = FALSE,  
#                      row.names = FALSE, col.names = FALSE)
# }

################################################################################

# write_fake_bim_fam <- function(n, m, bedfile) {
#   
#   # fam
#   fam <- data.frame(0L, paste0("ind_", 1:n), 0L, 0L, 0L, -9L,
#                     stringsAsFactors = FALSE)
#   famfile <- sub("\\.bed$", ".fam", bedfile)
#   write.table2(fam, famfile)
#   
#   # map
#   map <- data.frame(1L, paste0("snp_", 1:m), 0L, 0L,
#                     ifelse(cond <- (stats::runif(m) > 0.5), "A", "T"),
#                     ifelse(!cond, "A", "T"),
#                     stringsAsFactors = FALSE)
#   bimfile <- sub("\\.bed$", ".bim", bedfile)
#   write.table2(map, bimfile)
# }

################################################################################

#' Write PLINK files
#'
#' Function to write bed/bim/fam files from a pcadapt or an lfmm file.
#' Files shouldn't already exists.
#'
#' @param file A [mmapchar][mmapchar-class] object associated 
#'   with a pcadapt or lfmm file.
#' @param is.pcadapt a boolean value.
#'
#' @return The input `bedfile` path.
#' 
#' @export
#' 
writeBed <- function(file, is.pcadapt) {
  
  # Map file
  mmap <- mmapcharr::mmapchar(file, code = CODE_0123)
  ## Get dimensions
  n <- `if`(is.pcadapt, ncol(mmap), nrow(mmap))
  p <- `if`(is.pcadapt, nrow(mmap), ncol(mmap))
  
  # Get path to new bed file
  bedfile <- paste0(file, ".bed")
  if (file.exists(bedfile)) {
    message("The bed file already exists. Returning..")
  } else {
    writebed(bedfile, mmap, getInverseCode(), is.pcadapt)
  }
  
  structure(normalizePath(bedfile), n = n, p = p, class = "pcadapt_bed")
}

################################################################################