#' vcfR-based converter
#'
#' \code{vcf2pcadapt} uses the package \code{vcfR} to extract the genotype information from a vcf file and exports it
#' under the format required by \code{pcadapt}.
#'
#' @param input a character string specifying the name of the file to be
#' converted.
#' @param output a character string indicating the name of the output file.
#' @param allele.sep a vector of characters indicating what delimiters are used to separate alleles.
#' 
#' @examples
#' ## see also ?pcadapt for examples
#'
#' @importFrom vcfR read.vcfR extract.gt
#' @importFrom utils write.table
#'
#' @export
#'
vcf2pcadapt <- function(input, output = "tmp.pcadapt", allele.sep = c("/", "|")){
  if (file.exists(output)){
    file.remove(output)
  }
  obj.vcf <- vcfR::read.vcfR(input, verbose = FALSE)
  geno <- vcfR::extract.gt(obj.vcf)
  pos <- vcfR::getPOS(obj.vcf)
  nIND <- ncol(geno)
  nSNP <- nrow(geno)
  res <- vcf_convert(string_geno = geno, output = output, allele_sep = allele.sep)
  skipped <- sum(res)
  if (skipped > 0){
    cat(skipped, "variant(s) have been discarded as they are not SNPs.\n")
  } else {
    cat("No variant got discarded.\n")
  }
  print_convert(input = input, output = output, M = nSNP, N = nIND, pool = 0)
  write.table(pos[res == 0], "positions.txt", 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
  return(res)
}
