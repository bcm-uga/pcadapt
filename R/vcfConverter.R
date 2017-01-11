
#' Convert genotype information obtained with vcfR
#'
#' \code{convert.line} converts outputs of \code{extract.gt} to the format \code{pcadapt}.
#'
#' @param hap.block a genotype matrix obtained with `vcfR`. 
#' @param allele.sep a character indicating which type of delimiter is used to separate the alleles.
#'
#' @examples
#' ## see also ?pcadapt for examples
#'
#' @keywords internal
#'
#' @export
#'
convert.line <- function(hap.block,allele.sep="/"){
  geno.block <- array(9,dim=dim(hap.block))
  if (allele.sep == "/"){
    geno.block[which(hap.block=="0/0",arr.ind = TRUE)] <- 0
    geno.block[which(hap.block=="0/1",arr.ind = TRUE)] <- 1
    geno.block[which(hap.block=="1/0",arr.ind = TRUE)] <- 1
    geno.block[which(hap.block=="1/1",arr.ind = TRUE)] <- 2
    mask <- !as.logical(apply(hap.block,MARGIN=1,FUN=function(x){sum(!(x %in% c("0/0","0/1","1/0","1/1","./.",NA)))}))
  } else if (allele.sep == "|"){
    geno.block[which(hap.block=="0|0",arr.ind = TRUE)] <- 0
    geno.block[which(hap.block=="0|1",arr.ind = TRUE)] <- 1
    geno.block[which(hap.block=="1|0",arr.ind = TRUE)] <- 1
    geno.block[which(hap.block=="1|1",arr.ind = TRUE)] <- 2
    mask <- !as.logical(apply(hap.block,MARGIN=1,FUN=function(x){sum(!(x %in% c("0|0","0|1","1|0","1|1",".|.",NA)))}))      
  }
  filtered.geno <- geno.block[mask,]
  return(filtered.geno)
}

#' vcfR-based converter
#'
#' \code{vcf2pcadapt} uses the package \code{vcfR} to extract the genotype information from a vcf file and exports it
#' under the format required by \code{pcadapt}.
#'
#' @param geno a genotype matrix obtained with `vcfR`.
#' @param output.file a character string indicating the name of the output file.
#' @param allele.sep a character indicating which type of delimiter is used to separate the alleles.
#' @param blocksize an integer indicating the number of SNPs to be processed in the same chunk.
#' @param console.count an integer indicating the number of SNPs by which the progress bar should increment.
#'
#' @examples
#' ## see also ?pcadapt for examples
#'
#' @importFrom utils flush.console
#'
#' @keywords internal
#'
#' @export
#'
vcf2pcadapt <- function(geno,output.file="tmp.pcadapt",allele.sep="/",blocksize=10000,console.count=10000){
  if (file.exists(output.file)){
    file.remove(output.file)
  }
  nIND <- ncol(geno)
  nSNP <- nrow(geno)
  init.block <- 1
  if (nSNP > blocksize){
    end.block <- blocksize
  } else {
    end.block <- nSNP
  }
  skipped <- 0
  while (init.block < nSNP && end.block <= nSNP){
    snp.lines <- geno[init.block:end.block,]
    hap2geno <- t(convert.line(snp.lines,allele.sep = allele.sep))
    if (nSNP > blocksize){
      skipped <- skipped + blocksize - ncol(hap2geno)
    } else {
      skipped <- nSNP - ncol(hap2geno)
    }
    write(hap2geno,file = output.file,append = TRUE,ncolumns = nIND)  
    if (end.block%%console.count == 0){
      cat('\r',end.block,"SNPs have been processed.")
      flush.console() 
    }
    if (end.block+blocksize <= nSNP){
      init.block <- end.block + 1
      end.block <- end.block + blocksize
    } else {
      init.block <- end.block + 1
      end.block <- nSNP
    }
  }
  cat('\r',end.block,"variants have been processed.\n",skipped,"variants have been discarded as they are not SNPs.")
}
