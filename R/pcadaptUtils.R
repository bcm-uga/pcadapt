#' File Converter
#'
#' \code{read.pcadapt} converts \code{.vcf} and \code{.ped} files to an appropriate
#' type of file readable by \code{pcadapt}. You may find the converted file in the
#' current directory.
#'
#' @param input.filename a character string specifying the name of the file to be
#' converted if \code{local.env = FALSE}. If \code{local.env = TRUE}, \code{input.filename} refers
#' to the genotype matrix in the local environment.
#' @param type a character string specifying the type of data to be converted to the
#' \code{pcadapt} format. Supported formats are: \code{ped}, \code{vcf}, \code{lfmm}.
#' @param local.env a logical value indicating whether the input has to be read from the local
#' environment or from the working directory.
#' @param allele.sep a character string specifying the type of allele separator used in VCF files. Set to "/" by default, but can
#' be switched to "|" otherwise.
#' @param blocksize an integer specifying the number of markers to be processed in the mean time.
#'
#' @useDynLib pcadapt wrapper_converter
#' @importFrom utils tail
#' @importFrom vcfR read.vcfR extract.gt
#'
#' @export
#'
read.pcadapt <- function(input.filename,type,local.env=FALSE,allele.sep="/",blocksize=10000){
  if (local.env == FALSE){
    ## Check if input.filename is a character string ##
    if (class(input.filename) != "character"){
      stop("Argument input.filename has to be a character string.")  
    }
    ## Check if input.filename exists ##    
    if (!file.exists(input.filename)){
      stop(paste0("File ",input.filename," does not exist."))
    } 
    ## Check if argument type is missing ##
    if (missing(type)){
      stop("Argument type is missing.")
    }
    ## Check the file type ##
    if (class(type) != "character" || (!(type %in% c("vcf","ped","lfmm","pcadapt","vcfR")))){
      stop("Incorrect type.")
    }
    ## wrapper_converter ##
    if (type == "ped"){
      .C("wrapper_converter",as.character(input.filename),as.integer(0),PACKAGE = "pcadapt")
    } else if (type == "vcf"){
      ## The former vcf converter has been replaced ##
      #cat("If conversion fails, try type='vcfR' instead.\n")
      #.C("wrapper_converter",as.character(input.filename),as.integer(1),PACKAGE = "pcadapt")
      obj.vcf <- vcfR::read.vcfR(input.filename)
      geno <- vcfR::extract.gt(obj.vcf)
      vcf2pcadapt(geno,output.file = "tmp.pcadapt",allele.sep = allele.sep,blocksize = blocksize,console.count = blocksize)
    } else if (type == "lfmm"){
      .C("wrapper_converter",as.character(input.filename),as.integer(2),PACKAGE = "pcadapt")
    } else if (type == "vcfR"){
      obj.vcf <- vcfR::read.vcfR(input.filename)
      geno <- vcfR::extract.gt(obj.vcf)
      vcf2pcadapt(geno,output.file = "tmp.pcadapt",allele.sep = allele.sep,blocksize = blocksize,console.count = blocksize)
    }
    
    ##
    split.name <- unlist(unlist(strsplit(input.filename,"[.]")))
    if ((tail(split.name, n=1) %in% c("ped","vcf","lfmm","pcadapt")) && (length(split.name) > 1)){
      aux <- NULL
      for (k in (1:(length(split.name)-1))){
        aux <- paste0(aux,split.name[k],".")
      }
      aux <- paste0(aux,"pcadapt")
    } else if ((tail(split.name, n=1) != "lfmm") && type=="lfmm"){
      aux <- paste0(input.filename,".pcadapt")
    } else {
      aux <- input.filename
    }
    if (type == "vcfR"){
      aux <- "tmp.pcadapt"
    }
  } else if (local.env == TRUE){
    if (class(input.filename) %in% c("array","matrix","data.frame")){
      if (!(ncol(input.filename)>0) || !(nrow(input.filename)>0)){
        stop("Invalid input genotype matrix.")
      }
    } else {
      stop("Invalid input genotype matrix.")
    }
    if (type == "lfmm"){
      write.table(t(input.filename),"tmp.pcadapt",col.names = FALSE,row.names = FALSE)
    } else if (type == "pcadapt"){
      write.table(input.filename,"tmp.pcadapt",col.names = FALSE,row.names = FALSE)
    } else if (type == "vcf"){
      stop("Set argument local.env to FALSE.")
    } else if (type == "ped"){
      stop("Set argument local.env to FALSE.")
    }
    aux <- "tmp.pcadapt"
  }
  return(aux)
}

#' Population colorization
#'
#' \code{get.score.color} allows the user to display individuals of the same
#' pre-defined population with the same color when using the option
#' \code{"scores"} in \code{pcadapt}.
#'
#' @param pop a list of integers or strings specifying which population the
#' individuals belong to.
#'
#' @examples
#' ## see also ?pcadapt for examples
#'
#' @importFrom grDevices rainbow
#'
#' @keywords internal
#'
#' @export
#'
get.score.color = function(pop){
  pop.split <- list()
  list.ref <- unlist(pop)
  nIND <- length(list.ref)
  idx <- 1
  while (length(list.ref)>0){
    col <- list.ref[1]
    pop.split[[idx]] <- which(pop==col)
    idx <- idx+1
    list.ref <- list.ref[list.ref != col]
  }
  color.individuals <- array(dim=nIND)
  for (k in 1:length(pop.split)){
    color.individuals[unlist(pop.split[k])] <- k
  }
  return(color.individuals)
}

#' Retrieve population names
#'
#' \code{get.pop.names} retrieves the population names from the population file.
#'
#' @param pop a list of integers or strings specifying which population the
#' individuals belong to.
#'
#' @examples
#' ## see also ?pcadapt for examples
#'
#' @importFrom grDevices rainbow
#'
#' @keywords internal
#'
#' @export
#'
get.pop.names = function(pop){
  aux <- pop[1]
  res <- aux
  for (i in 1:(length(pop))){
    if (!(pop[i] %in% res)){
      aux <- c(aux,pop[i])
      res <- c(pop[i],res)
    }
  }
  return(aux)
}

#' Get the principal component the most associated with a genetic marker
#'
#' \code{get.pc} returns a data frame such that each row contains the index of
#' the genetic marker and the principal component the most correlated with it.
#'
#' @param x an object of class `pcadapt`. 
#' @param list a list of integers corresponding to the indices of the markers of interest.
#'
#' @examples
#' ## see also ?pcadapt for examples
#'
#' @export
#'
get.pc <- function(x,list){
  rem.na <- which(!is.na(x$loadings[list,1]))
  v <- array(0,dim=length(list))
  v[rem.na] <- sapply(list[rem.na],FUN=function(l){which(x$loadings[l,]^2==max(x$loadings[l,]^2,na.rm=TRUE))})
  df <- cbind(list,lapply(v,as.numeric))
  colnames(df) <- c("SNP","PC")
  return(df)
}

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
