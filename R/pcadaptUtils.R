#' File Converter
#'
#' \code{read.pcadapt} converts \code{.vcf} and \code{.ped} files to an appropriate
#' type of file readable by \code{pcadapt}. You may find the converted file in the
#' current directory.
#'
#' @param input a character string specifying the name of the file to be
#' converted if \code{local.env = FALSE}.
#' @param type a character string specifying the type of data to be converted to the
#' \code{pcadapt} format. Supported formats are: \code{ped}, \code{vcf}, \code{lfmm}.
#' @param local.env a logical value indicating whether the input has to be read from the local
#' environment or from the working directory.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' @param pop.sizes a vector specifying the number of individuals for each pool.
#' @param allele.sep a character string specifying the type of allele separator used in VCF files. Set to "/" by default, but can
#' be switched to "|" otherwise.
#' @param blocksize an integer specifying the number of markers to be processed in the mean time.
#'
#' @importFrom utils tail
#' @importFrom vcfR read.vcfR extract.gt
#'
#' @export
#'
read.pcadapt <- function(input, 
                         type, 
                         local.env = FALSE, 
                         ploidy, 
                         pop.sizes = NULL,
                         allele.sep = "/",
                         blocksize = 1000){
  
  ## In version 3.1.0, argument output.filename has been removed ##
  if (!missing(local.env)){
    warning("Argument local.env is deprecated. Please refer to the latest vignette for further information.")
  }
  
  if (class(input) == "character"){
    ## Check if input exists ##    
    if (!file.exists(input)){
      stop(paste0("File ", input, " does not exist."))
    } 
    ## Check if argument type is missing ##
    if (missing(type)){
      stop("Argument type is missing.")
    }
    ## Check if file type is supported ##
    if (class(type) != "character" || (!(type %in% c("vcf", "ped", "lfmm", "pcadapt")))){
      stop("Incorrect type.")
    }
    
    ## x.type to x.pcadapt ##
    split.name <- unlist(unlist(strsplit(input,"[.]")))
    if ((tail(split.name, n = 1) %in% c("ped", "vcf", "lfmm", "pcadapt")) && (length(split.name) > 1)){
      aux <- NULL
      for (k in (1:(length(split.name) - 1))){
        aux <- paste0(aux, split.name[k], ".")
      }
      aux <- paste0(aux, "pcadapt")
    } else {
      aux <- paste0(input, ".pcadapt")
    }
    
    ## File converter ##
    if (type == "ped"){
      ped2pcadapt(input)
    } else if (type == "vcf"){
      obj.vcf <- vcfR::read.vcfR(input)
      geno <- vcfR::extract.gt(obj.vcf)
      vcf2pcadapt(geno, output.file = aux, allele.sep = allele.sep, blocksize = blocksize, console.count = blocksize)
    } else if (type == "lfmm"){
      stop("Not done yet.")
    }
  } else if ((class(input) %in% c("matrix", "data.frame", "array"))){
    if (!(ncol(input) > 0) || !(nrow(input) > 0)){
      stop("Invalid input genotype matrix.")
    }
    if (type == "lfmm"){
      tmp <- t(as.matrix(input))
      tmp[which(is.na(tmp))] <- 9
    } else if (type == "pcadapt"){
      tmp <- as.matrix(input)
      tmp[which(is.na(tmp))] <- 9
    } else if (type == "vcf"){
      stop("Incorrect type.")
    } else if (type == "ped"){
      stop("Incorrect type.")
    } else if (type == "pool"){
      if (missing(ploidy)){
        warning("Argument ploidy is missing, proceeding with 'ploidy = 2'...")
        p <- 2
      } else {
        p <- ploidy
      }
      if (missing(pop.sizes)){
        warning("Argument pop.sizes is missing, proceeding with 100 individuals per pool.")
        nPOOL <- nrow(input)
        s <- as.vector(rep(100, nPOOL))
      } else {
        s <- as.vector(pop.sizes)
      }
      tmp <- sample_geno_cpp(freq = as.matrix(input), ploidy = p, sample_size = s)
      tmp[which(is.na(tmp))] <- 9
    } 
    aux <- tmp
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
