################################################################################

#' File Converter
#'
#' \code{read.pcadapt} converts genotype matrices or files to an appropriate
#' format readable by \code{pcadapt}. You may find the converted file in the
#' current directory.
#'
#' @param input a genotype matrix or a character string specifying the name of the file to be
#' converted.
#' @param type.input a character string specifying the type of data to be converted to the
#' \code{pcadapt} format. Supported formats are: \code{ped}, \code{vcf}, \code{lfmm}.
#' Deprecated argument.
#' @param local.env deprecated argument.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' @param pop.sizes a vector specifying the number of individuals to be sampled for each pool.
#' @param allele.sep a vector of characters indicating what delimiters are used in VCF files. By default, only
#' "|" and "/" are recognized.
#' @param blocksize deprecated argument.
#'
#' @importFrom utils tail
#'
#' @export
#'
read.pcadapt <- function(input, 
                         type.input = c("pcadapt", "lfmm", "vcf", 
                                        "ped", "pool", "example"), 
                         type.output = c("bed", "matrix"),
                         pop.sizes = NULL,
                         allele.sep = c("/", "|"), 
                         ploidy = 2, 
                         local.env, 
                         blocksize){
  
  ## In version 3.1.0, argument local.env has been removed ##
  if (!missing(local.env)) {
    warning("Argument local.env is deprecated. Please refer to the latest vignette for further information.")
  }
  
  ## In version 3.1.0, argument local.env has been removed ##
  if (!missing(blocksize)) {
    warning("Argument blocksize is deprecated. Please refer to the latest vignette for further information.")
  }
  
  if (type.input == "pool") {
    return(read.pcadapt.pool(input, type.output, pop.sizes, ploidy))
  } 
  
  if (class(input) == "character") {
    
    ## Check if input exists ##    
    if ((type.input != "example") && !file.exists(input)) 
      stop(paste("File", input, "does not exist."))
    
    ## File converter ##
    is.pcadapt <- TRUE
    if (type.input == "ped") {
      check_file_size(input)
      tmp <- tempfile(fileext = ".pcadapt")
      ped2pcadapt(input = input, output = tmp)
      input <- tmp
    } else if (type.input == "vcf") {
      check_file_size(input)
      tmp <- tempfile(fileext = ".pcadapt")
      vcf2pcadapt(input = input, output = tmp, allele.sep = allele.sep)
      input <- tmp
    } else if (type.input == "lfmm") {
      is.pcadapt <- FALSE
    } else if (type.input == "example") {
      input <- system.file("extdata", "geno3pops.pcadapt", package = "pcadapt")
    }
    
    return(pcadapt2other(input, is.pcadapt, type.output))
    
  } else if (class(input) %in% c("matrix", "data.frame", "array")) {
    
    if (!(ncol(input) > 0) || !(nrow(input) > 0)){
      stop("Invalid input genotype matrix.")
    }
    if (type.input == "lfmm") {
      tmp <- t(as.matrix(input))
      storage.mode(tmp) <- "integer"
      # tmp[which(is.na(tmp))] <- 9
    } else if (type.input == "pcadapt") {
      tmp <- as.matrix(input)
      storage.mode(tmp) <- "integer"
      # tmp[which(is.na(tmp))] <- 9
    } else {
      stop("Incorrect type.input for matrices.")
    }
    aux <- tmp
  } else {
    stop("Input should be a file path or a matrix-like object.")
  }
  
  aux
}

################################################################################

check_file_size <- function(file) {
  if (file.size(file) > 2e9) {
    stop(paste0("It seems that the input file is quite large.\n",
                "For large 'vcf' or 'ped' files, ", 
                "please use PLINK to convert them in 'bed'."))
  }
}

################################################################################

pcadapt2other <- function(file, is.pcadapt = TRUE, type.return = "bed") {
  
  if (type.return == "matrix") {
    mmap <- mmapcharr::mmapchar(file, code = mmapcharr:::CODE_012)
    `if`(is.pcadapt, t(mmap[]), mmap[])
  } else { # bed
    writeBed(file, is.pcadapt)
  }
}

################################################################################

read.pcadapt.pool <- function(input, pop.sizes, ploidy) {
  
  if (class(input) == "character") {
    
    fs <- get_size_cpp(input)
    nPOOL <- fs[1]
    if (length(pop.sizes) == nPOOL) {
      s <- as.vector(pop.sizes)
    } else {
      p <- ploidy
    }
    if (missing(pop.sizes)) {
      warning("Argument pop.sizes is missing, proceeding with 100 individuals per pool.")
      s <- as.vector(rep(100, nPOOL))
    } else {
      if (length(pop.sizes) == nPOOL) {
        s <- as.vector(pop.sizes)
      } else {
        warning("Argument pop.sizes must be of length n where n is the number of pools, proceeding with 100 individuals per pool.")  
        s <- as.vector(rep(100, nPOOL))
      }
    }
    tmp.aux <- get.output.name(name = input, ext = "lfmm")
    sample_geno_file(input = input, output = tmp.aux, ploidy = p, sample_size = s)
    aux <- get.output.name(name = tmp.aux)
    otpt <- lfmm2pcadapt(input = tmp.aux, output = aux)
    
  } else if (class(input) %in% c("matrix", "data.frame", "array")) {
    
    if (is.null(pop.sizes)) {
      warning("Argument pop.sizes is missing, proceeding with 100 individuals per pool.")
      nPOOL <- nrow(input)
      s <- as.vector(rep(100, nPOOL))
    } else {
      nPOOL <- nrow(input)
      if (length(pop.sizes) == nPOOL){
        s <- as.vector(pop.sizes)
      } else {
        warning("Argument pop.sizes must be of length n where n is the number of pools, proceeding with 100 individuals per pool.")  
        s <- as.vector(rep(100, nPOOL))
      }
    }
    tmp <- sample_geno_matrix(freq = as.matrix(input), ploidy = p, sample_size = s)
    # tmp[which(is.na(tmp))] <- 9
    
  } else {
    stop("Input should be a file path or a matrix-like object.")
  }
}

################################################################################

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
  while (length(list.ref) > 0){
    col <- list.ref[1]
    pop.split[[idx]] <- which(pop == col)
    idx <- idx + 1
    list.ref <- list.ref[list.ref != col]
  }
  color.individuals <- array(dim = nIND)
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
      aux <- c(aux, pop[i])
      res <- c(pop[i], res)
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
get.pc = function(x, list){
  rem.na <- which(!is.na(x$zscores[list, 1]))
  v <- vector(mode = "numeric", length = length(list))
  v[rem.na] <- sapply(list[rem.na], FUN = function(h){which(x$zscores[h, ]^2 == max(x$zscores[h, ]^2, na.rm = TRUE))})
  df <- data.frame(SNP = list, PC = v)
  return(df)
}

#' Replace the file extension with the .pcadapt extension
#'
#' \code{get.output.name} returns a character string specifying the name of the output file.
#'
#' @param name a character string specifying the name of the file to be converted.
#' @param ext a character string specifying the extension of the output file.
#'
#' @examples
#' ## see also ?pcadapt for examples
#'
#' @export
#'
get.output.name = function(name, ext = "pcadapt"){
  split.name <- unlist(unlist(strsplit(name,"[.]")))  
  if (length(split.name) > 1){
    aux <- NULL
    for (k in 1:(length(split.name) - 1)){
      aux <- paste0(aux, split.name[k], ".")
    }
    if (tail(split.name, n = 1) %in% c("ped", "vcf", "lfmm")){
      aux <- paste0(aux, ext)    
    } else {
      aux <- paste0(aux, tail(split.name, n = 1), ".", ext)
    }
  } else {
    aux <- paste0(split.name[1], ".", ext) 
  }
  return(aux)
}