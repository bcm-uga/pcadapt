################################################################################

#' File Converter
#'
#' \code{read.pcadapt} converts genotype matrices or files to an appropriate
#'   format readable by \code{pcadapt}. For a file as input, you can return 
#'   either a matrix or convert it in bed/bim/fam files. For a matrix as input,
#'   this return a matrix. 
#'
#' @param input a genotype matrix or a character string specifying the name of 
#'   the file to be converted.
#' @param type.in a character string specifying the type of data to be converted 
#'   from. Supported formats are: 
#' @param local.env deprecated argument.
#' @param ploidy an integer specifying the ploidy of the individuals. 
#'   Default is 2.
#' @param pop.sizes a vector specifying the number of individuals to be sampled 
#'   for each pool.
#' @param allele.sep a vector of characters indicating what delimiters are used 
#'   in VCF files. By default, only "|" and "/" are recognized. So, this argument
#'   is only useful for \code{type.in = "vcf"}.
#' @param blocksize deprecated argument.
#'
#' @export
#'
read.pcadapt <- function(input, 
                         type.in = c("pcadapt", "lfmm", "vcf", 
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
  
  
  if (class(input) == "character") {
    file2other(input, type.in, type.out, pop.sizes, allele.sep, ploidy)
  } else if (class(input) %in% c("matrix", "data.frame", "array")) {
    matrix2other(input, type.in, type.out, pop.sizes, allele.sep, ploidy)
  } else {
    stop("Input should be a file path or a matrix-like object.")
  }
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

file2other <- function(input, type.in, type.out,
                       pop.sizes, allele.sep, ploidy) {
  
  ## Check if input exists ##    
  if ((type.in != "example") && !file.exists(input)) 
    stop(paste("File", input, "does not exist."))
  
  ## File converter ##
  is.pcadapt <- TRUE
  if (type.in == "ped") {
    check_file_size(input)
    tmp <- tempfile(fileext = ".pcadapt")
    ped2pcadapt(input = input, output = tmp)
    bedfile <- paste0(input, ".bed")
    input <- tmp
  } else if (type.in == "vcf") {
    check_file_size(input)
    tmp <- tempfile(fileext = ".pcadapt")
    vcf2pcadapt(input = input, output = tmp, allele.sep = allele.sep)
    bedfile <- paste0(input, ".bed")
    input <- tmp
  } else if (type.in == "lfmm") {
    is.pcadapt <- FALSE
  } else if (type.in == "example") {
    input <- system.file("extdata", "geno3pops.pcadapt", package = "pcadapt")
  } else if (type.in == "pool") {
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
    
    tmp <- tempfile(fileext = "lfmm")
    sample_geno_file(input = input, output = tmp, ploidy = p, sample_size = s)
    is.pcadapt <- FALSE
    bedfile <- paste0(input, ".bed")
    input <- tmp
  }
  
  if (type.out == "matrix") {
    mmap <- mmapcharr::mmapchar(input, code = mmapcharr:::CODE_012)
    `if`(is.pcadapt, t(mmap[]), mmap[])
  } else { # bed
    writeBed(input, is.pcadapt, bedfile)
  }
}

################################################################################

matrix2other <- function(input, type.in, type.out,
                         pop.sizes, allele.sep, ploidy) {
  
  if (!(ncol(input) > 0) || !(nrow(input) > 0))
    stop("Invalid input genotype matrix.")
  
  if (type.out != "matrix")
    stop("If the input is a matrix, then only an output matrix is supported.")
  
  if (type.in == "lfmm") {
    res <- as.matrix(input)
  } else if (type.in == "pcadapt") {
    res <- t(as.matrix(input))
  } else if (type.in == "pool") {
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
    res <- sample_geno_matrix(freq = as.matrix(input), ploidy = p, sample_size = s)
  } else {
    stop("Incorrect type.in for matrices.")
  }
  
  storage.mode(res) <- "integer"
  res
}

################################################################################