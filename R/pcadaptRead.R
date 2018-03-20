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
#' @param type a character string specifying the type of data to be converted 
#'   from. Converters from 'vcf' and 'ped' formats are not maintained anymore;
#'   if you have any issue with those, please use PLINK 1.9 to convert them
#'   to the 'bed' format.
#' @param type.out Either a bed file or a standard R matrix. 
#'   If the input is a matrix, then the output is automatically a matrix 
#'   (so that you don't need to specify this parameter). 
#'   If the input is a bed file, then the output is also a bed file.
#' @param allele.sep a vector of characters indicating what delimiters are used 
#'   in VCF files. By default, only "|" and "/" are recognized. 
#'   So, this argument is only useful for \code{type = "vcf"}.
#' @param local.env deprecated argument.
#' @param ploidy    deprecated argument.
#' @param pop.sizes deprecated argument.
#' @param blocksize deprecated argument.
#' 
#' @importFrom utils read.table
#'
#' @export
#'
read.pcadapt <- function(input, 
                         type = c("pcadapt", "lfmm", "vcf", "bed", "ped",
                                  "pool", "example"), 
                         ## only when input is a file path
                         type.out = c("bed", "matrix"),
                         ## only when type == "vcf"
                         allele.sep = c("/", "|"), 
                         ## deprecated arguments
                         pop.sizes, ploidy, local.env, blocksize){
  
  # In version 3.1.0, arguments 'local.env' and 'blocksize' have been removed
  if (!missing(local.env)) warning("Argument 'local.env' is deprecated.")
  if (!missing(blocksize)) warning("Argument 'blocksize' is deprecated.")
  # In version 4.0.0, arguments 'ploidy' and 'pop.sizes' have been removed
  if (!missing(pop.sizes)) warning("Argument 'pop.sizes' is deprecated.")
  if (!missing(ploidy))    warning("Argument 'ploidy' is deprecated.")
  
  type <- match.arg(type)
  
  if (class(input) == "character") {
    file2other(input, type, match.arg(type.out), match.arg(allele.sep))
  } else if (class(input) %in% c("matrix", "data.frame", "array")) {
    matrix2other(input, type)
  } else {
    stop("Input should be a file path or a matrix-like object.")
  }
}

################################################################################

check_file_size <- function(file) {
  
  if (file.size(file) > 2e9) {
    stop(paste0("It seems that the input file is quite large.\n",
                "For large 'vcf' or 'ped' files, ", 
                "please use PLINK (> 1.9) to convert them in 'bed'."))
  }
}

################################################################################

file2other <- function(input, type.in, type.out, allele.sep) {
  
  # Check if input exists  
  if ((type.in != "example") && !file.exists(input)) 
    stop(paste("File", input, "does not exist."))
  
  # File converter
  is.pcadapt <- TRUE
  if (type.in == "ped") {
    check_file_size(input)
    tmp <- tempfile(fileext = ".pcadapt")
    ped2pcadapt(input = input, output = tmp)
    input <- tmp
  } else if (type.in == "vcf") {
    check_file_size(input)
    tmp <- tempfile(fileext = ".pcadapt")
    vcf2pcadapt(input = input, output = tmp, allele.sep = allele.sep)
    input <- tmp
  } else if (type.in == "lfmm") {
    is.pcadapt <- FALSE
  } else if (type.in == "example") {
    input <- system.file("extdata", "geno3pops.pcadapt", package = "pcadapt")
  } else if (type.in == "pool") {
    return(structure(as.matrix(read.table(input)), class = "pcadapt_pool"))
  } else if (type.in == "bed") {
    p <- nrow(data.table::fread(sub("\\.bed$", ".bim", input)))
    n <- nrow(data.table::fread(sub("\\.bed$", ".fam", input)))
    return(structure(normalizePath(input), n = n, p = p, class = "pcadapt_bed"))
  }
  
  if (type.out == "matrix") {
    #mmap <- mmapcharr::mmapchar(input, code = mmapcharr:::CODE_012)
    mmap <- mmapcharr::mmapchar(input, code = CODE_012)
    structure(`if`(is.pcadapt, t(mmap[]), mmap[]), class = "pcadapt_matrix")
  } else { # bed
    writeBed(input, is.pcadapt)
  }
}

################################################################################

matrix2other <- function(input, type.in) {
  
  if (!(ncol(input) > 0) || !(nrow(input) > 0))
    stop("Invalid input genotype matrix.")

  if (type.in == "lfmm") {
    res <- as.matrix(input)
  } else if (type.in == "pcadapt") {
    res <- t(as.matrix(input))
  } else if (type.in == "pool") {
    return(structure(as.matrix(input), class = "pcadapt_pool"))
  } else {
    stop("Incorrect type.in for matrices.")
  }
  
  storage.mode(res) <- "integer"
  structure(res, class = "pcadapt_matrix")
}

################################################################################