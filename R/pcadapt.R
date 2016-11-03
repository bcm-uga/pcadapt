#' Principal Component Analysis for outlier detection
#'
#' \code{pcadapt} performs principal component analysis and computes p-values to test for outliers. The test for
#' outliers is based on the correlations between genetic variation and the first \code{K} principal components.
#' \code{pcadapt} also handles Pool-seq data for which the statistical analysis is
#' performed on the genetic markers frequencies. Returns an object of class \code{pcadapt}.
#'
#' @details First, a principal component analysis is performed on the scaled and centered genotype data. To account for missing
#' data, the correlation matrix between individuals is computed using only the markers available for each
#' pair of individuals. Depending on the specified \code{method}, different test statistics can be used.
#'
#' \code{mahalanobis} (default): the robust Mahalanobis distance is computed for each genetic marker using a robust
#' estimate of both mean and covariance matrix between the \code{K} vectors of z-scores.
#'
#' \code{communality}: the communality statistic measures the proportion of variance explained by the first \code{K} PCs.
#'
#' \code{componentwise}: returns a matrix of z-scores.
#'
#' To compute p-values, test statistics (\code{stat}) are divided by a genomic inflation factor (\code{gif}) when \code{method="mahalanobis"}.
#' When \code{method="communality"}, the test statistic is first multiplied by \code{K} and divided by the percentage of variance explained by the first \code{K} PCs
#' before accounting for genomic inflation factor. When using \code{method="mahalanobis"} or \code{"communality"}, the scaled statistics (\code{chi2_stat}) should follow
#' a chi-squared distribution with \code{K} degrees of freedom. When using \code{method="componentwise"}, the z-scores should follow a chi-squared distribution with \code{1}
#' degree of freedom. For Pool-seq data, \code{pcadapt} provides p-values based on the Mahalanobis distance for each SNP.
#'
#' @param input a character string specifying the name of the file to be processed with \code{pcadapt}.
#' @param K an integer specifying the number of principal components to retain.
#' @param method a character string specifying the method to be used to compute
#' the p-values. Three statistics are currently available, \code{"mahalanobis"},
#' \code{"communality"} and \code{"componentwise"}.
#' @param data.type a character string specifying the type of data being read, either a \code{genotype} matrix (\code{data.type="genotype"}),
#' or a matrix of allele frequencies (\code{data.type="pool"}).
#' @param min.maf a value between \code{0} and \code{0.45} specifying the threshold of minor allele frequencies above which p-values are computed.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' @param output.filename a character string specifying the names of the files created by \code{pcadapt}.
#' @param clean.files a logical value indicating whether the auxiliary files should be deleted or not.
#' @param transpose deprecated argument.
#' @param cover.matrix a matrix specifying the average coverage per genetic marker and per population.
#' @return  The returned value \code{x} is an object of class \code{pcadapt}.
#' 
#' @importFrom utils read.table write.table
#'
#' @useDynLib pcadapt wrapper_pcadapt
#'
#' @export
#'
pcadapt <- function(input,
                        K=5,
                        method="mahalanobis",
                        data.type="genotype",
                        min.maf=0.05,
                        ploidy=2,
                        output.filename="pcadapt_output",
                        clean.files=TRUE,
                        transpose,
                        cover.matrix=NULL){

  #############################################
  ########## test arguments and init ##########
  #############################################
  
  ## In version 3.0.3, argument transpose has been removed ##
  if (!missing(transpose)){
    stop("Argument transpose is deprecated. Please refer to the latest vignette for further information.")
  }
  
  if (data.type == "genotype"){
    input.filename <- 1
    # input file
    if ((class(input) == "character") && (!file.exists(input))){
      stop(paste0("File ",input," does not exist."))
    } else if ((class(input) == "character") && (file.exists(input))){
      input.filename <- input
    }
    
    if (class(input.filename) != "character" || (!file.exists(input.filename))){
      stop("Invalid argument. Make sure the file exists or that the data is in the workspace.")
    }
    
    if (!(class(K) %in% c("numeric","integer")) || K <= 0){
      stop("K has to be a positive integer.")
    }
    
    if (!(method %in% c("mahalanobis","communality","componentwise"))){
      warning("Unknown method. 'mahalanobis' will be used hence.")
      method <- "mahalanobis"
    }
    
    if (class(min.maf) != "numeric" || min.maf < 0 || min.maf > 0.45){
      warning("min.maf has to be a real number between 0 and 0.45. Default value will be used hence.")
      min.maf <- 0.05
    }
    
    if (!(ploidy %in% c(1,2))){
      stop("pcadapt only supports haploid and diploid data.")
    }
    
    if (class(output.filename) != "character"){
      warning("output.filename has to be a character string. Default value will be used hence.")
      output.filename <- "pcadapt_output"
    }
    
    #############################################
    ######### use the wrapped functions #########
    #############################################
    
    if (ploidy == 2){
      .C("wrapper_pcadapt",
         as.character(input.filename),
         as.integer(K),
         as.double(min.maf),
         as.integer(0),
         as.character(output.filename),
         PACKAGE = "pcadapt"
      );
    } else if (ploidy == 1){
      .C("wrapper_pcadapt",
         as.character(input.filename),
         as.integer(K),
         as.double(min.maf),
         as.integer(1),
         as.character(output.filename),
         PACKAGE = "pcadapt"
      );
    }
    
    res <- create.pcadapt(output.filename,K,method,data.type,min.maf)
    
    if (clean.files == TRUE){
      # Clean files
      # if (file.exists(input.filename)){
      #   file.remove(input.filename)
      # }
      file.remove(paste0(output.filename,".loadings"))
      file.remove(paste0(output.filename,".scores"))
      file.remove(paste0(output.filename,".zscores"))
      file.remove(paste0(output.filename,".maf"))
      file.remove(paste0(output.filename,".sigma"))
    }
  } else if (data.type == "pool"){
    if (class(input) %in% c("array","matrix","data.frame")){
      data <- input  
    } else if (class(input) == "character"){
      data <- read.table(input)  
    } else {
      stop("Invalid input argument.")
    }
    if (missing(K)){
      K <- nrow(data)-1
    }
    if (missing(cover.matrix)){
      cov.mat <- array(1,dim=dim(data))
      res <- create.pcadapt.pool(data,K=K,min.maf=min.maf,cover.matrix=cov.mat)
    } else {
      res <- create.pcadapt.pool(data,K=K,min.maf=min.maf,cover.matrix=cover.matrix)
    }
    
   }
  return(res)
}




