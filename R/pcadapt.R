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
#' @param input a genotype matrix or a character string specifying the name of the file to be processed with \code{pcadapt}.
#' @param K an integer specifying the number of principal components to retain.
#' @param method a character string specifying the method to be used to compute
#' the p-values. Three statistics are currently available, \code{"mahalanobis"},
#' \code{"communality"} and \code{"componentwise"}.
#' @param data.type a character string specifying the type of data being read, either a \code{genotype} matrix (\code{data.type="genotype"}),
#' or a matrix of allele frequencies (\code{data.type="pool"}).
#' @param min.maf a value between \code{0} and \code{0.45} specifying the threshold of minor allele frequencies above which p-values are computed.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' @param output.filename deprecated argument.
#' @param clean.files deprecated argument.
#' @param transpose deprecated argument.
#' 
#' @return The returned value \code{x} is an object of class \code{pcadapt}.
#' 
#' @useDynLib pcadapt
#'
#' @export
#'
pcadapt = function(input, 
                   K = 5, 
                   method = "mahalanobis", 
                   data.type = "genotype",
                   min.maf = 0.05, 
                   ploidy = 2,
                   output.filename,
                   clean.files,
                   transpose){
  
  #############################################
  ########## test arguments and init ##########
  #############################################
  
  if (missing(input)){
    appDir = system.file("shiny-examples/app-pcadapt", package = "pcadapt")
    if (appDir == "") {
      stop("Could not find Shiny app in pcadapt.", call. = FALSE)
    }
    shiny::runApp(appDir, display.mode = "normal")  
  } else {
    
    ## In version 3.1.0, argument output.filename has been removed ##
    if (!missing(output.filename)){
      warning("Argument output.filename is deprecated. Please refer to the latest vignette for further information.")
    }
    
    ## In version 3.1.0, argument clean.files has been removed ##
    if (!missing(clean.files)){
      warning("Argument clean.files is deprecated. Please refer to the latest vignette for further information.")
    }
    
    ## In version 3.0.3, argument transpose has been removed ##
    if (!missing(transpose)){
      stop("Argument transpose is deprecated. Please refer to the latest vignette for further information.")
    }
    
    if (data.type == "genotype"){
      if (!(class(K) %in% c("numeric", "integer")) || K <= 0){
        stop("K has to be a positive integer.")
      }
      
      if (!(method %in% c("mahalanobis", "communality", "componentwise"))){
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
      
      local <- NULL
      if (is.character(input) && !file.exists(input)){
        stop(paste0("File ", input, " does not exist."))
      } else if (is.character(input) && file.exists(input)){
        local <- FALSE
      }
      
      if ((class(input) %in% c("matrix", "data.frame", "array"))){
        local <- TRUE  
      } 
      
      if (!is.null(local)){
        if (local == FALSE){
          obj.pca <- create.pcadapt(input = input, 
                                    K = K,
                                    method = method, 
                                    min.maf = min.maf, 
                                    ploidy = ploidy, 
                                    type = 0)
        } else if (local == TRUE){
          obj.pca <- create.pcadapt(input = as.matrix(input), 
                                    K = K, 
                                    method = method, 
                                    min.maf = min.maf, 
                                    ploidy = ploidy, 
                                    type = 1)
        }
        class(obj.pca) <- 'pcadapt'
        attr(obj.pca, "K") <- K
        attr(obj.pca, "data.type") <- "genotype"
        attr(obj.pca, "method") <- method
        attr(obj.pca, "min.maf") <- min.maf
        return(obj.pca)
      } else {
        stop("Input class not supported.")
      }
    } else if (data.type == "pool"){
      stop('Option data.type = "pool" is deprecated. Use the read.pcadapt function instead. Usage:\n 
         geno <- read.pcadapt(input, type = "pool", local.env = TRUE)\n
         x <- pcadapt(input = geno, K = ...)')
      stop('')
    }
  }
}

# Shiny app
#
# \code{pcadapt} comes with a Shiny interface.
#
#' @export
run.pcadapt <- function() {
  appDir <- system.file("shiny-examples", "app-pcadapt", package = "pcadapt")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `pcadapt`.", call. = FALSE)
  }
  shiny::runApp(appDir, launch.browser = TRUE)
}
