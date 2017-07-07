#' Genotype matrix imputation
#'
#' \code{impute.pcadapt} imputes values based on medians.
#'
#' @param input a genotype matrix or a character string specifying the name of 
#' the file to be imputed.
#' @param pop a vector of integers or strings specifying which subpopulation the 
#' individuals belong to.
#' @param skip.return a logical value specifying whether the list of markers to 
#' be skipped should be returned or not. 
#' 
#' @return The returned value is a list containing the test statistics and the 
#' associated p-values.
#' 
#' @export
#'
impute.pcadapt = function(input, pop, skip.return = FALSE){
  if (is.character(input)){
    dt <- as.matrix(data.table::fread(input))
  } else if (class(input) %in% c("array", "matrix", "data.frame")){
    dt <- as.matrix(input)
  } else {
    stop("Wrong argument.")
  }
  
  if (missing(pop)){
    y <- impute_geno(dt)   
  } else if (!missing(pop)){
    pop.int <- assign.int.labels(pop)
    pop.names <- pcadapt::get.pop.names(pop.int)
    y <- impute_geno_pop(dt, pop.int, pop.names)   
  }
  if (skip.return == FALSE){
    return(list(x = y$x[y$skip == 0, ]))  
  } else if (skip.return == TRUE){
    return(list(x = y$x[y$skip == 0, ], ix = which(y$skip == 1))) 
  }
}

