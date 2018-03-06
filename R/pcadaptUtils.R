################################################################################

CODE_012 <- rep(NA_integer_, 256); CODE_012[49:51] <- 0:2

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

################################################################################

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

################################################################################

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

################################################################################