################################################################################

#' Population colorization
#'
#' \code{get.score.color} allows the user to display individuals of the same
#' predefined population with the same color when using the option
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
#' If several components have the same maximum squared z-score, the first one
#' is returned. A missing value is returned only when all component z-scores
#' are missing for a marker.
#'
#' @param x an object of class `pcadapt`. 
#' @param list a list of integers corresponding to the indices of the markers of interest.
#'
#' @export
#'
get.pc <- function(x, list) {
  if (length(list) == 0L) {
    return(data.frame(SNP = list, PC = integer()))
  }

  squared.zscores <- x$zscores[list, , drop = FALSE]^2
  pc <- apply(squared.zscores, MARGIN = 1, FUN = function(zscores) {
    if (all(is.na(zscores))) {
      return(NA_integer_)
    }
    which.max(zscores)
  })

  data.frame(SNP = list, PC = as.integer(pc))
}

################################################################################
