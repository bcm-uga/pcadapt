#' Neutral Distribution Estimation
#'
#' \code{hist_plot} plots the histogram of the chi-squared statistics, as well 
#' as the estimated null distribution.
#'
#' @param x an output from \code{outlier} containing the chi-squared statistics.
#' @param K an integer indicating which principal component the histogram will 
#' be associated with.
#'
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom stats dchisq
#' @importFrom ggplot2 ggplot geom_histogram geom_line aes aes_string ggtitle
#'
#' @export
hist_plot = function(x, K) {
  
  maf.idx <- (x$maf >= attr(x, "min.maf"))
  
  if (attr(x, "method") == "componentwise") {
    k <- K
    z <- x$chi2.stat[maf.idx, k]
  } else if (attr(x, "method") != "componentwise") {
    k <- attr(x, "K")
    z <- x$chi2.stat[maf.idx]
  }
  
  min.z <- floor(min(z[which(!is.na(z))]))
  max.z <- floor(max(z[which(!is.na(z))]) + 1)
  
  if (max.z > 1e5){
    stop("Can't display the histogram as the values are too high.")
  }
  
  t <- seq(min.z, max.z, length = length(z))
  
  data.frame(abs = t,
             ord = dchisq(t, df = k), 
             chi2 = z) %>%
    ggplot() + 
    geom_histogram(aes_string(x = "chi2", 
                              y = "..density.."), 
                   binwidth = 0.5, 
                   fill = "#B0E2FF", 
                   alpha = 0.6, 
                   colour = "black") + 
    geom_line(aes_string(x = "abs", y = "ord"), 
              col = "#4F94CD", 
              size = 1) + 
    ggplot2::ggtitle("Statistics distribution")
  
}
