#' pcadapt visualization tool
#'
#' \code{plot.pcadapt} is a method designed for objects of class \code{pcadapt}.
#' It provides a plotting utile for quick visualization of \code{pcadapt} 
#' objects. Different options are currently available : \code{"screeplot"}, 
#' \code{"scores"}, \code{"stat.distribution"}, \code{"manhattan"} and 
#' \code{"qqplot"}. \code{"screeplot"} shows the decay of the genotype matrix 
#' singular values and provides a figure to help with the choice of \code{K}.
#' \code{"scores"} plots the projection of the individuals onto the first two 
#' principal components. \code{"stat.distribution"} displays the histogram of 
#' the selected test statistics, as well as the estimated distribution for the 
#' neutral SNPs. \code{"manhattan"} draws the Manhattan plot of the p-values 
#' associated with the statistic of interest. \code{"qqplot"} draws a Q-Q plot 
#' of the p-values associated with the statistic of interest.
#'
#' @param x an object of class "pcadapt" generated with \code{pcadapt}.
#' @param ... \dots
#' @param option a character string specifying the figures to be displayed. If 
#' \code{NULL} (the default), all three plots are printed.
#' @param K an integer specifying the principal component of interest. \code{K} 
#' has to be specified only when using the \code{loadings} option.
#' @param i an integer indicating onto which principal component the individuals
#' are projected when the "scores" option is chosen.
#' Default value is set to \code{1}.
#' @param j an integer indicating onto which principal component the individuals 
#' are projected when the "scores" option is chosen.
#' Default value is set to \code{2}.
#' @param pop a list of integers or strings specifying which subpopulation the 
#' individuals belong to.
#' @param gg.col a list of colors to be used in the score plot.
#' @param snp.info a list containing the names of all genetic markers present in
#' the input.
#' @param chr.info a list containing the chromosome information for each marker.
#' @param chr an integer specifying the chromosome to display.
#' @param threshold for the \code{"qqplot"} option, it displays an additional 
#' bar which shows the \code{threshold} percent of SNPs with smallest p-values
#' and separates them from SNPs with higher p-values.
#' @param by.step an integer.
#' @param hline a numeric value specifying the number of standard deviations 
#' above which the z-scores are considered extreme.
#' @param plt.pkg a character string specifying the package to be used to 
#' display the graphical outputs. Use \code{"plotly"} for interactive plots, or 
#' \code{"ggplot"} for static plots.
#' 
#' @examples
#' ## see ?pcadapt for examples
#'
#' @method plot pcadapt
#'
#' @export
plot.pcadapt = function(x, 
                        ..., 
                        option = "manhattan", 
                        K = NULL, 
                        i = 1, 
                        j = 2, 
                        pop, 
                        gg.col,
                        snp.info = NULL,
                        chr.info = NULL,
                        chr = 1,
                        threshold = NULL,
                        by.step = 10,
                        hline = 3.0,
                        plt.pkg = "ggplot"
){
  if (!(option %in% c("screeplot", 
                      "scores", 
                      "manhattan", 
                      "qqplot", 
                      "stat.distribution"))){
    warning(paste("Plotting option", option, "not valid, options currently available are: screeplot, scores, manhattan, qqplot, stat.distribution."))
  } else if (attr(x, "method") == "introgression"){
    loc.anc.plotting(x, by.step = by.step, hline = hline, chr = chr)
  } else {
    if (option == "screeplot"){
      scree.plotting(x, K)
    } else if (option == "scores"){
      if (attr(x, "data.type") != "pool"){
        if (missing(pop)){
          score.plotting(x, i, j, plt.pkg = plt.pkg)
        } else {
          if (missing(gg.col)){
            score.plotting(x, i, j, pop, plt.pkg = plt.pkg)
          } else {
            score.plotting(x, i, j, pop, gg.col, plt.pkg = plt.pkg)
          }
        }
      } else {
        score.plotting(x, i, j, pop = 1:dim(x$scores)[1], plt.pkg = plt.pkg)
      }
    } else if (option == "stat.distribution"){
      if ((attr(x, "method") %in% c("mahalanobis", "communality")) == FALSE){
        if (is.null(K)){
          warning("K has to be specified.")
        } else {
          neutral.plotting(x, K)
        }
      } else {
        neutral.plotting(x, 1)
      }
    } else if (option == "manhattan"){
      if ((attr(x,"method") %in% c("mahalanobis", "communality")) == FALSE){
        if (is.null(K)){
          warning("K has to be specified.")
        } else {
          manhattan.plotting(x, K, snp.info, 
                             chr.info = chr.info, 
                             plt.pkg = plt.pkg)
        }
      } else {
        manhattan.plotting(x, K = 1, snp.info,
                           chr.info = chr.info, 
                           plt.pkg = plt.pkg)
      }
    } else if (option == "qqplot"){
      if ((attr(x, "method") %in% c("mahalanobis", "communality")) == FALSE){
        if (is.null(K)){
          warning("K has to be specified")
        } else{
          pvalqq.plotting(x, K, threshold = threshold)
        }
      } else {
        pvalqq.plotting(x, K = 1,threshold = threshold)
      }
    }
  }
}

#' Neutral Distribution Estimation
#'
#' \code{neutral.plotting} plots the histogram of the chi-squared statistics, 
#' as well as the estimated null distribution.
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
#' @importFrom ggplot2 ggplot geom_histogram geom_line theme aes_string element_rect ggtitle
#'
#' @export
neutral.plotting = function(x, K){
  idxmaf <- x$maf >= attr(x, "min.maf")
  if (attr(x, "method") == "componentwise"){
    df <- 1
    z <- x$chi2.stat[idxmaf, df]
  } else if (attr(x, "method") != "componentwise" && attr(x, "data.type") != "pool"){
    df <- attr(x, "K")
    z <- x$chi2.stat[idxmaf]
  }
  min.z <- floor(min(z[which(!is.na(z))]))
  max.z <- floor(max(z[which(!is.na(z))]) + 1)
  if (max.z > 1e5){
    stop("Can't display the histogram as the values are too high.")
  }
  xx <- seq(min.z,max.z,length = length(z))
  ggdf <- as.data.frame(cbind(xx, dchisq(xx, df = df), z))
  colnames(ggdf) <- c("abs", "ord", "chi2")
  p0 <- ggplot() + 
    geom_histogram(data = ggdf, 
                   aes_string(x = "chi2", y = "..density.."), 
                   binwidth = 0.5, 
                   fill = "#B0E2FF", 
                   alpha = 0.6, 
                   colour = "black") + 
    geom_line(data = ggdf, 
              aes_string(x = "abs", y = "ord"), 
              col = "#4F94CD", 
              size = 1) + 
    ggplot2::ggtitle("Statistics distribution")
  print(p0)
}

#' Local Ancestry
#'
#' \code{loc.anc.plotting} plots the z-scores derived from the statistics 
#' computed with `scan.intro`.
#'
#' @param x an output from \code{outlier} containing the chi-squared statistics.
#' @param by.step an integer specifying the density
#' @param hline a numeric value specifying the number of standard deviations 
#' above which the z-scores are considered extreme.
#' @param chr an integer specifying the chromosome to display. If `chr==0`,
#' all the chromosomes are concatenated.
#'
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom ggplot2 ggplot aes_string geom_area ggtitle labs geom_hline
#'
#' @export
loc.anc.plotting = function(x, by.step = by.step, hline = hline, 
                            chr = chr){
  #to display one out of ten points, thus reducing plotting time
  if (chr == 0){
    nb.chr <- length(x) / 2
    nb.stat <- 0
    wg.stat <- NULL
    pos.vline <- NULL
    for (p in 1:nb.chr){
      nb.stat <- nb.stat + length(x[[2 * p - 1]])
      if (p < nb.chr){
        pos.vline <- c(pos.vline, nb.stat)  
      }
      wg.stat <- c(wg.stat, x[[2 * p - 1]])
    }
    subset <- seq(1, nb.stat, by = by.step)
    xaxis <- (1:nb.stat)[subset]
    yaxis <- wg.stat[subset]
    sign.Y <- (yaxis > 0)
    anc.lab <- character(length = length(subset))
    anc.lab[!sign.Y] <- attr(x, "ancstrl.1")
    anc.lab[sign.Y] <- attr(x, "ancstrl.2")  
  } else {
    subset <- seq(1, length(x[[2 * chr]]), by = by.step)
    xaxis <- (x[[2 * chr]])[subset]
    yaxis <- (x[[2 * chr - 1]])[subset]
    sign.Y <- (yaxis > 0)
    anc.lab <- character(length = length(subset))
    anc.lab[!sign.Y] <- attr(x, "ancstrl.1")
    anc.lab[sign.Y] <- attr(x, "ancstrl.2")    
  }
  df <- data.frame(X = xaxis, Y = yaxis, Ancestral = anc.lab)
  plt.1 <- ggplot2::ggplot(df, aes_string("X", "Y")) +
    ggplot2::geom_area(aes_string(fill = "Ancestral")) +
    ggplot2::ggtitle(paste0("Excess of local ancestry")) +
    ggplot2::labs(x = "Position", y = "z-scores") +
    ggplot2::geom_hline(yintercept = hline) +
    ggplot2::geom_hline(yintercept = 0.0)
  if (chr == 0){
    plt.1 <- plt.1 + ggplot2::geom_vline(xintercept = pos.vline, linetype = 4)
  }
  print(plt.1)
}


#' Local Ancestry
#'
#' \code{loc.anc.plotting} plots the z-scores derived from the statistics 
#' computed with `scan.intro`.
#'
#' @param geno a scaled genotype matrix..
#' @param obj.svd an object with \code{u}, \code{d} and \code{v}.
#' @param begin an integer.
#' @param end an integer.
#' @param i an integer.
#' @param j an integer.
#' @param pop a list of integers or strings specifying which subpopulation the 
#' individuals belong to.
#'
#' @examples
#' ## see ?pcadapt for examples
#' 
#' @importFrom graphics legend plot
#' 
#' @export
#'
display.local.pca <- function(geno, obj.svd, begin = 1, end = nrow(obj.svd$v), i = 1, j = 2, pop){
  u <- cmpt_local_pca(geno, obj.svd$v, sigma = obj.svd$d, beg = begin, end = end)  
  if (missing(pop)){
    graphics::plot(u[, i], u[, j], pch = 19, 
                   xlab = paste0("PC", i), 
                   ylab = paste0("PC", j), 
                   main = paste("Window ranging from", begin, "to", end)
    )
  } else {
    n.pop <- length(unique(pop))
    plt <- grDevices::rainbow(n.pop)
    col.pop <- vector("character", length = length(pop))
    for (k in 1:n.pop){
      col.pop[which(pop == unique(pop)[k])] <- plt[which(unique(pop) == unique(pop)[k])]    
    }
    graphics::plot(u[, i], u[, j], col = col.pop, pch = 19, 
                   xlab = paste0("PC", i), 
                   ylab = paste0("PC", j),
                   main = paste("Window ranging from", begin, "to", end)
    )
    graphics::legend('bottomleft', legend = unique(pop), 
                     lty = 1, col = plt, bty = 'n', cex = .75)
  }
}

