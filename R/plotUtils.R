################################################################################

utils::globalVariables(c("chi2", "density", "ord", "x", "y", "PC_i", "PC_j", "Pop"))

################################################################################

#' Principal Components Analysis Scree Plot
#'
#' \code{scree_plot} plots the scree plot associated with the principal components 
#' analysis performed on the dataset. NB : \code{pcadapt} has to be run on the
#' dataset in order to get an output readable by \code{plot.screePlot}
#'
#' @param x an output from \code{pcadapt} containing the singular values.
#' @param K an integer specifying the number of components to plot.
#'
#' @keywords internal
#' 
#' @import ggplot2
#' 
#' @export
#'
scree_plot = function(x, K = NULL) {
  
  if (is.null(K)) K <- attr(x, "K")
  
  if (K < 2) stop("The scree plot is not available for K < 2.")
  
  ggplot(data.frame(x = 1:K, y = x$singular.values[1:K]^2),
         aes(x, y)) +
    theme_bw(13) +
    geom_point() + 
    geom_line() + 
    labs(x = "PC", y = "Proportion of explained variance",
         title = paste("Scree Plot - K =", K)) +
    scale_y_log10()
}

################################################################################

#' Principal Components Analysis Scores Plot
#'
#' \code{"score_plot"} plots the projection of the individuals onto the 
#' first two principal components.
#'
#' @param x an output from \code{pcadapt} containing the scores.
#' @param i an integer indicating onto which principal component the individuals
#'  are projected when the "scores" option is chosen.
#' Default value is set to \code{1}.
#' @param j an integer indicating onto which principal component the individuals 
#' are projected when the "scores" option is chosen.
#' Default value is set to \code{2}.
#' @param pop a list of integers or strings specifying which subpopulation the 
#' individuals belong to.
#' @param plt.pkg a character string specifying the package to be used to 
#' display the graphical outputs. Use \code{"plotly"} for interactive plots, or 
#' \code{"ggplot"} for static plots.
#' 
#' @keywords internal
#'
#' @import ggplot2
#' @importFrom magrittr "%>%"
#'
#' @export
#'
score_plot = function(x, i = 1, j = 2, pop, col, plt.pkg = "ggplot") {
  
  if (attr(x, "K") == 1) j <- 1
  
  if (i > attr(x, "K"))
    stop(paste0("i can't exceed ", attr(x, "K"), "."))
  
  if (j > attr(x, "K"))
    stop(paste0("j can't exceed ", attr(x, "K"), "."))
  
  df <- data.frame(PC_i = x$scores[, i], PC_j = x$scores[, j])    
  
  if (!missing(pop)) df$Pop <- pop
  
  if (plt.pkg == "ggplot") {
    
    res.plot <- ggplot(df, aes(PC_i, PC_j)) + 
      theme_bw(13) +
      geom_point() + 
      labs(x = paste0("PC", i), y = paste0("PC", j),
           title = paste0("Projection onto PC", i, " and PC", j))
    
    if (!missing(pop)) {
      res.plot <- res.plot + 
        geom_point(aes(color = factor(Pop)))
      if (missing(col)) {
        res.plot <- res.plot + scale_color_hue(name = "")
      } else {
        if (length(col) < length(unique(pop))) {
          pers.col <- c(col, rainbow(length(unique(pop)) - length(col)))
        } else if (length(col) == length(unique(pop))) {
          pers.col <- col
        } else if (length(col) > length(unique(pop))) {
          pers.col <- col[1:length(unique(pop))]
        }
        res.plot <- res.plot + 
          scale_color_manual(name = "", values = pers.col)
      }
    }
    
    res.plot
    
  } else if (plt.pkg == "plotly") {
    
    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Please install package 'plotly'.")
    
    if (missing(pop)) {
      plotly::plot_ly(df, 
                      x = ~PC_i, 
                      y = ~PC_j, 
                      text = ~paste('Ind:', 1:nrow(x$scores)),
                      mode = "markers", 
                      type = "scatter", 
                      hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = FALSE),      
                       yaxis = list(title = paste0("PC", j)))
    } else {
      plotly::plot_ly(df, 
                      x = ~PC_i, 
                      y = ~PC_j, 
                      color = factor(pop), 
                      colors = col,
                      text = ~paste('Ind:', 1:nrow(x$scores)),
                      mode = "markers", 
                      type = "scatter", 
                      hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = FALSE),      
                       yaxis = list(title = paste0("PC", j)))  
    }
    
  } else {
    stop('`plt.pkg` should be either "ggplot" or "plotly"')
  }
}

################################################################################

#' Neutral Distribution Estimation
#'
#' \code{hist_plot} plots the histogram of the chi-squared statistics, as well 
#' as the estimated null distribution.
#'
#' @param x an output from \code{outlier} containing the chi-squared statistics.
#' @param K an integer indicating which principal component the histogram will 
#' be associated with.
#' 
#' @keywords internal
#'
#' @import ggplot2
#' @importFrom stats density
#'
#' @export
#' 
hist_plot = function(x, K) {
  
  maf.idx <- (x$maf >= attr(x, "min.maf"))
  
  if (attr(x, "method") == "componentwise") {
    k <- K
    z <- x$chi2.stat[maf.idx, k]
  } else {
    k <- attr(x, "K")
    z <- x$chi2.stat[maf.idx]
  }
  
  min.z <- floor  (min(z, na.rm = TRUE))
  max.z <- ceiling(max(z, na.rm = TRUE))
  
  if (max.z > 1e5)
    stop("Can't display the histogram as the values are too high.")
  
  t <- seq(min.z, max.z, length = length(z))
  
  ggplot(data.frame(abs = t, ord = stats::dchisq(t, df = k), chi2 = z)) + 
    theme_bw(13) +
    geom_histogram(aes(chi2, after_stat(density)), binwidth = 0.5, 
                   fill = "#B0E2FF", alpha = 0.6, colour = "black") + 
    geom_line(aes(abs, ord), col = "#4F94CD", linewidth = 1) + 
    ggtitle("Statistics distribution")
}

################################################################################

#' Manhattan Plot
#'
#' \code{manhattan_plot} displays a Manhattan plot which represents the p-values 
#' for each SNP for a particular test statistic.
#'
#' @param x an object of class "pcadapt" generated with \code{pcadapt} 
#' containing the p-values of interest.
#' @param chr.info a list containing the chromosome information for each marker.
#' @param snp.info a list containing the names of all genetic markers present in
#' the input.
#' @param K an integer specifying which principal component to display when 
#' \code{method="componentwise"}.
#' @param plt.pkg a character string specifying the package to be used to 
#' display the graphical outputs. Use \code{"plotly"} for interactive plots, or 
#' \code{"ggplot"} for static plots.
#' 
#' @keywords internal
#'
#' @import ggplot2
#'
#' @export
#'
manhattan_plot = function(x, chr.info, snp.info, plt.pkg = "ggplot", K = 1) {
  
  if (attr(x, "method") == "componentwise") {
    if (K > attr(x, "K")) {
      stop(paste0("K can't exceed ", attr(x, "K")), ".")
    }
    notNA.idx <- !is.na(x$pvalues[, K])
  } else {
    notNA.idx <- !is.na(x$pvalues)
  }
  
  df <- data.frame(x = which(notNA.idx), 
                   y = -pchisq(x$chi2.stat[notNA.idx], df = attr(x, "K"),
                               lower.tail = FALSE, log.p = TRUE) / log(10))
  
  if (plt.pkg == "ggplot") {
    
    if (!is.null(chr.info)) {
      chr.int <- chr.info %% 2
      df$chr = chr.int[notNA.idx]
      res.plot <- ggplot(df, aes(x, y)) + 
        geom_point(aes(color = factor(chr))) + 
        guides(color = "none") +
        scale_color_manual(values = c("grey10", "grey60"))
    } else {
      res.plot <- ggplot(df, aes(x, y)) + geom_point()
    }
    res.plot + theme_bw(13) +
      labs(x = paste0("SNP (with MAF>", attr(x, "min.maf"), ")"), 
           y = "-log10(p-values)", title = "Manhattan Plot")
    
  } else if (plt.pkg == "plotly") {
    
    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Please install package 'plotly'.")
    
    txt <- if (is.null(snp.info)) df$x else snp.info[notNA.idx]
    
    if (is.null(chr.info)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~x, y = ~y,
                            text = ~paste('SNP:', txt),
                            marker = list(color = "#132B43"),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text")
      p1 <- plotly::layout(p0, title = "Manhattan Plot", 
                           xaxis = list(title = "Index", showgrid = FALSE),      
                           yaxis = list(title = "-log10(p-values)"),
                           showlegend = FALSE)
    } else {
      chr <- chr.info[notNA.idx] %% 2
      p0 <- plotly::plot_ly(df, x = ~x, y = ~y,
                            text = ~paste(paste0("Chr: ", chr.info[notNA.idx]), 
                                          sep = "<br>", 
                                          paste0("SNP: ", txt)),
                            color = as.character(chr),
                            colors = c("#56B1F7", "#132B43"),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text")
      p1 <- plotly::layout(p0, title = "Manhattan Plot", 
                           xaxis = list(title = "Index", showgrid = FALSE),      
                           yaxis = list(title = "-log10(p-values)"),
                           showlegend = FALSE)
    }
    p1
  }
}

################################################################################

#' p-values Q-Q Plot
#'
#' \code{qq_plot} plots a Q-Q plot of the p-values computed.
#'
#' @param x an output from \code{outlier} containing the p-values of interest.
#' @param K an integer specifying which principal component to display when \code{method="componentwise"}.
#'
#' @keywords internal
#'
#' @import ggplot2
#'
#' @export
#'
qq_plot = function(x, K = 1) {
  
  if (attr(x, "method") == "componentwise") {
    if (K > attr(x, "K")) {
      stop(paste0("K can't exceed ", attr(x, "K")), ".")
    }
    notNA.idx <- !is.na(x$pvalues[, K])
  } else {
    notNA.idx <- !is.na(x$pvalues)
  }
  
  lpval <- -pchisq(x$chi2.stat[notNA.idx], df = attr(x, "K"),
                   lower.tail = FALSE, log.p = TRUE) / log(10)
  sorted.lpval <- sort(lpval, decreasing = TRUE)
  expected.p <- stats::ppoints(length(lpval))
  
  ggplot(data.frame(x = -log10(expected.p), y = sorted.lpval)) +
    theme_bw(13) +
    geom_point(aes(x, y)) +
    geom_abline(color = "red") +
    labs(x = "Expected -log10(p-values)",
         y = "Observed -log10(p-values)", title = "Q-Q plot")
}

################################################################################

#' pcadapt visualization tool
#'
#' \code{plot.pcadapt} is a method designed for objects of class \code{pcadapt}.
#' It provides plotting options for quick visualization of \code{pcadapt} 
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
#' @param i an integer indicating onto which principal component the individuals
#' are projected when the "scores" option is chosen.
#' Default value is set to \code{1}.
#' @param j an integer indicating onto which principal component the individuals 
#' are projected when the "scores" option is chosen.
#' Default value is set to \code{2}.
#' @param pop a list of integers or strings specifying which subpopulation the 
#' individuals belong to.
#' @param col a list of colors to be used in the score plot.
#' @param chr.info a list containing the chromosome information for each marker.
#' @param snp.info a list containing the names of all genetic markers present in
#' the input.
#' @param plt.pkg a character string specifying the package to be used to 
#' display the graphical outputs. Use \code{"plotly"} for interactive plots, or 
#' \code{"ggplot"} for static plots.
#' @param K an integer specifying the principal component of interest. \code{K} 
#' has to be specified only when using the \code{"componentwise"} method.
#' 
#' @examples
#' ## see ?pcadapt for examples
#'
#' @method plot pcadapt
#'
#' @export
plot.pcadapt = function(x, 
                        ..., 
                        option = c("manhattan", "screeplot", "scores", 
                                   "qqplot", "stat.distribution"), 
                        i = 1, 
                        j = 2, 
                        pop, 
                        col,
                        chr.info = NULL,
                        snp.info = NULL,
                        plt.pkg = "ggplot",
                        K = NULL) {
  
  option <- match.arg(option)
  
  if (option == "screeplot") {
    scree_plot(x, K)
  } else if (option == "scores") {
    if (missing(pop)) {
      score_plot(x, i, j, plt.pkg = plt.pkg)
    } else {
      if (missing(col)) {
        score_plot(x, i, j, pop, plt.pkg = plt.pkg)
      } else {
        score_plot(x, i, j, pop, col, plt.pkg = plt.pkg)
      }
    }
    
  } else if (option == "stat.distribution") {
    if (attr(x, "method") == "componentwise") {
      if (is.null(K)) {
        stop("K has to be specified.")
      } else {
        hist_plot(x, K)
      }
    } else {
      hist_plot(x)
    }
  } else if (option == "manhattan") {
    if (attr(x, "method") == "componentwise") {
      if (is.null(K)) {
        stop("K has to be specified.")
      } else {
        manhattan_plot(
          x,
          chr.info = chr.info,
          snp.info = snp.info,
          plt.pkg = plt.pkg,
          K = K
        )
      }
    } else {
      manhattan_plot(x,
                     chr.info = chr.info,
                     snp.info = snp.info,
                     plt.pkg = plt.pkg)
    }
  } else if (option == "qqplot") {
    if (attr(x, "method") == "componentwise") {
      if (is.null(K)) {
        stop("K has to be specified")
      } else{
        qq_plot(x, K)
      }
    } else {
      qq_plot(x)
    }
  }
}

################################################################################
