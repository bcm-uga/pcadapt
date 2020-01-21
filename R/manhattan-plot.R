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
  
  if (plt.pkg == "ggplot") {
    if (!is.null(chr.info)) {
      chr.int <- chr.info %% 2
      df <- data.frame(x = which(notNA.idx), 
                       y = -as.numeric(pchisq(x$chi2.stat[notNA.idx], 
                                              df = attr(x, "K"), 
                                              lower.tail = FALSE, 
                                              log.p = TRUE) / log(10)),
                       chr = chr.int[notNA.idx])
      res.plot <- ggplot(df, aes_string("x", "y")) + 
        geom_point(aes(colour = factor(df$chr))) + 
        guides(colour = FALSE) +
        scale_color_manual(values = c("black", "grey"))
    } else {
      df <- data.frame(x = which(notNA.idx), 
                       y = -as.numeric(pchisq(x$chi2.stat[notNA.idx], 
                                              df = attr(x, "K"), 
                                              lower.tail = FALSE, 
                                              log.p = TRUE) / log(10)))
      res.plot <- ggplot(df, aes_string("x", "y")) + geom_point()
    }
    res.plot <- res.plot + ggtitle("Manhattan Plot") +
      labs(x = paste0("SNP (with mAF>", attr(x, "min.maf"), ")"), 
           y = "-log10(p-values)")
    print(res.plot)
  } else if (plt.pkg == "plotly") {
    
    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Please install package 'plotly'.")
    
    df <- data.frame(xx = seq_along(x$pvalues)[notNA.idx], 
                     yy = -as.numeric(pchisq(x$chi2.stat[notNA.idx], 
                                             df = attr(x, "K"), 
                                             lower.tail = FALSE, 
                                             log.p = TRUE) / log(10))[notNA.idx])
    if (is.null(snp.info)) {
      txt <- seq_along(x$pvalues)[notNA.idx]
    } else {
      txt <- snp.info[notNA.idx]
    }
    
    if (is.null(chr.info)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~xx, y = ~yy,
                            text = ~paste('SNP: ', txt),
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
      p0 <- plotly::plot_ly(df, x = ~xx, y = ~yy,
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
    print(p1)
  }
}
