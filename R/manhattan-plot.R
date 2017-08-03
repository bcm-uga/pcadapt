#' Manhattan Plot
#'
#' \code{manhattan.plotting} displays a Manhattan plot which represents the p-values 
#' for each SNP for a particular test statistic.
#'
#' @param x an object of class "pcadapt" generated with \code{pcadapt} 
#' containing the p-values of interest.
#' @param K an integer specifying the number of components to take into account 
#' in the scree plot.
#' @param snp.info a list containing the names of all genetic markers present in
#' the input.
#' @param chr.info a list containing the chromosome information for each marker.
#' @param plt.pkg a character string specifying the package to be used to 
#' display the graphical outputs. Use \code{"plotly"} for interactive plots, or 
#' \code{"ggplot"} for static plots.
#' 
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom ggplot2 qplot guides ggtitle
#' @importFrom plotly plot_ly layout
#'
#' @export
#'
manhattan.plotting = function(x, K, snp.info, chr.info, plt.pkg){
  if (K > attr(x, "K")){
    stop(paste0("K can't exceed ", attr(x, "K")), ".")
  }
  
  if (attr(x, "method") == "componentwise"){
    nan.idx <- !is.na(x$pvalues[, K])
    pval.K <- x$pvalues[nan.idx, K]
  } else {
    nan.idx <- !is.na(x$pvalues)
    pval.K <- x$pvalues[nan.idx]
  }
  
  if (plt.pkg == "ggplot"){
    if (!is.null(chr.info)){
      chr.int <- chr.info %% 2
      ggdf <- data.frame(x = which(nan.idx), 
                         y = -log10(pval.K),
                         chr = chr.int[nan.idx])
      res.plot <- ggplot2::ggplot(ggdf, aes_string("x", "y")) + 
        ggplot2::geom_point(aes(colour = factor(ggdf$chr))) + 
        ggplot2::guides(colour = FALSE) +
        ggplot2::scale_color_manual(values = c("#56B1F7", "#132B43"))
    } else {
      ggdf <- data.frame(x = which(nan.idx), 
                         y = -log10(pval.K))
      res.plot <- ggplot2::ggplot(ggdf, aes_string("x", "y")) +
        ggplot2::geom_point(color = "#132B43")
    }
    res.plot <- res.plot + ggplot2::ggtitle("Manhattan Plot") +
      ggplot2::labs(x = paste0("SNP (with mAF>", attr(x, "min.maf"), ")"), 
                    y = "-log10(p-values)")
    print(res.plot)
  } else if (plt.pkg == "plotly"){
    df <- data.frame(xx = (1:length(x$pvalues))[nan.idx], 
                     yy = -log10(x$pvalues)[nan.idx])
    if (is.null(snp.info)){
      txt <- (1:length(x$pvalues))[nan.idx]
    } else {
      txt <- snp.info[nan.idx]
    }
    
    if (is.null(chr.info)){
      p0 <- plotly::plot_ly(df, 
                            x = ~xx, y = ~yy,
                            text = ~paste('SNP: ', txt),
                            marker = list(color = "#132B43"),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text")
      p1 <- plotly::layout(p0, title = "Manhattan Plot", 
                           xaxis = list(title = "Index", showgrid = F),      
                           yaxis = list(title = "-log10(p-values)"),
                           showlegend = FALSE)
    } else {
      chr <- chr.info[nan.idx] %% 2
      p0 <- plotly::plot_ly(df, x = ~xx, y = ~yy,
                            text = ~paste(paste0("Chr: ", chr.info[nan.idx]), 
                                          sep = "<br>", 
                                          paste0("SNP: ", txt)),
                            color = as.character(chr),
                            colors = c("#56B1F7", "#132B43"),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text")
      p1 <- plotly::layout(p0, title = "Manhattan Plot", 
                           xaxis = list(title = "Index", showgrid = F),      
                           yaxis = list(title = "-log10(p-values)"),
                           showlegend = FALSE)
    }
    print(p1)
  }
}
