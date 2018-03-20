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
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom ggplot2 ggplot ggtitle labs geom_point guides aes aes_string geom_point scale_color_hue
#' @importFrom magrittr "%>%"
#'
#' @export
#'
score_plot = function(x, i = 1, j = 2, pop, col, plt.pkg = "ggplot"){
  
  if (attr(x, "K") == 1) {
    j <- 1
  }
  
  if (i > attr(x, "K")) {
    stop(paste0("i can't exceed ", attr(x, "K"), "."))
  }
  
  if (j > attr(x, "K")) {
    stop(paste0("j can't exceed ", attr(x, "K"), "."))
  }
  
  df <- data.frame(PC_i = x$scores[, i], PC_j = x$scores[, j])    
  
  if (!missing(pop)) {
    df$Pop <- pop
  }
  
  if (plt.pkg == "ggplot") {
    res.plot <- ggplot2::ggplot(df, aes_string("PC_i", "PC_j")) + 
      ggplot2::geom_point() + 
      ggplot2::ggtitle(paste0("Projection onto PC", i, " and PC", j)) +
      ggplot2::labs(x = paste0("PC", i), y = paste0("PC", j))
    
    if (!missing(pop)) {
      res.plot <- res.plot + 
        ggplot2::geom_point(aes(colour = factor(df$Pop)))
      if (missing(col)) {
        res.plot <- res.plot + ggplot2::scale_color_hue(name = "")
      } else {
        if (length(col) < length(unique(pop))) {
          pers.col <- c(col, rainbow(length(unique(pop)) - length(col)))
        } else if (length(col) == length(unique(pop))) {
          pers.col <- col
        } else if (length(col) > length(unique(pop))) {
          pers.col <- col[1:length(unique(pop))]
        }
        res.plot <- res.plot + 
          ggplot2::scale_color_manual(name = "", values = pers.col)
      }
    }
    print(res.plot)
  } else if (plt.pkg == "plotly") {
    if (missing(pop)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~PC_i, 
                            y = ~PC_j, 
                            text = ~paste('Ind: ', 1:nrow(x$scores)),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = F),      
                       yaxis = list(title = paste0("PC", j)))
    } else if (!missing(pop)) {
      p0 <- plotly::plot_ly(df, 
                            x = ~PC_i, 
                            y = ~PC_j, 
                            color = pop, 
                            text = ~paste('Ind: ', 1:nrow(x$scores)),
                            mode = "markers", 
                            type = "scatter", 
                            hoverinfo = "text") %>%
        plotly::layout(title = paste0("Projection onto PC", i, " and PC", j), 
                       xaxis = list(title = paste0("PC", i), showgrid = F),      
                       yaxis = list(title = paste0("PC", j)))  
    }
    print(p0)
  }
}