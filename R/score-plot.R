#' Principal Components Analysis Scores Plot
#'
#' \code{"score.plotting"} plots the projection of the individuals onto the 
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
#'
#' @export
#'
score.plotting = function(x, i = 1, j = 2, pop, gg.col, plt.pkg = "ggplot"){
  
  if (attr(x, "K") == 1){
    warning("K = 1, option not available since two principal components have to be computed at least.")
  } else {
    
    if (i > attr(x, "K")){
      stop(paste0("i can't exceed ", attr(x, "K"), "."))
    }
    
    if (j > attr(x, "K")){
      stop(paste0("j can't exceed ", attr(x, "K"), "."))
    }
    
    ggdf <- data.frame(PC_i = x$scores[, i], PC_j = x$scores[, j])  
    
    if (plt.pkg == "ggplot"){
      res.plot <- ggplot2::ggplot(ggdf, aes_string("PC_i", "PC_j")) + 
        ggplot2::geom_point() + 
        ggplot2::ggtitle(paste0("Projection onto PC", i, " and PC", j)) +
        ggplot2::labs(x = paste0("PC", i), y = paste0("PC", j))
      
      if (!missing(pop)){
        pop.to.int <- get.score.color(pop)
        popnames <- get.pop.names(pop)
        ggdf$Pop <- pop.to.int
        res.plot <- res.plot + ggplot2::geom_point(aes(colour = factor(ggdf$Pop)))
        if (missing(gg.col)){
          res.plot <- res.plot + ggplot2::scale_color_hue(name = " ", labels = popnames)
        } else {
          if (length(gg.col) < length(popnames)){
            pers.col <- c(gg.col, rainbow(length(popnames) - length(gg.col)))
          } else if (length(gg.col) == length(popnames)){
            pers.col <- gg.col
          } else if (length(gg.col) > length(popnames)){
            pers.col <- gg.col[1:length(popnames)]
          }
          res.plot <- res.plot + ggplot2::scale_color_manual(name = " ", labels = popnames, values = pers.col)
        }
      }
      print(res.plot)
    } else if (plt.pkg == "plotly"){
      if (missing(pop)){
        p0 <- plotly::plot_ly(ggdf, x = ~PC_i, y = ~PC_j, 
                              text = ~paste('Ind: ', 1:nrow(x$scores)),
                              mode = "markers", 
                              type = "scatter", 
                              hoverinfo = "text")
        p1 <- plotly::layout(p0, title = paste0("Projection onto PC", i, " and PC", j), 
                             xaxis = list(title = paste0("PC", i), showgrid = F),      
                             yaxis = list(title = paste0("PC", j)))
      } else if (!missing(pop)){
        p0 <- plotly::plot_ly(ggdf, x = ~PC_i, y = ~PC_j, 
                              color = pop, text = ~paste('Ind: ', 1:nrow(x$scores)),
                              mode = "markers", 
                              type = "scatter", 
                              hoverinfo = "text")
        p1 <- plotly::layout(p0, title = paste0("Projection onto PC", i, " and PC", j), 
                             xaxis = list(title = paste0("PC", i), showgrid = F),      
                             yaxis = list(title = paste0("PC", j)))  
      }
    }
  }
}