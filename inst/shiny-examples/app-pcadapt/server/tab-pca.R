output$pcaPlot <- plotly::renderPlotly({
  inFile <- input$file1
  inPop <- input$file2
  inSNP <- input$file3
  file <- inFile$datapath
  K <- input$K  
  
  if (is.null(inFile)){
    return(NULL)
  }
  if (is.null(r.pop()$inPop)){
    df.scores <- data.frame(xx = r.x()$x$scores[, r.ij()$i], yy = r.x()$x$scores[, r.ij()$j])
    plotly::plot_ly(df.scores, x = ~xx, y = ~yy,
                    mode = "markers",
                    type = "scatter",
                    hoverinfo = "text") %>%
      plotly::layout(xaxis = list(title = paste0("PC", r.ij()$i), showgrid = F),
                     yaxis = list(title = paste0("PC", r.ij()$j)))
  } else {
    pop <- as.character(read.table(r.pop()$inPop$datapath)[, 1])
    df.scores <- data.frame(xx = r.x()$x$scores[, r.ij()$i], yy = r.x()$x$scores[, r.ij()$j],
                            pop = pop, ind = 1:nrow(r.x()$x$scores))
    plotly::plot_ly(df.scores, x = ~xx, y = ~yy,
                    color = ~pop, text = ~paste('Ind: ', 1:nrow(r.x()$x$scores)),
                    mode = "markers",
                    type = "scatter",
                    hoverinfo = "text") %>%
      plotly::layout(xaxis = list(title = paste0("PC", r.ij()$i), showgrid = F),
                     yaxis = list(title = paste0("PC", r.ij()$j)))
  }
})
