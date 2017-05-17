output$screePlot <- plotly::renderPlotly({
  inFile <- input$file1
  K <- input$K
  nSNP <- length(r.x()$x$pvalues)
  if (is.null(inFile)){
    return(NULL)
  }
  df <- data.frame(xx = 1:K, 
                   yy = (r.x()$x$singular.values[1:K]) ^ 2 / nSNP)
  
  plotly::plot_ly(df, x = ~xx, y = ~yy,
                  mode = "lines+markers",
                  type = "scatter",
                  hoverinfo = "text") %>%
    plotly::layout(xaxis = list(title = "K", showgrid = F, autotick = FALSE,
                                dtick = 1),
                   yaxis = list(title = "Proportion of variance explained"))
  
})