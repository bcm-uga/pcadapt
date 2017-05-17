output$histPlot <- plotly::renderPlotly({
  inFile <- input$file1
  file <- inFile$datapath
  
  if (is.null(inFile)){
    return(NULL)
  }
  
  df <- data.frame(xx = r.x()$x$pvalues)
  plotly::plot_ly(df, x = ~xx, type = "histogram") %>%
    plotly::layout(xaxis = list(title = "p-values"),
                   yaxis = list(title = "Counts"))
})