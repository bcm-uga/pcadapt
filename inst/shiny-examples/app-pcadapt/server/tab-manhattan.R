#output$distPlot <- plotly::renderPlotly({
output$distPlot <- renderPlot({
  inFile <- input$file1
  if (is.null(inFile) || input$package != "ggplot2"){
    return(NULL)
  }
  # df <- data.frame(xx = (1:length(r.x()$x$pvalues))[r.x()$x$maf >= input$min.maf], 
  #                  yy = -log10(r.x()$x$pvalues)[r.x()$x$maf >= input$min.maf])
  # if (is.null(input$file3)){
  #   txt <- (1:length(r.x()$x$pvalues))[r.x()$x$maf >= input$min.maf]
  # } else {
  #   txt <- as.character(read.table(input$file3$datapath)[, 1])
  # }
  if (input$package == "ggplot2"){
    plot(r.x()$x, option = "manhattan")
  }
  
  # if (r.plotPackage()$plotPackage == "plotly"){
  #   #plot(r.x()$x, option = "manhattan", plt.pkg = "plotly")
  # plotly::plot_ly(df, x = ~xx, y = ~yy,
  #                 text = ~paste('SNP: ', txt),
  #                 mode = "markers",
  #                 type = "scatter",
  #                 hoverinfo = "text") %>%
  #   plotly::layout(xaxis = list(title = "Index", showgrid = F),
  #                  yaxis = list(title = "-log10(p-values)"))
  # } else {
  #   plot(r.x()$x, option = "manhattan")
  # }
})

output$distPlotly <- plotly::renderPlotly({
  inFile <- input$file1
  if (is.null(inFile) || (input$package != "plotly")){
    return(NULL)
  }
  if (input$package == "plotly"){
    plot(r.x()$x, option = "manhattan", plt.pkg = "plotly")
  }
})

