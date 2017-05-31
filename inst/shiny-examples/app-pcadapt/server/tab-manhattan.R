makePlot <- eventReactive(input$plotmanhattan_ggplot, {
  inFile <- input$file1
  if (is.null(inFile) || input$package != "ggplot2"){
    return(NULL)
  } else if (input$package == "ggplot2" && is.null(input$file_chr)){
    plot(r.x()$x, option = "manhattan")
  } else if (input$package == "ggplot2" && !is.null(input$file_chr)){
    chr.info <- read.table(r.chr()$inCHR$datapath)[, 1]
    plot(r.x()$x, option = "manhattan", chr.info = chr.info)
  }
})

makePlotly <- eventReactive(input$plotmanhattan_plotly, {
  inFile <- input$file1
  if (is.null(inFile) || (input$package != "plotly")){
    return(NULL)
  } else if (input$package == "plotly" && is.null(input$file_chr)){
    
    plot(r.x()$x, option = "manhattan", plt.pkg = "plotly")
  } else if (input$package == "plotly" && !is.null(input$file_chr)){
    chr.info <- read.table(r.chr()$inCHR$datapath)[, 1]
    plot(r.x()$x, option = "manhattan", plt.pkg = "plotly", chr.info = chr.info)
  }
})

output$distPlot <- renderPlot({
  makePlot()
})

output$distPlotly <- plotly::renderPlotly({
  makePlotly()
})

