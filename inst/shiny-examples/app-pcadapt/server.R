library(pcadapt)
library(plotly)
library(ggplot2)
library(shiny)
library(rjson)

options(shiny.maxRequestSize = 2000 * 1024^2) 

shiny::shinyServer(function(input, output) {
  
  r.x <- reactive({
    inFile <- input$file1
    file <- inFile$datapath
    K <- input$K  
    ploidy <- input$ploidy
    min.maf <- input$min.maf

    if (is.null(inFile)){
      return(NULL)
    }
    obj.pcadapt <- pcadapt::pcadapt(file, K = K, ploidy = ploidy, min.maf = min.maf) 
    list(x = obj.pcadapt)
  })
  
  r.ij <- reactive({
    list(i = input$i, j = input$j)
  })
  
  r.pop <- reactive({
    list(inPop = input$file2)
  })
  
  r.ID <- reactive({
    list(inSNP = input$file3)
  })
  
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
      layout(xaxis = list(title = "K", showgrid = F, autotick = FALSE,
                          dtick = 1),      
             yaxis = list(title = "Proportion of variance explained"))
  })
  
  output$distPlot <- plotly::renderPlotly({
    inFile <- input$file1
    if (is.null(inFile)){
      return(NULL)
    }
    df <- data.frame(xx = (1:length(r.x()$x$pvalues))[r.x()$x$maf >= input$min.maf], 
                     yy = -log10(r.x()$x$pvalues)[r.x()$x$maf >= input$min.maf])
    if (is.null(input$file3)){
      txt <- (1:length(r.x()$x$pvalues))[r.x()$x$maf >= input$min.maf]
    } else {
      txt <- as.character(read.table(input$file3$datapath)[, 1])
    }
    plotly::plot_ly(df, x = ~xx, y = ~yy,
                    text = ~paste('SNP: ', txt),
                    mode = "markers", 
                    type = "scatter", 
                    hoverinfo = "text") %>%
      layout(xaxis = list(title = "Index", showgrid = F),      
             yaxis = list(title = "-log10(p-values)"))

  })
  
  output$histPlot <- plotly::renderPlotly({
    inFile <- input$file1
    file <- inFile$datapath
    
    if (is.null(inFile)){
      return(NULL)
    }

    df <- data.frame(xx = r.x()$x$pvalues)
    plotly::plot_ly(df, x = ~xx, type = "histogram") %>%
      layout(xaxis = list(title = "p-values"),
             yaxis = list(title = "Counts"))
  })
  
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
        layout(xaxis = list(title = paste0("PC", r.ij()$i), showgrid = F),
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
        layout(xaxis = list(title = paste0("PC", r.ij()$i), showgrid = F),
               yaxis = list(title = paste0("PC", r.ij()$j)))
    }
  })
  
  output$outlierTable <- DT::renderDataTable({
    inFile <- input$file1
    file <- inFile$datapath
    if (is.null(inFile)){
      return(NULL)
    }
    
    nSNP <- length(r.x()$x$pvalues)
    sorted.snp <- sort(r.x()$x$pvalues, index.return = TRUE)
    pc <- get.pc(r.x()$x, 1:nSNP)
    
    if (is.null(r.ID()$inSNP)){
      ID <- rep(NA, nSNP)
      df <- data.frame(Rank = 1:nSNP, ID = ID, Index = sorted.snp$ix, 
                       pvalue  = sorted.snp$x, PC = as.integer(pc[, 2]))
    } else {
      ID <- read.table(r.ID()$inSNP$datapath, header = FALSE)[,1]
      df <- data.frame(Rank = 1:nSNP, ID = ID, Index = sorted.snp$ix, 
                       pvalue  = sorted.snp$x, PC = as.integer(pc[, 2]))
    }
    

    dt <- DT::datatable(df, class = 'cell-border stripe', rownames = FALSE)
  })
  
  output$rcommand <- renderText({
    inFile <- input$file1
    if (is.null(inFile)){
      return(NULL)
    }
    paste0("object.pcadapt <- pcadapt(x = ", input$file1$name,
           ", K = ", input$K, ", ploidy = ", input$ploidy,
           ", min.maf = ", input$min.maf, ")")
  })
  
})
