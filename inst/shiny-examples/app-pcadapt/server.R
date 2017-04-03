library(pcadapt)
library(plotly)
library(ggplot2)
library(shiny)
library(rjson)

shiny::shinyServer(function(input, output) {
  
  output$distPlot <- plotly::renderPlotly({
    inFile <- input$file1
    inPop <- input$file2
    file <- inFile$datapath
    K <- input$K  
    ploidy <- input$ploidy
    min.maf <- input$min.maf
    i <- input$i
    j <- input$j
    
    if (is.null(inFile)){
      return(NULL)
    }
    
    x <- pcadapt::pcadapt(file, K = K, ploidy = ploidy, min.maf = min.maf)
    
    if (input$opt == "Manhattan"){
      df <- data.frame(xx = 1:length(x$pvalues), yy = -log10(x$pvalues))
      plotly::plot_ly(df, x = ~xx, y = ~yy,
                      text = ~paste('SNP: ', 1:length(x$pvalues)),
                      mode = "markers", type = "scatter", hoverinfo = "text") %>%
        layout(p,                       
               title = "Manhattan Plot", 
               xaxis = list(   
                 title = "Index",
                 showgrid = F),      
               yaxis = list(          
                 title = "-log10(p-values)")   
        )
    } else if (input$opt == "Histogram"){
      df <- data.frame(xx = x$pvalues)
      plotly::plot_ly(df, x = ~xx, type = "histogram") %>%
        layout(p,                       
               title = "Histogram", 
               xaxis = list(   
                 title = "p-values",
                 showgrid = F),
               yaxis = list(          
                 title = "Counts")
        )
    } else if (input$opt == "PCA"){ 
      if (is.null(inPop)){
        df.scores <- data.frame(xx = x$scores[, i], yy = x$scores[, j])
        plotly::plot_ly(df.scores, x = ~xx, y = ~yy, 
                        mode = "markers", type = "scatter", hoverinfo = "text") %>%
          layout(p,                       
                 title = "PCA", 
                 xaxis = list(   
                   title = paste0("PC", i),
                   showgrid = F),      
                 yaxis = list(          
                   title = paste0("PC", j))      
          )
        
      } else {
        pop <- as.character(read.table(inPop$datapath)[, 1])
        df.scores <- data.frame(xx = x$scores[, i], yy = x$scores[, j],
                                pop = pop, ind = 1:nrow(x$scores))
        plotly::plot_ly(df.scores, x = ~xx, y = ~yy, 
                        color = ~pop, text = ~paste('Ind: ', 1:nrow(x$scores)),
                        mode = "markers", type = "scatter", hoverinfo = "text") %>%
          layout(p,                       
                 title = "PCA", 
                 xaxis = list(   
                   title = paste0("PC", i),
                   showgrid = F),      
                 yaxis = list(          
                   title = paste0("PC", j))      
          )
      }
    }
  })
  
})
