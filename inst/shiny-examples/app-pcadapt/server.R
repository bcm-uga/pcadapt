library(pcadapt)
library(plotly)
library(ggplot2)
library(shiny)
library(rjson)

options(shiny.maxRequestSize = 2000 * 1024^2) 
options(warn = -1)

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
  
  r.chr <- reactive({
    list(inCHR = input$file_chr)
  })
  
  r.plotPackage <- reactive({
    list(plotPackage = input$plotPackage)
  })
  
  r.file_intrg_pop <- reactive({
    if (is.null(input$file_intrg_pop)){
      return(NULL)
    }
    list(path = input$file_intrg_pop$datapath)
  })
  
  output$file_intrg_pop <- reactive({
    return(!is.null(r.file_intrg_pop()))
  })
  
  outputOptions(output, 'file_intrg_pop', suspendWhenHidden = FALSE)
  
  source(file.path("server", "tab-screeplot.R"), local = TRUE)$value
  source(file.path("server", "tab-pca.R"), local = TRUE)$value
  source(file.path("server", "tab-manhattan.R"), local = TRUE)$value  
  source(file.path("server", "tab-histogram.R"), local = TRUE)$value
  source(file.path("server", "tab-outlier.R"), local = TRUE)$value
  source(file.path("server", "tab-rcommand.R"), local = TRUE)$value
  
})