library(shiny)
library(plotly)

shiny::shinyUI(fluidPage(
  
  # Application title
  titlePanel("pcadapt"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose pcadapt file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".pcadapt")
      ),
      
      fileInput("file2", "Choose population file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".pop",
                  ".txt")
      ),
      
      selectInput("opt", label = "Option", 
                  choices = c("PCA", "Manhattan")),
      
      numericInput("K", label = "K", value = 2),
      numericInput("i", label = "i", value = 1),
      numericInput("j", label = "j", value = 2)

    ),
    
    mainPanel(
       plotly::plotlyOutput("distPlot")
    )
  )
))
