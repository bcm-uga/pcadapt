library(shiny)
library(plotly)
library(rjson)
suppressPackageStartupMessages(library(shinyjs))
library(shinyAce)
library(shinyBS)
library(shinydashboard)
library(shinythemes)

#source(file.path("ui", "helpers.R"))

ui = fluidPage(
  tagList(
    #shinythemes::themeSelector(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
      tags$script(type="text/javascript", src = "busy.js")
    )
    
  ), 
  div(class = "busy",
      p("Calculation in progress.."),
      img(src = "https://media.giphy.com/media/DIQqlqU0SeKXe/giphy.gif")
  ),
  navbarPage(
    theme = shinythemes::shinytheme("flatly"),
    title = "pcadapt",
    tabPanel("Local Adaptation",
             sidebarPanel(
               fileInput("file1", "Choose pcadapt file",
                         
                         multiple = FALSE,
                         accept = c(
                           "text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv",
                           ".pcadapt")
               ),
               downloadLink("sampleDataFile", "Example data file"),
               fileInput("file2", "Choose population file",
                         accept = c(
                           "text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv",
                           ".pop",
                           ".txt",
                           ".fam")
               ),
               
               fileInput("file3", "Choose SNP file",
                         accept = c(
                           "text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv",
                           ".txt",
                           ".snp")
               ),
               
               numericInput("K", label = "K", value = 2, min = 1),
               numericInput("ploidy", label = "ploidy", value = 2, min = 1, max = 2),
               numericInput("min.maf", label = "min.maf", value = 0.05, min = 0.0, max = 0.45),
               width = 4
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel(strong("Screeplot"), plotlyOutput("screePlot")),
                 tabPanel(title = strong("PCA"),
                          br(),
                          numericInput("i", label = "i", value = 1, min = 1),
                          numericInput("j", label = "j", value = 2, min = 1),
                          fixedRow(column(10, plotly::plotlyOutput("pcaPlot")))),
                 tabPanel(strong("Manhattan Plot"), plotly::plotlyOutput("distPlot")),
                 tabPanel(strong("p-values histogram"), plotly::plotlyOutput("histPlot")),
                 tabPanel(strong("Ranked SNPs"), DT::dataTableOutput("outlierTable")),
                 tabPanel(strong("R command"), verbatimTextOutput("rcommand"))
               ), width = 8
             )
    ),
    tabPanel("Introgression")
  )
)



# header <- dashboardHeader(
#   title = "pcadapt"
# )
# 
# body <- dashboardBody(
#   useShinyjs(),
#   theme = "bootstrap.css",
#   tagList(
#     tags$head(
#       tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
#       tags$script(type="text/javascript", src = "busy.js")
#     )
#   ),
#   div(class = "busy",  
#       p("Calculation in progress.."), 
#       img(src = "https://media.giphy.com/media/DIQqlqU0SeKXe/giphy.gif")
#   ),
#   sidebarLayout(
#     sidebarPanel(
#       fileInput("file1", "Choose pcadapt file",
#                 accept = c(
#                   "text/csv",
#                   "text/comma-separated-values,text/plain",
#                   ".csv",
#                   ".pcadapt")
#       ),
# 
#       fileInput("file2", "Choose population file",
#                 accept = c(
#                   "text/csv",
#                   "text/comma-separated-values,text/plain",
#                   ".csv",
#                   ".pop",
#                   ".txt",
#                   ".fam")
#       ),
# 
#       fileInput("file3", "Choose SNP file",
#                 accept = c(
#                   "text/csv",
#                   "text/comma-separated-values,text/plain",
#                   ".csv",
#                   ".txt",
#                   ".snp")
#       ),
# 
#       numericInput("K", label = "K", value = 2, min = 1),
#       numericInput("ploidy", label = "ploidy", value = 2, min = 1, max = 2),
#       numericInput("min.maf", label = "min.maf", value = 0.05, min = 0.0, max = 0.45),
#       width = 3
#     ),
# 
#     mainPanel(
#       tabBox(
#         tabPanel(strong("Screeplot"), plotlyOutput("screePlot")),
#         tabPanel(title = strong("PCA"), 
#                  br(),
#                  numericInput("i", label = "i", value = 1, min = 1),
#                  numericInput("j", label = "j", value = 2, min = 1),
#                  fixedRow(column(10, plotly::plotlyOutput("pcaPlot")))),
#         tabPanel(strong("Manhattan Plot"), plotly::plotlyOutput("distPlot")),
#         tabPanel(strong("p-values histogram"), plotly::plotlyOutput("histPlot")),
#         tabPanel(strong("Ranked SNPs"), DT::dataTableOutput("outlierTable")),
#         tabPanel(strong("R command"), verbatimTextOutput("rcommand"))
#       ), width = 9
#     )
#   )
# )
# 
# dashboardPage(
#   dashboardSidebar(disable = TRUE),
#   body,
#   skin = "blue"
# )