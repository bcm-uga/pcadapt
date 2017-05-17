library(shiny)
library(plotly)
library(rjson)
suppressPackageStartupMessages(library(shinyjs))
library(shinyAce)
library(shinyBS)
library(shinydashboard)
library(shinythemes)

source(file.path("ui", "helpers.R"))

ui = fluidPage(
  tagList(
    #shinythemes::themeSelector(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
      tags$script(type="text/javascript", src = "busy.js")
    )
    
  ), 
  div(class = "busy",
      p("Calculation in progress..."),
      img(src = "https://media.giphy.com/media/DIQqlqU0SeKXe/giphy.gif")
  ),
  navbarPage(
    theme = shinythemes::shinytheme("flatly"),
    title = "pcadapt",
    tabPanel("Local Adaptation",
             sidebarPanel(
               fileInput(
                 "file1",
                 div("Choose pcadapt file",
                     #helpPopup('Not working atm'), br(),
                     div(tags$a(href = "https://drive.google.com/uc?export=download&id=0B9o4VIJJSodfOTNOelRKNm9MQ2M", 
                                h6("download example file"),
                                target = "_blank")
                     )
                 ),
                 multiple = FALSE,
                 accept = c(
                   'text/csv',
                   'text/comma-separated-values',
                   '.csv',
                   ".pcadapt"
                 )
               ),
               fileInput(
                 "file2",
                 div("Choose population file",
                     #helpPopup('Not working atm'), br(),
                     div(tags$a(href = "https://drive.google.com/uc?export=download&id=0B9o4VIJJSodfaUwwM3JZdVZTRWs", 
                                h6("download example file"), 
                                target = "_blank")
                     )
                 ),
                 multiple = FALSE,
                 accept = c(
                   'text/csv',
                   'text/comma-separated-values',
                   '.csv',
                   '.pop',
                   '.txt',
                   '.fam')
               ),
               fileInput(
                 "file3",
                 div("Choose SNP file",
                     #helpPopup('Not working atm'), br(),
                     div(tags$a(href = "https://drive.google.com/uc?export=download&id=0B9o4VIJJSodfbmJJQjRPcGRBRUU", 
                                h6("download example file"), 
                                target = "_blank")
                     )
                 ),
                 multiple = FALSE,
                 accept = c(
                   'text/csv',
                   'text/comma-separated-values',
                   '.csv',
                   '.txt',
                   '.snp')
               ),
               numericInput("K", label = "K", value = 2, min = 1),
               numericInput("ploidy", label = "ploidy", value = 2, min = 1, max = 2),
               numericInput("min.maf", label = "min.maf", value = 0.05, min = 0.0, max = 0.45, step = 0.05),
               width = 4
             ),
             mainPanel(
               tabsetPanel(
                 source(file.path("ui", "tab-screeplot.R"), local = TRUE)$value,
                 source(file.path("ui", "tab-pca.R"), local = TRUE)$value,
                 source(file.path("ui", "tab-manhattan.R"), local = TRUE)$value,
                 source(file.path("ui", "tab-histogram.R"), local = TRUE)$value,
                 source(file.path("ui", "tab-outlier.R"), local = TRUE)$value,
                 source(file.path("ui", "tab-rcommand.R"), local = TRUE)$value
               ), width = 8
             )
    ),
    tabPanel("Introgression")
  )
)