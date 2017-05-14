library(shiny)
library(plotly)
library(rjson)
suppressPackageStartupMessages(library(shinyjs))
library(shinyAce)
library(shinyBS)
library(shinydashboard)

header <- dashboardHeader(
  title = "pcadapt"
)

body <- dashboardBody(
  useShinyjs(),
  theme = "bootstrap.css",
  tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),

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
      numericInput("i", label = "i", value = 1, min = 1),
      numericInput("j", label = "j", value = 2, min = 1)
    ),

    mainPanel(
      tabBox(
        tabPanel("PCA", plotlyOutput("pcaPlot")),
        tabPanel("Manhattan Plot", plotlyOutput("distPlot")),
        tabPanel("p-values histogram", plotlyOutput("histPlot"))
      )
    )
  )
)

dashboardPage(
  header,
  dashboardSidebar(disable = TRUE),
  body,
  skin = "blue"
)