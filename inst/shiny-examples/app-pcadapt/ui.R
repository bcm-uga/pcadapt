library(shiny)
library(plotly)
library(rjson)
suppressPackageStartupMessages(library(shinyjs))
library(shinyAce)
library(shinyBS)
library(shinydashboard)
library(shinythemes)

source(file.path("ui", "helpers.R"))

header <- dashboardHeader(
  title = "pcadapt"
)

body <- dashboardBody(  
  tagList(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
      tags$script(type="text/javascript", src = "busy.js")
    )
    
  ), 
  div(class = "busy",
      p("Calculation in progress..."),
      img(src = "https://media.giphy.com/media/DIQqlqU0SeKXe/giphy.gif")
  ),
  tabItems(
    tabItem("local-adaptation",
            fluidPage(
              theme = shinythemes::shinytheme("flatly"),
              title = " ",
              source(file.path("ui", "tab-local-adaptation.R"), local = TRUE)$value
            )
    ),
    tabItem("introgression",
            fluidPage(
              theme = shinythemes::shinytheme("flatly"),
              title = " ",
              source(file.path("ui", "tab-introgression.R"), local = TRUE)$value  
            )
    )
  )
)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Local Adaptation", tabName = "local-adaptation", icon = icon("bar-chart-o")),
    menuItem("Introgression", icon = icon("th"), tabName = "introgression", badgeLabel = "dev", badgeColor = "red"), 
    menuItem("Vignette", icon = icon("book"), href = "https://bioshock38.github.io/pcadapt/articles/pcadapt.html"), 
    menuItem("Github", icon = icon("github"), href = "https://github.com/BioShock38/pcadapt"),
    menuItem("CRAN", icon = icon("code"), href = "https://cran.r-project.org/web/packages/pcadapt/index.html") 
  )  
)

dashboardPage(
  header,
  sidebar,
  body
)