tabPanel(
  strong("Manhattan Plot"), 
  selectInput("package", label = "package", choices = c("ggplot2", "plotly")),
  conditionalPanel("input.package == 'ggplot2'",
                   mainPanel(
                     plotOutput("distPlot"), width = 12
                   )),
  conditionalPanel("input.package == 'plotly'",
                   mainPanel(
                     plotlyOutput("distPlotly"), width = 12
                   )),
  value = 3
)

