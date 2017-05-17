tabPanel(
  strong("Manhattan Plot"), 
  selectInput("package", label = "package", choices = c("ggplot2", "plotly")),
  # radioButtons("radio", label = h3("pkg"),
  #              choices = list("ggplot2" = 1, 
  #                             "plotly" = 2, 
  #                             selected = 1)
  # ),
  conditionalPanel("input.package == 'ggplot2'",
                   mainPanel(
                     plotOutput("distPlot"), width = 12
                   )),
  conditionalPanel("input.package == 'plotly'",
                   mainPanel(
                     plotlyOutput("distPlotly"), width = 12
                   ))
  
  
  
  
)

