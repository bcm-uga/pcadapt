tabPanel(
  strong("Manhattan Plot"), 
  fluidRow(
    column(width = 6,
           selectInput("package", label = "package", choices = c("ggplot2", "plotly")),
           conditionalPanel(condition = "input.package == 'ggplot2'", {
             div(actionButton("plotmanhattan_ggplot", "Plot", icon = icon("bar-chart")))  
           }),
           conditionalPanel(condition = "input.package == 'plotly'", {
             div(actionButton("plotmanhattan_plotly", "Plot", icon = icon("bar-chart")))  
           })
    )
  ),
  tags$style(type='text/css', "#plotmanhattan_ggplot { width:100%; margin-top: 15px;}"),
  tags$style(type='text/css', "#plotmanhattan_plotly { width:100%; margin-top: 15px;}"),
  br(),
  conditionalPanel("input.package == 'ggplot2'", 
                   fluidPage(
                     mainPanel(
                       plotOutput("distPlot"), width = 12
                     )
                   )
  ),
  conditionalPanel("input.package == 'plotly'",
                   fluidPage(
                     mainPanel(
                       plotly::plotlyOutput("distPlotly"), width = 12
                     )
                   )
  ),
  value = 3
)

