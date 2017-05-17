tabPanel(
  strong("Manhattan Plot"), 
  fluidRow(column(6, 
                  fileInput( 
                    "file_chr",
                    div("Add chromosome file",
                        div(tags$a(href = "https://drive.google.com/uc?export=download&id=0B9o4VIJJSodfbk9jVXJYS1p4MGM", 
                                   h6("download example file"),
                                   target = "_blank")
                        )
                    ),
                    multiple = FALSE,
                    accept = c(
                      'text/csv',
                      'text/comma-separated-values',
                      '.csv',
                      '.txt'
                    )
                  )
  ),
  column(6, 
         selectInput("package", label = "package", choices = c("ggplot2", "plotly"))
  )
  
  ),
  conditionalPanel("input.package == 'ggplot2'",
                   mainPanel(
                     plotOutput("distPlot"), width = 12
                   )),
  conditionalPanel("input.package == 'plotly'",
                   mainPanel(
                     plotlyOutput("distPlotly"), width = 12
                   ))
  
  
  
  
)

