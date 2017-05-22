tabPanel(title = strong("PCA"),
         fluidRow(
           column(6, numericInput("i", label = "x-axis", value = 1, min = 1)),
           column(6, numericInput("j", label = "y-axis", value = 2, min = 1))
         ),
         fixedRow(
           column(10, 
                  plotly::plotlyOutput("pcaPlot")
           )
         ), value = 2
)