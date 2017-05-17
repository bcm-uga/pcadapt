tabPanel(title = strong("PCA"),
         numericInput("i", label = "i", value = 1, min = 1),
         numericInput("j", label = "j", value = 2, min = 1),
         fixedRow(
           column(10, 
                  plotly::plotlyOutput("pcaPlot")
           )
         )
)