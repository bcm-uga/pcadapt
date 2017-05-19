tabPanel("Local Adaptation",
         column(
           fluidRow(
             box(width = NULL, status = "primary",
                 conditionalPanel(condition = "true",
                                  fileInput(
                                    "file1",
                                    div("Choose pcadapt file",
                                        #helpPopup('Not working atm'),
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
                                  numericInput("K", label = "K", value = 2, min = 1),
                                  numericInput("ploidy", label = "ploidy", value = 2, min = 1, max = 2),
                                  numericInput("min.maf", label = "min.maf", value = 0.05, min = 0.0, max = 0.45, step = 0.05)
                 )
             )
           ), 
           fluidRow(
             conditionalPanel(condition = "input.conditionedPanels == 2 || input.conditionedPanels == 3 || input.conditionedPanels == 5",
                              box(width = NULL, status = "primary",
                                  conditionalPanel(condition = "input.conditionedPanels == 2",
                                                   fileInput(
                                                     "file2",
                                                     div("Add population file",
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
                                                   )
                                                   
                                  ),
                                  conditionalPanel(condition = ("input.conditionedPanels == 3"),
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
                                                       '.txt',
                                                       '.chr'
                                                     )
                                                   )
                                  ),
                                  conditionalPanel(condition = ("input.conditionedPanels == 5"),
                                                   fileInput(
                                                     "file3",
                                                     div("Add SNP file",
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
                                                   )
                                  )
 
                              )
             )
           ),
           width = 4
         ),
         mainPanel(
           box(width = NULL, status = "primary",
               tabsetPanel(
                 source(file.path("ui", "tab-screeplot.R"), local = TRUE)$value,
                 source(file.path("ui", "tab-pca.R"), local = TRUE)$value,
                 source(file.path("ui", "tab-manhattan.R"), local = TRUE)$value,
                 source(file.path("ui", "tab-histogram.R"), local = TRUE)$value,
                 source(file.path("ui", "tab-outlier.R"), local = TRUE)$value,
                 source(file.path("ui", "tab-rcommand.R"), local = TRUE)$value,
                 id = "conditionedPanels"
               )  
           ), width = 8
         )
)