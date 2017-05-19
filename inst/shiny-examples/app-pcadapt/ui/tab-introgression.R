tabPanel("Introgression",
         column(width = 4,
                fluidRow(
                  box(width = NULL, status = "primary",
                      fileInput(
                        "file_intrg_geno",
                        div("Choose pcadapt file",
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
                        "file_intrg_pop",
                        div("Choose population file",
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
                  )
                ), 
                fluidRow(
                  box(width = NULL, status = "primary",
                      numericInput("K_intrg", label = "K", value = 2, min = 1),
                      numericInput("ploidy_intrg", label = "ploidy", value = 2, min = 1, max = 2),
                      numericInput("min.maf_intrg", label = "min.maf", value = 0.05, min = 0.0, max = 0.45, step = 0.05)
                  )
                ),
                fluidRow(
                  conditionalPanel(condition = "output.file_intrg_pop == true",
                                   box(width = NULL, status = "primary",
                                       selectInput("ancstrl1", "Ancestral 1", choices = c("2", "Tricho"), 
                                                   selected = NULL, 
                                                   multiple = FALSE,
                                                   selectize = TRUE,
                                                   width = NULL,
                                                   size = NULL),
                                       selectInput("ancstrl2", "Ancestral 2", choices = c("2", "Tricho"), 
                                                   selected = NULL, 
                                                   multiple = FALSE,
                                                   selectize = TRUE,
                                                   width = NULL,
                                                   size = NULL)
                                       
                                   )
                  )
                )
         ),
         mainPanel(width = 8
         )
)