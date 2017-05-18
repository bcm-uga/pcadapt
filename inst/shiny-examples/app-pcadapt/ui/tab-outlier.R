tabPanel(
  strong("Ranked SNPs"), 
  DT::dataTableOutput("outlierTable"),
  value = 5
)