output$rcommand <- renderText({
  inFile <- input$file1
  if (is.null(inFile)){
    return(NULL)
  }
  paste0("object.pcadapt <- pcadapt(x = '", input$file1$name,
         "', K = ", input$K, ", ploidy = ", input$ploidy,
         ", min.maf = ", input$min.maf, ")")
})