output$outlierTable <- DT::renderDataTable({
  inFile <- input$file1
  file <- inFile$datapath
  if (is.null(inFile)){
    return(NULL)
  }
  
  nSNP <- length(r.x()$x$pvalues)
  sorted.snp <- sort(r.x()$x$pvalues, index.return = TRUE)
  pc <- get.pc(r.x()$x, 1:nSNP)
  
  if (is.null(r.ID()$inSNP)){
    ID <- rep(NA, nSNP)
    df <- data.frame(Rank = 1:nSNP, ID = ID, Index = sorted.snp$ix, 
                     pvalue  = sorted.snp$x, PC = as.integer(pc[, 2]))
  } else {
    ID <- read.table(r.ID()$inSNP$datapath, header = FALSE)[,1]
    df <- data.frame(Rank = 1:nSNP, ID = ID, Index = sorted.snp$ix, 
                     pvalue  = sorted.snp$x, PC = as.integer(pc[, 2]))
  }
  
  dt <- DT::datatable(df, class = 'cell-border stripe', rownames = FALSE)
})