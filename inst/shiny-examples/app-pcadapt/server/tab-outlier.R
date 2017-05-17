output$outlierTable <- DT::renderDataTable({
  inFile <- input$file1
  file <- inFile$datapath
  if (is.null(inFile)){
    return(NULL)
  }
  
  sorted.snp <- sort(r.x()$x$pvalues, index.return = TRUE)
  nna.idx <- which(!is.na(r.x()$x$pvalues))
  nSNP <- length(r.x()$x$pvalues[nna.idx]) 
  pc <- get.pc(r.x()$x, nna.idx[sorted.snp$ix])
  
  if (is.null(r.ID()$inSNP)){
    ID <- rep(NA, nSNP)
    df <- data.frame(Rank = 1:nSNP, ID = ID, Index = sorted.snp$ix, 
                     pvalue  = sorted.snp$x, PC = pc[, 2])
  } else {
    ID <- read.table(r.ID()$inSNP$datapath, header = FALSE)[,1]
    ID <- ID[nna.idx]
    df <- data.frame(Rank = 1:nSNP, ID = ID[sorted.snp$ix], Index = sorted.snp$ix, 
                     pvalue  = sorted.snp$x, PC = pc[, 2])
  }
  
  dt <- DT::datatable(df, class = 'cell-border stripe', rownames = FALSE)
})