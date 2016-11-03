#' pcadapt visualization tool
#'
#' \code{plot.pcadapt} is a method designed for objects of class \code{pcadapt}.
#' It provides a plotting utile for quick visualization of \code{pcadapt} objects.
#' Different options are currently available : \code{"screeplot"}, \code{"scores"}, \code{"stat.distribution"},
#' \code{"manhattan"} and \code{"qqplot"}.
#' \code{"screeplot"} shows the decay of the genotype matrix singular values and provides
#' a figure to help with the choice of \code{K}.
#' \code{"scores"} plots the projection of the individuals onto the first two principal components.
#' \code{"stat.distribution"} displays the histogram of the selected test statistics, as well as
#' the estimated distribution for the neutral SNPs.
#' \code{"manhattan"} draws the Manhattan plot of the p-values associated with the statistic of interest.
#' \code{"qqplot"} draws a Q-Q plot of the p-values associated with the statistic of interest.
#'
#' @param x an object of class "pcadapt" generated with \code{pcadapt}.
#' @param ... \dots
#' @param option a character string specifying the figures to be displayed. If \code{NULL} (the default), all three plots are printed.
#' @param K an integer specifying the principal component of interest. \code{K} has to be specified only when using the \code{loadings} option.
#' @param i an integer indicating onto which principal component the individuals are projected when the "scores" option is chosen.
#' Default value is set to \code{1}.
#' @param j an integer indicating onto which principal component the individuals are projected when the "scores" option is chosen.
#' Default value is set to \code{2}.
#' @param pop a list of integers or strings specifying which subpopulation the individuals belong to.
#' @param threshold for the \code{"qqplot"} option, it displays an additional bar which shows the \code{threshold} percent of SNPs with smallest p-values
#' and separates them from SNPs with higher p-values.
#'
#' @examples
#' ## see ?pcadapt for examples
#'
#' @method plot pcadapt
#'
#' @export
plot.pcadapt = function(x,...,option="manhattan",K=NULL,i=1,j=2,pop,threshold=NULL){
  if (!(option %in% c("screeplot","scores","manhattan","qqplot","stat.distribution"))){
    warning(paste("Plotting option",option,"not valid, options currently available are: screeplot, scores, manhattan, qqplot, stat.distribution."))
  } else {
    if (option == "screeplot"){
      scree.plotting(x,K)
    } else if (option == "scores"){
      if (attr(x,"data.type")!="pool"){
        if (missing(pop)){
          score.plotting(x,i,j)
        } else {
          score.plotting(x,i,j,pop)
        }
      } else {
        score.plotting(x,i,j,pop=1:dim(x$scores)[1])
      }
    } else if (option == "stat.distribution"){
      if ((attr(x,"method") %in% c("mahalanobis","communality")) == FALSE){
        if (is.null(K)){
          warning("K has to be specified.")
        } else {
          neutral.plotting(x,K)
        }
      } else {
        neutral.plotting(x,1)
      }
    } else if (option == "manhattan"){
      if ((attr(x,"method") %in% c("mahalanobis","communality")) == FALSE){
        if (is.null(K)){
          warning("K has to be specified.")
        } else {
          manhattan.plotting(x,K)
        }
      } else {
        manhattan.plotting(x,K=1)
      }
    } else if (option == "qqplot"){
      if ((attr(x,"method") %in% c("mahalanobis","communality")) == FALSE){
        if (is.null(K)){
          warning("K has to be specified")
        } else{
          pvalqq.plotting(x,K,threshold=threshold)
        }
      } else {
        pvalqq.plotting(x,K=1,threshold=threshold)
      }
    }
  }
}

#' Principal Components Analysis Scores Plot
#'
#' \code{"score.plotting"} plots the projection of the individuals onto the first two principal components.
#'
#' @param x an output from \code{pcadapt} containing the scores.
#' @param i an integer indicating onto which principal component the individuals are projected when the "scores" option is chosen.
#' Default value is set to \code{1}.
#' @param j an integer indicating onto which principal component the individuals are projected when the "scores" option is chosen.
#' Default value is set to \code{2}.
#' @param pop a list of integers or strings specifying which subpopulation the individuals belong to.
#'
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom ggplot2 ggplot ggtitle labs geom_point guides aes aes_string geom_point scale_color_hue
#'
#' @export
#'
score.plotting = function(x,i=1,j=2,pop){
  
  if (attr(x,"K")==1){
    warning("K=1, option not available since two principal components have to be computed at least.")
  } else {
    
    if (i == j){
      stop("j has to be different from i.")
    }
    
    if (i>attr(x,"K")){
      stop(paste0("i can't exceed ",attr(x,"K"),"."))
    }
    
    if (j>attr(x,"K")){
      stop(paste0("j can't exceed ",attr(x,"K"),"."))
    }
    if (missing(pop)){
      ggdf <- as.data.frame(cbind(x$scores[,i],x$scores[,j]))  
      colnames(ggdf) <- c("PC_i","PC_j")
      res.plot <- ggplot2::ggplot(ggdf,aes_string("PC_i","PC_j")) +
        geom_point() 
    } else {
      pop.to.int <- get.score.color(pop)
      popnames <- get.pop.names(pop)
      ggdf <- as.data.frame(cbind(x$scores[,i],x$scores[,j],pop.to.int)) 
      colnames(ggdf) <- c("PC_i","PC_j","Pop")
      res.plot <- ggplot2::ggplot(ggdf,aes_string("PC_i","PC_j")) +
        ggplot2::geom_point(aes(colour=factor(ggdf$Pop))) +
        ggplot2::scale_color_hue(name=" ",labels=popnames)
    }
    res.plot <- res.plot + ggplot2::ggtitle(paste0("Projection onto PC",i," and PC",j)) +
      ggplot2::labs(x=paste0("PC",i),y=paste0("PC",j))
    print(res.plot)
  }
}

#' Manhattan Plot
#'
#' \code{manhattanPlot} displays a Manhattan plot which represents the p-values for each SNP for a particular test statistic.
#'
#' @param x an object of class "pcadapt" generated with \code{pcadapt} containing the p-values of interest.
#'
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom ggplot2 qplot guides ggtitle
#'
#' @export
#'
manhattan.plotting = function(x,K){
  if (K > attr(x,"K")){
    stop(paste0("K can't exceed ",attr(x,"K")),".")
  }
  if (attr(x,"method")=="componentwise"){
    pval.K <- x$pvalues[!is.na(x$pvalues[,K]),K]
  } else {
    pval.K <- x$pvalues[!is.na(x$pvalues)]
  }
  p0 <- ggplot2::qplot(1:length(pval.K),-log10(pval.K),col="red",xlab=paste0("SNP (with mAF>",attr(x,"min.maf"),")"),ylab="-log10(p-values)") +
    ggplot2::guides(colour=FALSE) + 
    ggplot2::ggtitle("Manhattan Plot")
  print(p0)
}

#' Principal Components Analysis Scree Plot
#'
#' \code{scree.plotting} plots the scee plot associated with the principal components analysis performed on the dataset.
#' NB : \code{pcadapt} has to be run on the dataset in order to get an output readable by \code{plot.screePlot}
#'
#' @param x an output from \code{pcadapt} containing the singular values.
#' @param K an integer specifying the number of components to take into account in the scree plot.
#'
#' @examples
#' ## see ?fastpcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom ggplot2 qplot geom_line guides ggtitle
#'
#' @export
#'
scree.plotting = function(x,K){
  if (is.null(K)){m <- attr(x,"K")}
  else {m <- K}
  if (m<2){warning("K = 1, the scree plot is thus composed of a unique point.")}
  nSNP <- length(x$maf)
  p0 <- ggplot2::qplot(x=1:m,y=(x$singular.values[1:m])^2/nSNP,col="red",xlab="PC",ylab="Proportion of explained variance") + 
    ggplot2::geom_line() + ggplot2::guides(colour=FALSE) +
    ggplot2::ggtitle(paste("Scree Plot - K =",m))
  print(p0)
}

#' p-values Q-Q Plot
#'
#' \code{pvalqq.plotting} plots a Q-Q plot of the p-values computed.
#'
#' @param x an output from \code{outlier} containing the p-values of interest.
#' @param K an integer specifying which principal component to display when \code{method="componentwise"}.
#' @param threshold a real number between \code{0} et \code{1}.
#'
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom ggplot2 qplot geom_vline guides geom_abline ggtitle aes
#'
#' @export
#'
pvalqq.plotting = function(x,K,threshold){
  if (attr(x,"method")=="componentwise"){
    sorted.pval <- sort(x$pvalues[x$maf>=attr(x,"min.maf"),K])
  } else {
    sorted.pval <- sort(x$pvalues[x$maf>=attr(x,"min.maf")])
  }
  p <- length(sorted.pval)
  expected.p <- 1:p/p
  p0 <- ggplot2::qplot(-log10(expected.p),-log10(sorted.pval),col="red",xlab="Expected -log10(p-values)",ylab="Observed -log10(p-values)") + 
    ggplot2::geom_abline()
  if (!missing(threshold)){
    q <- floor(threshold*p)
    pval.thresh <- expected.p[q]
    p0 <- p0 + ggplot2::geom_vline(aes(xintercept = -log10(pval.thresh)),colour="blue")
  }
  p0 <- p0 + ggplot2::ggtitle("Q-Q plot") + ggplot2::guides(colour=FALSE)
  print(p0)
}

#' Neutral Distribution Estimation
#'
#' \code{neutral.plotting} plots the histogram of the chi-squared statistics, as well as the estimated null distribution.
#'
#' @param x an output from \code{outlier} containing the chi-squared statistics.
#' @param K an integer indicating which principal component the histogram will be associated with.
#'
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom stats dchisq
#' @importFrom ggplot2 ggplot geom_histogram geom_line theme aes_string element_rect ggtitle
#'
#' @export
neutral.plotting = function(x,K){
  idxmaf <- x$maf>=attr(x,"min.maf")
  if (attr(x,"method")=="componentwise"){
    df <- 1
    z <- x$chi2.stat[idxmaf,df]
  } else if (attr(x,"method") != "componentwise" && attr(x,"data.type") != "pool"){
    df <- attr(x,"K")
    z <- x$chi2.stat[idxmaf]
  }
  min.z <- floor(min(z[which(!is.na(z))]))
  max.z <- floor(max(z[which(!is.na(z))])+1)
  if (max.z > 1e5){
    stop("Can't display the histogram as the values are too high.")
  }
  xx <- seq(min.z,max.z,length=length(z))
  ggdf <- as.data.frame(cbind(xx,dchisq(xx,df=df),z))
  colnames(ggdf) <- c("abs","ord","chi2")
  p0 <- ggplot() + 
    geom_histogram(data=ggdf,aes_string(x="chi2",y="..density.."),binwidth=0.5,fill="#B0E2FF",alpha=0.6,colour="black") + 
    geom_line(data=ggdf,aes_string(x="abs",y="ord"),col="#4F94CD",size=1) + 
    ggplot2::ggtitle("Statistics distribution")
  print(p0)
}

