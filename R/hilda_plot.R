#' Plot mutation signatures from HiLDA output
#'
#' @param hildaResult a rjags class output by HiLDA
#' @param sigOrder the order of signatures if needed (default: NULL)
#' @param colorList a vector of color for mutational exposures barplots
#' @param ... Other arguments passed on to methods
#' @return a plot object containing all mutational signatures
#'
#' @examples
#' 
#' inputFile <- system.file("extdata/hildaLocal.rdata", package="HiLDA")
#' hildaLocal <- readRDS(inputFile)
#' hildaPlotSignature(hildaLocal)
#'
#' @importFrom cowplot plot_grid
#' @importFrom methods is
#' @export

hildaPlotSignature <- function(hildaResult, sigOrder=NULL, colorList = NULL,
                               ...) {
    if(is(hildaResult) != "rjags") {
        stop("Not an output object from running the HiLDA tests.")
    }
    
    numSig <- dim(hildaResult$BUGSoutput$mean$pStates1)[3]
    
    if (is.null(sigOrder)) {
        sigOrder <- seq_len(numSig)
    }
    
    if(length(sigOrder) != numSig) {
        stop("The order of signatures has more or less samples.")
    }

    
    if (is.null(colorList)) {
        colorList <- hcl(h=seq(15, 375, length=numSig + 1),
                         l=65, c=100)[seq_len(numSig)]
    }

    if(length(colorList) != numSig) {
        stop("The length of the color list has more or less samples.")
    }
    
    plotList <- vector("list", numSig)

    feature <- c(length(hildaResult$BUGSoutput$mean$pStates1[,,1]),
                 rep(ncol(hildaResult$BUGSoutput$mean$pStates2[,,1]),
                     nrow(hildaResult$BUGSoutput$mean$pStates2[,,1])))

    trDir <- "pStates3" %in% names(hildaResult$BUGSoutput$sims.list)

    if(trDir == TRUE){
        feature <- c(feature, 2)
    }

    numBases <- 1 + dim(hildaResult$BUGSoutput$mean$pStates2)[1]

    for (i in sigOrder) {
        tempSig <- matrix(0, nrow=numBases, ncol=6)
        if (trDir == TRUE) {
            tempSig[1, ] <- hildaResult$BUGSoutput$mean$pStates1[, , i]
            tempSig[2:(length(feature) - 1), seq_len(feature[2])] <-
                hildaResult$BUGSoutput$mean$pStates2[, , i]
            tempSig[length(feature), seq_len(2)] <-
                hildaResult$BUGSoutput$mean$pStates3[, , i]
        } else {
            tempSig[1, ] <- hildaResult$BUGSoutput$mean$pStates1[, , i]
            tempSig[2:length(feature), seq_len(feature[2])] <-
                hildaResult$BUGSoutput$mean$pStates2[, , i]
        }

        plotList[[which(sigOrder == i)]] <-
            visPMS(tempSig, numBases=numBases,
                   trDir=trDir, isScale=TRUE, ...) +
            theme(panel.border=
                      element_rect(color=colorList[which(sigOrder ==i)],
                                   fill=NA, size=2, linetype = 1),
                  plot.margin=unit(c(0.01, 0.01, 0.01, 0.01), "npc"))
    }

    cowplot::plot_grid(plotlist=plotList, ncol=1)
}



#' Read the raw mutation data with the mutation feature vector format,
#'   estimate and plot both mutation signatures and their fractions
#'
#' @param inputG a MutationFeatureData S4 class output by the pmsignature.
#' @param hildaResult a rjags class output by HiLD.
#' @param sigOrder the order of signatures if needed (default: NULL).
#' @param refGroup the samples in the reference group (default: NULL).
#' @param sortSampleNum whether to sort plots by number of mutations
#'        (default: TRUE).
#' @param refName the name of reference group (default: Control)
#' @param altName the name of the other group (default: Case)
#' @param charSize the size of the character on the signature plot (default: 3)
#' @return a list of a signature plot and a barplot of mutational exposures
#'
#' @examples
#'
#' load(system.file("extdata/sample.rdata", package="HiLDA"))
#' inputFile <- system.file("extdata/hildaLocal.rdata", package="HiLDA")
#' hildaLocal <- readRDS(inputFile)
#'
#' hildaBarplot(G, hildaLocal, refGroup=1:4)
#'
#' @importFrom cowplot plot_grid
#' @importFrom tidyr gather_
#' @importFrom grid unit.c
#' @importFrom forcats fct_relevel fct_reorder
#' @importFrom stats aggregate
#' @importFrom methods is
#' @export


hildaBarplot <- function(inputG, hildaResult, sigOrder=NULL, refGroup,
                         sortSampleNum=TRUE, refName="Control", altName="Case",
                         charSize=3) {
    numSig <- dim(hildaResult$BUGSoutput$mean$pStates1)[3]

    if(is(inputG)[1] != "MutationFeatureData") {
        stop("Not an output object from reading in the data using HiLDA.")
    }
    
    if(is(hildaResult) != "rjags") {
        stop("Not an output object from running the HiLDA tests.")
    }
    
    if (is.null(sigOrder)) {
        sigOrder <- seq_len(numSig)
    }

    if(length(sigOrder) != numSig) {
        stop("The order of signatures has more or less samples.")
    }
    
    if (length(refGroup) >= length(inputG@sampleList)) {
        stop(paste("There are more reference samples than the total samples"))
    }

    membership <- data.frame(sample=forcats::fct_reorder(inputG@sampleList,
                             seq_len(length(inputG@sampleList))),
                             rbind(hildaResult$BUGSoutput$mean$p1[, sigOrder],
                                   hildaResult$BUGSoutput$mean$p2[, sigOrder]))

    colnames(membership)[-1] <- paste0("Sig", seq_len(numSig))
    countData <- data.frame(t(inputG@countData))
    numMutations <- aggregate(countData$X3, by=list(Category=countData$X2),
                              FUN=sum)

    if (length(refGroup) > length(inputG@sampleList)) {
        stop(paste("The length of groups is greater than the sample size!"))
    }

    membership$grp <- altName
    membership$grp[refGroup] <- refName

    if (sortSampleNum == TRUE) {
        membership$sample <- forcats::fct_reorder2(membership$sample,
                                                   membership$grp,
                                                   numMutations$x, .desc=TRUE)
    }

    sigPlot <- hildaPlotSignature(hildaResult, sigOrder, charSize=charSize) +
        theme(plot.margin=grid::unit.c(unit(1, "lines") + unit(0.05, "npc"),
                                       unit(0, "npc"), unit(0.045, "npc"),
                                       unit(0, "npc")))

    propPlot <- ggplot(tidyr::gather_(membership, key_col = "sig", 
                                      value_col = "frac", 
                            gather_cols = c(paste0("Sig", seq_len(numSig)))),
                       aes(x=.data$sample, y=.data$frac, fill=.data$sig)) +
        geom_bar(stat="identity", width=0.8) + ylab("Proportions") +
        facet_grid(~grp, scales="free_x", space="free_x") +
        theme_bw() +
        theme(legend.position="none", 
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(), 
              axis.title.x=element_blank(),
              panel.border=element_blank(), 
              panel.grid.major.x=element_blank())

    return(list(sigPlot=sigPlot, propPlot=propPlot))
}



#' Read the raw mutation data with the mutation feature vector format,
#'   estimate and plot both mutation signatures and their fractions
#'
#' @param inputG a MutationFeatureData S4 class output by the pmsignature.
#' @param hildaResult a rjags class output by HiLDA.
#' @param sigOrder the order of signatures if needed (default: NULL).
#' @param charSize the size of the character on the signature plot (default: 3)
#' @return a list of the signature plot and the mean difference plot.
#'
#' @examples
#'
#' load(system.file("extdata/sample.rdata", package="HiLDA"))
#' inputFile <- system.file("extdata/hildaLocal.rdata", package="HiLDA")
#' hildaLocal <- readRDS(inputFile)
#'
#' hildaDiffPlot(G, hildaLocal)
#'
#' @importFrom cowplot plot_grid
#' @importFrom grid unit.c
#' @import ggplot2
#' @importFrom methods is
#'
#' @export


hildaDiffPlot <- function(inputG, hildaResult, sigOrder=NULL, charSize=3) {
    numSig <- dim(hildaResult$BUGSoutput$mean$pStates1)[3]

    if (is.null(sigOrder)) {
        sigOrder <- seq_len(numSig)
    }

    if(is(inputG)[1] != "MutationFeatureData") {
        stop("Not an output object from reading in the data using HiLDA.")
    }
    
    if(is(hildaResult) != "rjags") {
        stop("Not an output object from running the HiLDA tests.")
    }
    
    membership <- data.frame(signature=paste("signature", seq_len(numSig)),
                             hildaLocalResult(hildaResult)[sigOrder,
                                                           c(3, 5, 7)])
    colnames(membership)[-1] <- c("lower", "med", "upper")


    sigPlot <- hildaPlotSignature(hildaResult, sigOrder, charSize=charSize) +
        theme(plot.margin=grid::unit.c(unit(1, "lines") + unit(0.05, "npc"),
                                       unit(0, "npc"), unit(0.045, "npc"),
                                       unit(0, "npc")))

    diffPlot <- ggplot(membership,
                       aes(x=rev(.data$signature), y=.data$med, 
                           color=.data$signature)) +
        geom_pointrange(aes(ymin=.data$lower, ymax=.data$upper)) +
        coord_flip() +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        ylab("Difference in mean exposures (%) by HiLDA") +
        theme_minimal() +
        theme(legend.position="none", legend.title=element_blank(),
              axis.title.y=element_blank(), axis.text.y=element_blank(),
              axis.ticks.y=element_blank(), strip.text=element_blank())

    return(list(sigPlot=sigPlot, diffPlot=diffPlot))
}
