#' Plot mutation signatures from pmsignature output
#'
#' @param inputParam a estimatedParameters S4 class output by the pmsignature.
#' @param sigOrder the order of signatures if needed (default: NULL).
#' @param colorList a list of color to highlight the signatures (default: NULL).
#' @param ... Other arguments passed on to methods.
#' @return a plot object containing all mutational signatures
#'
#' @examples
#'
#' load(system.file("extdata/sample.rdata", package="HiLDA"))
#' Param <- pmgetSignature(G, K = 3)
#' pmPlotSignature(Param)
#'
#' @importFrom cowplot plot_grid
#' @importFrom grDevices hcl
#' @importFrom methods is
#' @export

pmPlotSignature <- function(inputParam, sigOrder=NULL, colorList=NULL, ...) {
    numSig <- inputParam@signatureNum
    numBases <- inputParam@flankingBasesNum

    if (is.null(sigOrder)) {
        sigOrder <- seq_len(numSig)
    }

    if (is(inputParam)[1] != "EstimatedParameters") {
        stop("The inputParam object is not an EstimatedParameters object")
    }
    
    if (is.null(colorList)) {
        colorList <- hcl(h=seq(15, 375, length=numSig + 1),
                         l=65, c=100)[seq_len(numSig)]
    }


    plotList <- vector("list", numSig)

    ## store the signature plots into a
    for (i in sigOrder) {
        plotList[[which(sigOrder == i)]] <-
            visPMS(inputParam@signatureFeatureDistribution[i, , ],
                   numBases=numBases,
                   trDir=inputParam@transcriptionDirection,
                   isScale=TRUE, ...) +
            theme(panel.border=
                      element_rect(color=colorList[which(sigOrder == i)],
                                   fill=NA, size=2, linetype = 1),
                  plot.margin=unit(c(0.01, 0.01, 0.01, 0.01), "npc"))
    }

    cowplot::plot_grid(plotlist=plotList, ncol=1)
}

#' Plot both mutation signatures and their mutational exposures from
#' pmsignature output
#'
#' @param inputG a MutationFeatureData S4 class output by the pmsignature.
#' @param inputParam a estimatedParameters S4 class output by the pmsignature.
#' @param sigOrder the order of signatures if needed (default: NULL).
#' @param refGroup the samples in the reference group (default: NULL).
#' @param sortSampleNum whether to sort by number of mutations (default: TRUE).
#' @param refName the name of reference group (default: Control).
#' @param altName the name of the other group (default: Case).
#' @param charSize the size of the character on the signature plot (default: 3).
#' @return a list of a signature plot and a barplot of mutational exposures
#'
#'
#' @examples
#'
#' load(system.file("extdata/sample.rdata", package="HiLDA"))
#' Param <- pmgetSignature(G, K = 3)
#'
#' pmPlots <- pmBarplot(G, Param, refGroup=1:4)
#' cowplot::plot_grid(pmPlots$sigPlot, pmPlots$propPlot, rel_widths = c(1,3))
#'
#'
#' @importFrom tidyr gather_
#' @importFrom grid unit.c
#' @importFrom forcats fct_relevel fct_reorder
#' @importFrom methods is
#' @export

pmBarplot <- function(inputG, inputParam, sigOrder=NULL, refGroup=NULL,
                      sortSampleNum=TRUE, refName="Control", altName="Case",
                      charSize=3) {
    numSig <- inputParam@signatureNum

    if (is.null(sigOrder)) {
        sigOrder <- seq_len(numSig)
    }
    
    if (is(inputParam)[1] != "EstimatedParameters") {
        stop("The inputParam object is not an EstimatedParameters object")
    }
    
    if(is(inputG)[1] != "MutationFeatureData") {
        stop("Not an output object from reading in the data using HiLDA.")
    }
    
    membership <- data.frame(sample=forcats::fct_reorder(inputParam@sampleList,
                                    seq_len(length(inputParam@sampleList))),
                             inputParam@sampleSignatureDistribution[, sigOrder])
    colnames(membership)[-1] <- paste0("Sig", seq_len(numSig))
    countData <- data.frame(t(inputG@countData))
    numMutations <- aggregate(countData$X3, by=list(Category=countData$X2),
                              FUN=sum)

    if (is.null(refGroup)) {
        if (sortSampleNum == TRUE) {
            membership$sample <- forcats::fct_reorder(membership$sample,
                                                      numMutations$x,
                                                      .desc=TRUE)
        }
        sigPlot <- pmPlotSignature(inputParam, sigOrder, charSize=charSize) +
            theme(plot.margin=unit(c(0.045, 0, 0.045, 0), "npc"))

        propPlot <- ggplot(tidyr::gather_(membership, key_col = "sig", 
                                          value_col = "frac", 
                                          gather_cols = 
                            c(paste0("Sig",seq_len(inputParam@signatureNum)))),
                           aes(x=.data$sample, y=.data$frac, fill=.data$sig)) +
            geom_bar(stat="identity", width=0.8) + ylab("Proportions") +
            theme_bw() +
            theme(legend.position="none",
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.title.x=element_blank(),
                  panel.border=element_blank(),
                  panel.grid.major.x=element_blank())
    } else {
        if (length(refGroup) > length(inputParam@sampleList)) {
            stop(paste("The length of groups is greater than the sample size!"))
        }

        membership$grp <- altName
        membership$grp[refGroup] <- refName

        if (sortSampleNum == TRUE) {
            membership$sample <- forcats::fct_reorder2(membership$sample,
                                                       membership$grp,
                                                       numMutations$x,
                                                       .desc=TRUE)
        }

        sigPlot <- pmPlotSignature(inputParam, sigOrder, charSize=charSize) +
            theme(plot.margin=grid::unit.c(unit(1, "lines") + unit(0.05, "npc"),
                                           unit(0, "npc"),
                                           unit(0.045, "npc"), unit(0, "npc")))

        propPlot <- ggplot(tidyr::gather_(membership, key_col = "sig", 
                                          value_col = "frac", 
                                          gather_cols = 
                            c(paste0("Sig",seq_len(numSig)))),
                           aes(x=.data$sample, y=.data$frac, fill=.data$sig)) +
            geom_bar(stat="identity", width=0.8) + ylab("Proportions") +
            facet_grid(~grp, scales="free_x", space="free_x") +
            theme_bw() + theme(legend.position="none",
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            panel.border=element_blank(),
            panel.grid.major.x=element_blank())
    }

    return(list(sigPlot=sigPlot, propPlot=propPlot))
}



#' Plot both mutation signatures and their mutational exposures from
#' pmsignature output for more than two groups
#'
#' @param inputG a MutationFeatureData S4 class output by the pmsignature.
#' @param inputParam a estimatedParameters S4 class output by the pmsignature.
#' @param sigOrder the order of signatures if needed (default: NULL).
#' @param groupIndices a vector of group indicators (default: NULL).
#' @param sortSampleNum an indictor variable on whether samples are sorted by
#'        the number of mutations (default: TRUE).
#' @param charSize the size of the character on the signature plot (default: 3)
#' @return a list of the signature plot and the mean difference plot.
#'
#' @examples
#'
#' load(system.file("extdata/sample.rdata", package="HiLDA"))
#' Param <- pmgetSignature(G, K = 3)
#'
#' pmPlots <- pmMultiBarplot(G, Param, groupIndices=c(1, rep(2,3), rep(3,6)))
#' cowplot::plot_grid(pmPlots$sigPlot, pmPlots$propPlot, rel_widths = c(1,3))
#'
#' @importFrom cowplot plot_grid
#' @importFrom tidyr gather_
#' @importFrom grid unit.c
#' @importFrom forcats fct_relevel fct_reorder
#' @importFrom methods is
#' @export

pmMultiBarplot <- function(inputG, inputParam, sigOrder=NULL, groupIndices,
                           sortSampleNum=TRUE, charSize=3) {
    numSig <- inputParam@signatureNum

    if (is.null(sigOrder)) {
        sigOrder <- seq_len(numSig)
    }

    if(is(inputG)[1] != "MutationFeatureData") {
        stop("Not an output object from reading in the data using HiLDA.")
    }
    
    if (is(inputParam)[1] != "EstimatedParameters") {
        stop("The inputParam object is not an EstimatedParameters object")
    }
    
    if (length(unique(groupIndices)) == 1) {
        stop(paste("More than one group is required!"))
    }

    membership <- data.frame(sample=forcats::fct_reorder(inputParam@sampleList,
                            seq_len(length(inputParam@sampleList))),
                            inputParam@sampleSignatureDistribution[, sigOrder])
    colnames(membership)[-1] <- paste0("Sig", seq_len(numSig))
    countData <- data.frame(t(inputG@countData))
    numMutations <- aggregate(countData$X3, by=list(Category=countData$X2),
                              FUN=sum)

    if (length(groupIndices) > length(inputParam@sampleList)) {
        stop(paste("The length of groups is greater than the sample size!"))
    }

    membership$grp <- groupIndices

    if (sortSampleNum == TRUE) {
        membership$sample <- forcats::fct_reorder2(membership$sample,
                                                   membership$grp,
                                                   numMutations$x,
                                                   .desc=TRUE)
    }

    sigPlot <- pmPlotSignature(inputParam, sigOrder, charSize=charSize) +
        theme(plot.margin=grid::unit.c(unit(1, "lines") + unit(0.05, "npc"),
                                       unit(0, "npc"), unit(0.045, "npc"),
                                       unit(0, "npc")))

    propPlot <- ggplot(tidyr::gather_(membership, key_col = "sig", 
                                      value_col = "frac", 
                        c(paste0("Sig",seq_len(numSig)))),
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
