## ----style, echo = FALSE, results = 'asis'---------------------------------
library(BiocStyle)

## --------------------------------------------------------------------------
library(HiLDA)
inputFile <- system.file("extdata/esophageal.mp.txt.gz", package="HiLDA")
G <- hildaReadMPFile(inputFile, numBases=5, trDir=TRUE)

## --------------------------------------------------------------------------
load(system.file("extdata/sample.rdata", package = "HiLDA"))
class(G)

## --------------------------------------------------------------------------
set.seed(123)
hildaGlobal <- hildaTest(inputG=G, numSig=3, localTest=FALSE, 
                         refGroup=1:4, nIter=1000)
hildaLocal <- hildaTest(inputG=G, numSig=3, localTest=TRUE, 
                        refGroup=1:4, nIter=1000)

## --------------------------------------------------------------------------
Param <- pmgetSignature(G, K = 3)

## --------------------------------------------------------------------------
set.seed(123)
hildaGlobal <- hildaTest(inputG=G, numSig=3, useInits = Param,
                         localTest=TRUE, refGroup=1:4, nIter=1000)
hildaLocal <- hildaTest(inputG=G, numSig=3, useInits = Param,
                        localTest=TRUE, refGroup=1:4, nIter=1000)

## --------------------------------------------------------------------------
hildaRhat(hildaGlobal)
hildaRhat(hildaLocal)

## --------------------------------------------------------------------------
pmPlots <- pmBarplot(G, Param, refGroup=1:4, sigOrder=c(1,3,2))
cowplot::plot_grid(pmPlots$sigPlot, pmPlots$propPlot, rel_widths = c(1,3))

