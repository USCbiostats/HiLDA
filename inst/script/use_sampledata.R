## load the library
library(HiLDA)

## import sample data
inputFile <- system.file("extdata/sample.rdata", package="HiLDA")
load(inputFile)

## run global and local test from HiLDA
hildaGlobal <- hildaGlobalTest(inputG=G, numSig=3, refGroup=1:4, nIter=2000)
hildaLocal <- hildaLocalTest(inputG=G, numSig=3, refGroup=1:4, nIter=2000)

## get signatures
K <- 3
Param <- pmsignature::getPMSignature(G, K=K)

## run global and local test from HiLDA
hildaGlobal <- hildaGlobalTest(inputG=G, numSig=3, useInits=Param,
                               refGroup=1:4, nIter=2000)
hildaLocal <- hildaLocalTest(inputG=G, numSig=3, useInits=Param,
                             refGroup=1:4, nIter=2000)

## visualize signatures from pmsignature
pmPlots <- pmBarplot(G, Param, refGroup=1:4)
cowplot::plot_grid(pmPlots$sigPlot, pmPlots$propPlot, rel_widths = c(1,3))

## visualize signatures from pmsignature
hildaPlots <- hildaBarplot(G, hildaLocal, refGroup=1:4)
cowplot::plot_grid(pmPlots$sigPlot, pmPlots$propPlot, rel_widths = c(1,3))

## examine the posterior distribution of mean differences
hildaGlobalResult(hildaGlobal)
hildaLocalResult(hildaLocal)

## visualize exposure differences from HiLDA
hildaDiffPlots <- hildaDiffPlot(G, hildaLocal)
cowplot::plot_grid(hildaDiffPlots$sigPlot, hildaDiffPlots$diffPlot,
                   rel_widths = c(1,3))


