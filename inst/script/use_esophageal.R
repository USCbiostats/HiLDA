## load the library
library(HiLDA)
library(pmsignature)

## import sample data
inputFile <- system.file("extdata/esophageal.mp.txt.gz", package="HiLDA")
G <- hildaReadMPFile(inputFile, numBases=5, trDir=TRUE)

## get signatures
K <- 4
Param <- getPMSignature(G, K=K)
