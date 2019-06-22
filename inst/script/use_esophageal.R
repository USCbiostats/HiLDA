## load the library
library(HiLDA)

## read in the data
inputFile <- system.file("extdata/esophageal.mp.txt.gz", package="HiLDA")

## perform Shiraishi et al's method, pmsignature, to extract signatures
G <- hildaReadMPFile(inputFile, numBases=5, trDir=TRUE)

