#' Apply HiLDA to statistically testing the burdens of mutation signatures
#' between two groups
#'
#' @param inputG a MutationFeatureData S4 class output by the pmsignature.
#' @param numSig an integer number of the number of mutational signatures.
#' @param refGroup the reference group.
#' @param useInits a EstimatedParameters S4 class output by the pmsignature
#'        (default: NULL)
#' @param sigOrder the indice indicating the samples in the reference group
#'        (default: NULL)
#' @param nIter number of total iterations per chain (default: 2000).
#' @param nBurnin length of burn (default: 0).
#' @param ... Other arguments passed on to methods.
#' @return the output jags file
#' @examples
#'
#' load(system.file("extdata/sample.rdata", package="HiLDA"))
#'
#' ## with initial values
#' hildaLocal <- hildaLocalTest(inputG=G, numSig=3, refGroup=1:4, nIter=1000)
#'
#' @importFrom R2jags jags
#' @export


hildaLocalTest <- function(inputG, numSig, refGroup, useInits=NULL,
                           sigOrder=NULL, nIter=2000, nBurnin=0, ...) {
    if (is.null(sigOrder)) {
        sigOrder <- seq_len(numSig)
    }

    # reshape the counts of input data
    countWide <- as.data.frame(t(inputG@countData))
    colnames(countWide) <- c("type", "sample", "count")
    countLong <- reshape(countWide, idvar="sample", timevar="type",
                         direction="wide")
    countLong <- countLong[order(countLong[, 1]), ]
    countLong[is.na(countLong)] <- 0
    countLong <- countLong[, -1]

    # generate the known data for MCMC
    Num1 <- length(refGroup)
    Num2 <- length(inputG@sampleList) - Num1
    rownames(countLong) <- seq_len(Num1 + Num2)
    mutationN <- as.vector(rowSums(countLong))
    caseGroup <- setdiff(seq_len(Num1 + Num2), refGroup)

    nFlanking <- inputG@flankingBasesNum - 1
    nFeature <- length(inputG@possibleFeatures)

    xG1 <- structure(.Data=rep(0, Num1 * max(mutationN[refGroup]) * nFeature),
                     .Dim=c(max(mutationN[refGroup]), nFeature, Num1))
    xG2 <- structure(.Data=rep(0, Num2 * max(mutationN[caseGroup]) * nFeature),
                     .Dim=c(max(mutationN[caseGroup]), nFeature, Num2))

    for (i in refGroup) {
        xG1[seq_len(sum(countLong[i, ])), , which(refGroup == i)] <-
    t(inputG@featureVectorList[, rep(seq_len(ncol(inputG@featureVectorList)),
                                             countLong[i, ])])
    }

    for (i in caseGroup) {
        xG2[seq_len(sum(countLong[i, ])), , which(caseGroup == i)] <-
    t(inputG@featureVectorList[, rep(seq_len(ncol(inputG@featureVectorList)),
                                             countLong[i, ])])
    }

    ## set up the MCMC
    N1 <- mutationN[refGroup]
    N2 <- mutationN[caseGroup]

    jdata <- list(I1=Num1, I2=Num2, K=numSig, numStates=6,
                  numFlank=nFlanking, N1=N1, N2=N2,
                  xG1=xG1, xG2=xG2)

    varSaveTr <- c("pStates1", "pStates2", "pStates3",
                   "p1", "p2", "alpha", "beta")
    varSave <- c("pStates1", "pStates2", "p1", "p2", "alpha", "beta")

    # initial values for signatures
    if (is.null(useInits) == FALSE) {
        sig <- abind::abind(lapply(sigOrder, function(x)
            useInits@signatureFeatureDistribution[x, , ]), along=3)
        if(useInits@transcriptionDirection == TRUE){
            inits_tr <- list(list(pStates1=array(sig[1, , sigOrder],
                                                 dim=c(1, 6, numSig)),
                                  pStates2=array(sig[2:nFeature, seq_len(4),
                                      sigOrder], dim=c(nFlanking, 4, numSig)),
                                  pStates3=array(sig[nFeature, seq_len(2),
                                      sigOrder], dim=c(1, 2, numSig))),
                             list(pStates1=array(sig[1, , sigOrder],
                                                 dim=c(1, 6, numSig)),
                                  pStates2=array(sig[2:nFeature, seq_len(4),
                                      sigOrder], dim=c(nFlanking, 4, numSig)),
                                  pStates3=array(sig[nFeature, seq_len(2),
                                      sigOrder], dim=c(1, 2, numSig))))
        } else {
            inits <- list(list(pStates1=array(sig[1, , sigOrder],
                                              dim=c(1, 6, numSig)),
                               pStates2=array(sig[2:nFeature, seq_len(4),
                                   sigOrder], dim=c(nFlanking, 4, numSig))),
                          list(pStates1=array(sig[1, , sigOrder],
                                              dim=c(1, 6, numSig)),
                               pStates2=array(sig[2:nFeature, seq_len(4),
                                   sigOrder], dim=c(nFlanking, 4, numSig))))
        }
    }


    if (is.null(useInits) == FALSE) {
        ## which means using inits from Param
        if (inputG@transcriptionDirection == TRUE) {
            modelFit <- R2jags::jags(model.file=system.file("local_tr.txt",
                package="HiLDA"), data=jdata, parameters.to.save=varSaveTr,
                inits=inits_tr, n.chains=2, n.iter=nIter, n.burnin=nBurnin)
        } else {
            modelFit <- R2jags::jags(model.file=system.file("local.txt",
                package="HiLDA"), data=jdata, parameters.to.save=varSave,
                inits=inits, n.chains=2, n.iter=nIter, n.burnin=nBurnin)
        }
    } else {
        if (inputG@transcriptionDirection == TRUE) {
            modelFit <- R2jags::jags(model.file=system.file("local_tr.txt",
                package="HiLDA"), data=jdata, parameters.to.save=varSaveTr,
                n.chains=2, n.iter=nIter, n.burnin=nBurnin)
        } else {
            modelFit <- R2jags::jags(model.file=system.file("local.txt",
                package="HiLDA"), data=jdata, parameters.to.save=varSave,
                n.chains=2, n.iter=nIter, n.burnin=nBurnin)
        }
    }

    return(modelFit)
}
