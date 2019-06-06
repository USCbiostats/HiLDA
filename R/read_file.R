#' Read the raw mutation data of Mutation Position Format.
#'
#' @description
#' The mutation position format is tab-delimited text file, where
#' the 1st-5th columns shows sample names, chromosome names,
#' coordinates, reference bases (A, C, G, or T) and
#' the alternate bases (A, C, G, or T), respectively. An example is as follows;
#'
#' ---
#'
#' sample1 chr1 100 A C
#'
#' sample1 chr1 200 A T
#'
#' sample1 chr2 100 G T
#'
#' sample2 chr1 300 T C
#'
#' sample3 chr3 400 T C
#'
#' ---
#'
#' Also, this function usually can accept compressed files (e.g., by gzip, bzip2
#' and so on) when using recent version of R.
#' Currently, only UCSC hg19 (BSgenome.Hsapiens.UCSC.hg19) is supported.
#'
#' @param infile the path for the input file for the mutation data of Mutation
#' Position Format.
#' @param numBases the number of upstream and downstream flanking bases
#' (including the mutated base) to take into account.
#' @param trDir the index representing whether transcription direction is
#' considered or not.
#' The gene annotation information is given by UCSC knownGene
#' (TxDb.Hsapiens.UCSC.hg19.knownGene object) When trDir is TRUE, the mutations
#' located in intergenic region are excluded from the analysis.
#' @param type this argument can take either "independent", "full", or "custom".
#' @param bs_genome this argument specifies the reference genome (e.g., B
#' Sgenome.Mmusculus.UCSC.mm10 can be used for the mouse genome).
#' See https://bioconductor.org/packages/release/bioc/html/BSgenome.html for the
#' available genome list
#' @param txdb_transcript this argument specified the transcript database
#' (e.g., TxDb.Mmusculus.UCSC.mm10.knownGene can be used for the mouse genome).
#' See https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html
#' for details.
#' @return The output is an instance of MutationFeatureData S4 class (which
#' stores summarized information on mutation data). This will be typically used
#' as the initial values for the global test and the local test.
#'
#' @importFrom utils read.table
#' @examples
#' inputFile <- system.file("extdata/esophageal.mp.txt.gz", package="HiLDA")
#' G <- hildaReadMPFile(inputFile, numBases=5, trDir=TRUE)
#'
#'
#' @export
hildaReadMPFile <- function(infile, numBases = 3, trDir = FALSE,
                            type = "independent", bs_genome = NULL,
                            txdb_transcript = NULL) {

    if (is.null(bs_genome) == TRUE) {
      bs_genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    }

    if (is.null(txdb_transcript) == TRUE) {
      txdb_transcript <-
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    }

    if (type == "independent") {
      fdim <- c(6, rep(4, numBases - 1), rep(2, as.integer(trDir)))
    } else if (type == "full") {
      fdim <- c(6 * 4^(numBases - 1) * 2^(as.integer(trDir)))
    } else {
      stop('for reading mutation position format, the type argument has to be
           "independent" or "full"')
    }

    if (numBases %% 2 != 1) {
      stop("numBases should be odd numbers")
    }
    centerInd <- (numBases + 1) / 2

    mutFile <- read.table(infile, sep="\t", header=FALSE,
                          stringsAsFactors = FALSE)

    chrInfo <- mutFile[,2]
    posInfo <- mutFile[,3]
    ref_base <- Biostrings::DNAStringSet(mutFile[,4])
    alt_base <- Biostrings::DNAStringSet(mutFile[,5])
    sampleName_str <- as.character(mutFile[,1])

    gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chrInfo,
                                                             start = posInfo,
                                                             end = posInfo),
                                                  ignore.strand = TRUE)

    ranges <- GenomicRanges::resize(gr, numBases, fix = "center")
    context <- Biostrings::getSeq(bs_genome, ranges)


    ## check the consistency between the input reference base and the obtained
    ## base using hg19 reference genome.
    removeInd <- which(XVector::subseq(context, start = centerInd,
                                       end = centerInd) != ref_base)
    if (length(removeInd) > 0) {
      warning(paste("The central bases are inconsistent in", length(removeInd),
                    "mutations. We have removed them."))
      context <- context[-removeInd]
      ref_base <- ref_base[-removeInd]
      alt_base <- alt_base[-removeInd]
      sampleName_str <- sampleName_str[-removeInd]
      chrInfo <- chrInfo[-removeInd]
      posInfo <- posInfo[-removeInd]
    }

    # check the characters on alternative base
    alphabetFreq <- Biostrings::alphabetFrequency(alt_base)
    removeInd <- which(rowSums(alphabetFreq[,seq_len(4)]) != 1)
    if (length(removeInd) > 0) {
      warning(paste("The characters other than (A, C, G, T) are included in
                    alternate bases of", length(removeInd),
                    "mutations. We have removed them."))
      context <- context[-removeInd]
      ref_base <- ref_base[-removeInd]
      alt_base <- alt_base[-removeInd]
      sampleName_str <- sampleName_str[-removeInd]
      chrInfo <- chrInfo[-removeInd]
      posInfo <- posInfo[-removeInd]
    }

    # check the characters on flanking bases
    alphabetFreq <- Biostrings::alphabetFrequency(context)
    removeInd <- which(alphabetFreq[,"A"] + alphabetFreq[,"C"] +
                         alphabetFreq[,"G"] + alphabetFreq[,"T"] != numBases)
    if (length(removeInd) > 0) {
      context <- context[-removeInd]
      ref_base <- ref_base[-removeInd]
      alt_base <- alt_base[-removeInd]
      sampleName_str <- sampleName_str[-removeInd]
      chrInfo <- chrInfo[-removeInd]
      posInfo <- posInfo[-removeInd]
      warning(paste("The characters other than (A, C, G, T) are included in
                    flanking bases of", length(removeInd),
                    "mutations. We have removed them."))
    }

    # check the characters on alternative base
    alphabetFreq <- Biostrings::alphabetFrequency(alt_base)
    removeInd <- which(ref_base == alt_base)
    if (length(removeInd) > 0) {
      warning(paste("The reference base and alternative bases are equal for",
                    length(removeInd), "mutations. We have removed them."))
      context <- context[-removeInd]
      ref_base <- ref_base[-removeInd]
      alt_base <- alt_base[-removeInd]
      sampleName_str <- sampleName_str[-removeInd]
      chrInfo <- chrInfo[-removeInd]
      posInfo <- posInfo[-removeInd]
    }


    revCompInd <- which(as.character(XVector::subseq(context,
                                                     start = centerInd,
                                                     end = centerInd)) %in%
                          c("A", "G"))
    context[revCompInd] <- Biostrings::reverseComplement(context[revCompInd])
    ref_base[revCompInd] <- Biostrings::reverseComplement(ref_base[revCompInd])
    alt_base[revCompInd] <- Biostrings::reverseComplement(alt_base[revCompInd])

    # Obtaining transcription strand information using GenomicRanges packages
    if (trDir == TRUE) {
      gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chrInfo,
                                                               start = posInfo,
                                                               end = posInfo),
                                                    ignore.strand = TRUE)
      txdb <- txdb_transcript

      txdb_bed <- GenomicFeatures::transcripts(txdb)
      gr_txdb <- GenomicRanges::findOverlaps(gr, txdb_bed,
                                             ignore.strand = FALSE)
      gr_strand <- cbind(
        S4Vectors::queryHits(gr_txdb),
        as.character(S4Vectors::as.factor(
          BiocGenerics::strand(txdb_bed[S4Vectors::subjectHits(gr_txdb)]))))
      ugr_strand <- unique(gr_strand[gr_strand[, 2] != "*" ,], MARGIN=1)

      rmdup_ugr_strand <- ugr_strand[!duplicated(ugr_strand[, 1]), ]
      txdb_plus_gr_ind <- as.integer(
        rmdup_ugr_strand[rmdup_ugr_strand[, 2] == "+", 1])
      txdb_minus_gr_ind <-
        as.integer(rmdup_ugr_strand[rmdup_ugr_strand[, 2] == "-", 1])

      strandInfo <- rep("*", length(gr))
      strandInfo[setdiff(txdb_plus_gr_ind, revCompInd)] <- "+"
      strandInfo[intersect(txdb_plus_gr_ind, revCompInd)] <- "-"
      strandInfo[setdiff(txdb_minus_gr_ind, revCompInd)] <- "-"
      strandInfo[intersect(txdb_minus_gr_ind, revCompInd)] <- "+"

      warning(paste("Out of", length(context), "mutations, we could obtain
                    transcription direction information for",
                    length(txdb_plus_gr_ind) + length(txdb_minus_gr_ind),
                    "mutation. Other mutations are removed."))

      context <- context[strandInfo != "*"]
      ref_base <- ref_base[strandInfo != "*"]
      alt_base <- alt_base[strandInfo != "*"]
      sampleName_str <- sampleName_str[strandInfo != "*"]
      chrInfo <- chrInfo[strandInfo != "*"]
      posInfo <- posInfo[strandInfo != "*"]
      strandInfo <- strandInfo[strandInfo != "*"]
    }

    if (trDir == FALSE) {
      strandInfo <- NULL
    }

    mutFeatures <- getMutationFeatureVector(context, ref_base, alt_base,
                                            strandInfo, numBases, type)

    suSampleStr <- sort(unique(sampleName_str))
    lookupSampleInd <- seq_len(length(suSampleStr))
    names(lookupSampleInd) <- suSampleStr
    sampleIDs <- lookupSampleInd[sampleName_str]


    featStr <- apply(mutFeatures, 1, paste0, collapse=",")

    suFeatStr <- sort(unique(featStr))
    lookupFeatInd <- seq_len(length(suFeatStr))
    names(lookupFeatInd) <- suFeatStr

    rawCount <- data.frame(sample = sampleIDs, mutInds = lookupFeatInd[featStr])

    tableCount <- table(rawCount)
    w <- which(tableCount > 0, arr.ind=TRUE)
    procCount <- cbind(w[,2], w[,1], tableCount[w])


    if (length(fdim) == 1) {
      mutFeatList <- matrix(as.integer(suFeatStr), length(suFeatStr), 1)
    } else {
      mutFeatList <- t(vapply(suFeatStr,
                              function(x) as.numeric(unlist(strsplit(x, ","))),
                              numeric(numBases + as.integer(trDir))))
    }

    rownames(mutFeatList) <- NULL
    rownames(procCount) <- NULL

    if (trDir == FALSE) {
      strandInfo_for_class <- rep(NA, length(chrInfo))
    } else {
      strandInfo_for_class <- strandInfo
    }

    return(new(Class = "MutationFeatureData",
               type = type,
               flankingBasesNum = as.integer(numBases),
               transcriptionDirection = trDir,
               possibleFeatures = as.integer(fdim),
               featureVectorList = t(mutFeatList),
               sampleList = suSampleStr,
               countData = t(procCount),
               mutationPosition =
                 data.frame(chr = chrInfo, pos = posInfo,
                            ref = ref_base, alt = alt_base,
                            strand = strandInfo_for_class,
                            context = context,
                            sampleID = unname(sampleIDs),
                            mutID = unname(lookupFeatInd[featStr]),
                            stringsAsFactors = FALSE)
    )
    )

  }



#' Get mutation feature vector from context sequence data and reference and
#' alternate allele information
#'
#' @param context the context sequence data around the mutated position. This
#' shoud be Biostrings::DNAStringSet class
#' @param ref_base the reference bases at the mutated position.
#' @param alt_base the alternate bases at the mutated position.
#' @param strandInfo transcribed strand information at the mutated position.
#'   (this is optional)
#' @param numBases the number of flanking bases around the mutated position.
#' @param type the type of mutation feature vecotr (should be "independent" or
#' "full").
#' @return a mutation featuer vector
#'
#'
#'
#' @export
getMutationFeatureVector <- function(context, ref_base, alt_base,
                                     strandInfo = NULL, numBases, type) {

    trDir <- !is.null(strandInfo)
    if (type == "independent") {
      fdim <- c(6, rep(4, numBases - 1), rep(2, as.integer(trDir)))
    } else if (type == "full") {
      fdim <- c(6 * 4^(numBases - 1) * 2^(as.integer(trDir)))
    } else {
      stop('the type argument has to be "independent" or "full"')
    }

    if (numBases %% 2 != 1) {
      stop("numBases should be odd numbers")
    }
    centerInd <- (numBases + 1) / 2

    if (type == "independent") {

      mutFeatures <- matrix(0, length(ref_base), length(fdim))

      centralBaseList <- data.frame(ref_base=rep(c("C", "T"), each=3),
                                    alt_base=c("A", "G", "T", "A", "C", "G"),
                                    stringsAsFactors = FALSE)

    ## rewrite the part from pmsignature

    for (i in seq_len(6)) {
      mutFeatures[which(ref_base == centralBaseList$ref_base[i] &
                        alt_base == centralBaseList$alt_base[i]), 1] <- i
    }

    columnInd <- 2
    for (baseInd in seq_len(numBases)) {
      if (baseInd == centerInd) {
        next
      }

      baseList <- c("A", "C", "G", "T")

      ## rewrite the part from pmsignature
      for (base in seq_along(baseList)) {
        mutFeatures[which(XVector::subseq(context,
                                          start=baseInd,
                                          end=baseInd) == baseList[base]),
                    columnInd] <- base
      }

      columnInd <- columnInd + 1
    }

    if (trDir ==TRUE) {
      mutFeatures[which(strandInfo == "+"), length(fdim)] <- 1
      mutFeatures[which(strandInfo == "-"), length(fdim)] <- 2
    }

    } else {

    mutFeatures <- matrix(1, length(ref_base), length(fdim))

    tempDigits <- 1
    for (i in seq_len((numBases - 1) / 2)) {

      baseInd <- numBases + 1 - i

      mutFeatures[which(XVector::subseq(context, start=baseInd, end=baseInd)
                        == "C"), 1] <-
        mutFeatures[which(XVector::subseq(context, start=baseInd, end=baseInd)
                          == "C"), 1] + tempDigits * 1

      mutFeatures[which(XVector::subseq(context, start=baseInd, end=baseInd)
                        == "G"), 1] <-
        mutFeatures[which(XVector::subseq(context, start=baseInd, end=baseInd)
                          == "G"), 1] + tempDigits * 2

      mutFeatures[which(XVector::subseq(context, start=baseInd, end=baseInd)
                        == "T"), 1] <-
        mutFeatures[which(XVector::subseq(context, start=baseInd, end=baseInd)
                          == "T"), 1] + tempDigits * 3

      tempDigits <- tempDigits * 4
      baseInd <- i;
      mutFeatures[which(XVector::subseq(context, start=baseInd, end=baseInd)
                        == "C"), 1] <-
        mutFeatures[which(XVector::subseq(context, tart=baseInd, end=baseInd)
                          == "C"), 1] + tempDigits * 1

      mutFeatures[which(XVector::subseq(context, start=baseInd, end=baseInd)
                        == "G"), 1] <-
        mutFeatures[which(XVector::subseq(context, start=baseInd, end=baseInd)
                          == "G"), 1] + tempDigits * 2

      mutFeatures[which(XVector::subseq(context, start=baseInd, end=baseInd)
                        == "T"), 1] <-
        mutFeatures[which(XVector::subseq(context, start=baseInd, end=baseInd)
                          == "T"), 1] + tempDigits * 3

      tempDigits <- tempDigits * 4
    }

    mutFeatures[which(ref_base == "C" & alt_base == "G"), 1] <-
      mutFeatures[which(ref_base == "C" & alt_base == "G"), 1] +
      tempDigits * 1
    mutFeatures[which(ref_base == "C" & alt_base == "T"), 1] <-
      mutFeatures[which(ref_base == "C" & alt_base == "T"), 1] +
      tempDigits * 2
    mutFeatures[which(ref_base == "T" & alt_base == "A"), 1] <-
      mutFeatures[which(ref_base == "T" & alt_base == "A"), 1] +
      tempDigits * 3
    mutFeatures[which(ref_base == "T" & alt_base == "C"), 1] <-
      mutFeatures[which(ref_base == "T" & alt_base == "C"), 1] +
      tempDigits * 4
    mutFeatures[which(ref_base == "T" & alt_base == "G"), 1] <-
      mutFeatures[which(ref_base == "T" & alt_base == "G"), 1] +
      tempDigits * 5

    if (trDir ==TRUE) {
      tempDigits <- tempDigits * 6
      mutFeatures[which(strandInfo == "-"), 1] <-
        mutFeatures[which(strandInfo == "-"), 1] + tempDigits * 1
    }

  }

  return(mutFeatures)

}

