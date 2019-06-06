#' @title visualize probabisitic mutaiton signature for the independent model
#' @description Generate visualization of mutation signatures for the model with
#'              substitution patterns and flanking bases represented by the
#'              indepenent representation.
#'
#' @param vF a matrix for mutation signature
#' @param numBases the number of flanking bases
#' @param baseCol the colour of the bases (A, C, G, T, plus/minus strand)
#' @param trDir the index whether the strand direction is plotted or not
#' @param charSize the size of the character
#' @param isScale the index whether the height of the flanking base is changed
#'        or not
#' @param alpha the parameter for the Renyi entropy (applicable only if the
#'        isScale is TRUE)
#' @param charLimit the limit of char size
#' @return a plot of the input mutational signature
#'
#' @import ggplot2
#'
#' @examples
#'
#' load(system.file("extdata/sampleParam.rdata", package="HiLDA"))
#'
#' sig <- slot(Param, "signatureFeatureDistribution")[1,,]
#' visPMS(sig, numBases = 5, isScale = TRUE)
#'
#'
#' @export
visPMS <- function(vF, numBases, baseCol=NA, trDir=FALSE, charSize=5,
                       isScale=FALSE, alpha=2, charLimit=0.25) {

    if (is.na(baseCol)) {
      gg_color_hue6 <- hcl(h=seq(15, 375, length=7), l=65, c=100)[seq_len(6)]
      baseCol <- c(gg_color_hue6[c(3, 5, 2, 1, 6)])
    }

    centerBase <- (1 + numBases) / 2

    v1 <- vF[1,seq_len(6)]
    V2 <- vF[2:(numBases),seq_len(4)]
    A <- matrix(0, numBases, 4)
    B <- matrix(0, 4, 4)

    if (trDir == TRUE) {
      v3 <- vF[(numBases + 1),seq_len(2)]
    }

    for (l in seq_len(numBases)) {
      if (l < centerBase) {
        A[l, ] <- V2[l, ]
      } else if (l > centerBase) {
        A[l, ] <- V2[l - 1, ]
      }
    }
    A[centerBase,2] <- sum(v1[seq_len(3)])
    A[centerBase,4] <- sum(v1[4:6])

    B[2, c(1, 3, 4)] <- v1[seq_len(3)] / sum(v1[seq_len(3)])
    B[4, c(1, 2, 3)] <- v1[4:6] / sum(v1[4:6])

    renyi <- function(p, tAlpha=alpha) {
      if (tAlpha == 1) {
        return(- sum(p * log2(p), na.rm=TRUE))
      } else {
        return( log(sum(p^tAlpha)) / (1 - tAlpha))
      }
    }

    if (isScale == FALSE) {
      fheight <- rep(1, numBases)
    } else {
      fheight <- 0.5 * (2 - apply(A, MARGIN=1, FUN=renyi))
    }

    ##  collecting data for ggplot
    x_start <- c()
    x_end <- c()
    y_start <- c()
    y_end <- c()
    text_x <- c()
    text_y <- c()
    text_lab <- c()
    text_col <- c()
    rectType <- c()
    num2base <- c("A", "C", "G", "T")

    # flanking bases
    tempStartX <- 0
    for (i in seq_len(numBases)) {
      x_start <- c(x_start, tempStartX + c(0, cumsum(A[i,seq_len(3)])))
      x_end <- c(x_end, tempStartX + cumsum(A[i,seq_len(4)]))
      y_start <- c(y_start, rep(0, 4))
      y_end <- c(y_end, rep(fheight[i], 4))
      rectType <- c(rectType, c("A", "C", "G", "T"))
      for (j in seq_len(4)) {
        tempPos <- c(0, cumsum(A[i,seq_len(4)]))
        if (A[i,j] > charLimit && fheight[i] > charLimit) {
          text_x <- c(text_x, tempStartX + 0.5 * (tempPos[j] + tempPos[j + 1]))
          text_y <- c(text_y, 0.5 * (0 + fheight[i]))
          text_lab <- c(text_lab, num2base[j])
          text_col <- c(text_col, "w")
        }
      }
      tempStartX <- tempStartX + 1.25
    }

    ## alternative bases from C
    tempStartX <- (centerBase - 1) * 1.25
    x_start <- c(x_start, rep(tempStartX, 4))
    x_end <- c(x_end, rep(tempStartX + A[centerBase, 2], 4))
    y_start <- c(y_start, 2 + c(0, cumsum(B[2,seq_len(3)])))
    y_end <- c(y_end, 2 + cumsum(B[2,seq_len(4)]))
    rectType <- c(rectType, c("A", "C", "G", "T"))

    tempPos <- c(0, cumsum(B[2,seq_len(4)]))
    for (j in seq_len(4)) {
      if (A[centerBase, 2] > charLimit && B[2,j] > charLimit) {
        text_x <- c(text_x, tempStartX + 0.5 * A[centerBase, 2])
        text_y <- c(text_y, 2 + 0.5 * (tempPos[j] + tempPos[j + 1]))
        text_lab <- c(text_lab, num2base[j])
        text_col <- c(text_col, "w")
      }
    }

    ## alternative bases from T
    tempStartX <- tempStartX + A[centerBase, 2]
    x_start <- c(x_start, rep(tempStartX, 4))
    x_end <- c(x_end, rep(tempStartX + A[centerBase, 4], 4))
    y_start <- c(y_start, 2 + c(0, cumsum(B[4,seq_len(3)])))
    y_end <- c(y_end, 2 + cumsum(B[4,seq_len(4)]))
    rectType <- c(rectType, c("A", "C", "G", "T"))

    tempPos <- c(0, cumsum(B[4,seq_len(4)]))
    for (j in seq_len(4)) {
      if (A[centerBase, 4] > charLimit && B[4,j] > charLimit) {
        text_x <- c(text_x, tempStartX + 0.5 * A[centerBase, 4])
        text_y <- c(text_y, 2 + 0.5 * (tempPos[j] + tempPos[j + 1]))
        text_lab <- c(text_lab, num2base[j])
        text_col <- c(text_col, "w")
      }
    }

    if (trDir == TRUE) {

      # draw direction bias
      x_start <- c(x_start, (numBases - 1) * 1.25 + 0.24)
      x_end <- c(x_end, (numBases - 1) * 1.25 + 0.49)
      y_start <- c(y_start, 2)
      y_end <- c(y_end, 2 + v3[1])
      rectType <- c(rectType, c("+"))

      if (v3[1] > 0.125) {
        text_x <- c(text_x, (numBases - 1) * 1.25 + 0.5 * (0.24 + 0.49))
        text_y <- c(text_y, 2 + 0.5 * v3[1])
        text_lab <- c(text_lab, "+")
        text_col <- c(text_col, "w")
      }


      x_start <- c(x_start, (numBases - 1) * 1.25 + 0.51)
      x_end <- c(x_end, (numBases - 1) * 1.25 + 0.76)
      y_start <- c(y_start, 2)
      y_end <- c(y_end, 2 + v3[2])
      rectType <- c(rectType, c("-"))

      if (v3[2] > 0.125) {
        text_x <- c(text_x, (numBases - 1) * 1.25 + 0.5 * (0.51 + 0.76))
        text_y <- c(text_y, 2 + 0.5 * v3[2])
        text_lab <- c(text_lab, "-")
        text_col <- c(text_col, "w")
      }

    }

    ## arrow
    xs <- c(1 / 3, 2 / 3, 2 / 3, 5 / 6, 1 / 2, 1 / 6, 1 / 3, 1 / 3) +
          (centerBase - 1) * 1.25
    ys <- c(1 / 4, 1 / 4, 1 / 2, 1 / 2, 3 / 4, 1 / 2, 1 / 2, 1 / 4) + 1
    vs <- rep("arrow", length(xs))

    arrow_poly <- data.frame(x=xs, y=ys, v=vs)
    rect_data <- data.frame(x_start=x_start, x_end=x_end,
                            y_start=y_start, y_end=y_end, rectType=rectType)
    text_data <- data.frame(x=text_x, y=text_y,
                            label=text_lab, text_col=text_col)

    ggplot() +
      geom_rect(data=rect_data, aes(xmin=x_start, xmax=x_end,
                                      ymin=y_start, ymax=y_end,
                                      fill=rectType)) +
      geom_polygon(data=arrow_poly, aes(x=x, y=y, fill=v))  +
      geom_text(data=text_data, aes(label=label, x=x, y=y,
                                      colour=text_col), size=charSize) +
      scale_colour_manual(values=c("#FFFFFF")) +
      scale_fill_manual(
        values=c("A"=baseCol[1], "C"=baseCol[2], "G"=baseCol[3],
                   "T"=baseCol[4], "+"=baseCol[5], "-"=baseCol[6],
                   arrow="#A8A8A8")) +
      guides(fill=FALSE) +
      guides(colour=FALSE) +
      guides(size=FALSE) +
      theme(axis.text=element_blank(),
            axis.ticks=element_blank(),
            panel.background=element_blank(),
            panel.grid=element_blank(),
            axis.title=element_blank(),
            line = element_blank())

}
