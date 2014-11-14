#' Earth Mover's Distance algorithm for differential analysis of genomics data.
#'
#' \code{\link{calculate_emd}} will usually be the only function needed.
#'
#'
#' @import emdist
#' @import BiocParallel
#' @import matrixStats
#' @import ggplot2
#' @name emdomics-package
#' @docType package
NULL


#' @export
#' @title Earth Mover's Distance for differential analysis of genomics data
#' @description This is the main user interface to the \pkg{EMDomics} package, and
#' will usually the only function needed.
#'
#' The algorithm is used to compare genomics data between two groups, refered to
#' herein as "group A" and "group B". Usually the data will be gene expression
#' values from array-based or sequence-based experiments, but data from other
#' types of experiments can also be analyzed (i.e. copy number variation).
#'
#' Traditional methods like Significance Analysis of Microarrays (SAM) and Linear
#' Models for Microarray Data (LIMMA) use significance tests based on summary
#' statistics (mean and standard deviation) of the two distributions. This
#' approach tends to give non-significant results if the two distributions are
#' highly heterogeneous, which can be the case in many biological circumstances
#' (e.g sensitive vs. resistant tumor samples).
#'
#' The Earth Mover's Distance algorithm instead computes the "work" needed
#' to transform one distribution into the other, thus capturing possibly
#' valuable information relating to the overall difference in shape between
#' two heterogeneous distributions.
#'
#' The EMD-based algorithm implemented in \pkg{EMDomics} has two main steps.
#' First, a matrix (e.g. of expression data) is divided into data for "group A"
#' and "group B", and the EMD score is calculated using the two groups for each
#' gene in the data set. Next, the labels for group A and group B are randomly
#' permuted a specified number of times, and an EMD score for each permutation is
#' calculated. The median of the permuted scores for each gene is used as
#' the null distribution, and the False Discovery Rate (FDR) is computed for
#' a range of permissive to restrictive significance thresholds. The threshold
#' that minimizes the FDR is defined as the q-value, and is used to interpret
#' the significance of the EMD score analogously to a p-value (e.g. q-value
#' < 0.05 = significant.)
#'
#' Note that q-values of 0 are adjusted to 1/(nperm+1). For this reason, the
#' \code{nperm} parameter should not be too low (the default of 100 is
#' reasonable).
#'
#' @param data A matrix containing genomics data (e.g. gene expression levels).
#' The rownames should contain gene identifiers, while the column names should
#' contain sample identifiers.
#' @param samplesA A vector of sample names identifying samples in \code{data}
#' that belong to "group A". The names must corresponding to column names
#' in \code{data}.
#' @param samplesB A vector of sample names identifying samples in \code{data}
#' that belong to "group B". The names must corresponding to column names
#' in \code{data}.
#' @param binSize The bin size to be used when generating histograms of
#' the data for "group A" and "group B". Defaults to 0.2.
#' @param nperm An integer specifying the number of randomly permuted EMD
#' scores to be computed. Defaults to 100.
#' @param verbose Boolean specifying whether to display progress messages.
#' @return The function returns an \code{\link{EMDomics}} object.
#' @examples
#' # 100 genes, 100 samples
#' dat <- matrix(rnorm(10000), nrow=100, ncol=100)
#' rownames(dat) <- paste("gene", 1:100, sep="")
#' colnames(dat) <- paste("sample", 1:100, sep="")
#'
#' # "group A" = first 50, "group B" = second 50
#' groupA <- colnames(dat)[1:50]
#' groupB <- colnames(dat)[51:100]

#' results <- calculate_emd(dat, groupA, groupB, nperm=10)
#' head(results$emd)
#' @seealso \code{\link{EMDomics}} \code{\link[emdist]{emd2d}}
calculate_emd <- function(data, samplesA, samplesB, binSize=0.2,
                            nperm=100, verbose=TRUE) {

  # transpose and coerce to df (for bplapply)
  data.df <- as.data.frame(t(data))
  sample_names <- rownames(data.df)

  idxA <- match(samplesA, sample_names)
  idxB <- match(samplesB, sample_names)

  # computes EMD score for a single gene
  emd_gene <- function(geneData, idxA, idxB) {

    dataA <- geneData[idxA]
    dataB <- geneData[idxB]

    bins <- seq(floor(min(c(dataA, dataB))),
                ceiling(max(c(dataA, dataB))),
                by=binSize )

    histA <- hist(dataA, breaks=bins, plot=FALSE)
    histB <- hist(dataB, breaks=bins, plot=FALSE)

    densA <- as.matrix(histA$density)
    densB <- as.matrix(histB$density)

    emdist::emd2d(densA, densB)

  }

  # computes log2 fold change
  fc <- function(geneData, idxA, idxB) {

    dataA <- geneData[idxA]
    dataB <- geneData[idxB]

    meanA <- mean(dataA)
    meanB <- mean(dataB)

    log2(2^meanA/2^meanB)
  }

  # ---------- emd ------------

  # calculate emd for each gene
  if (verbose)
    message("Calculating emd...", appendLF=FALSE)

  emd <- unlist(BiocParallel::bplapply(data.df, emd_gene, idxA, idxB))

  emd <- as.matrix(emd)
  colnames(emd) <- "emd"

  if (verbose)
    message("done.")

  if (verbose)
    message("Calculating fold change...", appendLF=FALSE)

  fc <- unlist(BiocParallel::bplapply(data.df, fc, idxA, idxB))

  fc <- as.matrix(fc)
  colnames(fc) <- "fc"

  if (verbose)
    message("done.")


  # --------------- permuted emd scores ----------------

  sample_count <- length(samplesA)+length(samplesB)

  # matrix to hold permuted emd values
  emd.perm <- matrix(nrow=ncol(data.df), ncol=nperm)
  rownames(emd.perm) <- colnames(data.df)
  colnames(emd.perm) <- as.character(1:nperm)

  for (i in 1:nperm) {

    msg <- paste("Calculating permuted emd #", i, " of ",
                 nperm, "...", sep="")

    if (verbose)
      message(msg, appendLF=FALSE)

    # permute samples
    idx.perm <- sample(1:sample_count, replace=FALSE)
    data.perm <- data.df[idx.perm, ]

    # calculate emd for permuted samples
    emd.perm[, i] <- unlist(BiocParallel::bplapply(data.perm, emd_gene,
                                                   idxA, idxB))

    if (verbose)
      message("done.")

  }

  # ------------------ q-values --------------------

  if (verbose)
    message("Calculating q-values...", appendLF=FALSE)

  perm.medians <- matrixStats::rowMedians(emd.perm)

  # generate thresholds and qval matrix
  thr_upper <- ceiling(max(emd))
  thr <- seq(thr_upper, 0, by = -0.001)
  qvals <- matrix(1, nrow=nrow(emd), ncol=length(thr))

  colnames(qvals) <- thr
  rownames(qvals) <- rownames(emd)

  # calculate fdr at each threshold
  j <- 0
  for (d in thr) {

    j <- j+1

    # calculate true discoveries at this threshold
    idx <- which(emd > d)
    n.signif <- length(idx)
    genes.signif <- rownames(emd)[idx]

    # calculate false discoveries at this threshold
    idx <- which(perm.medians > d)
    n.fd <- length(idx)

    fdr <- n.fd/n.signif
    qvals[genes.signif, j] <- fdr

  }

  # final q-value = smallest fdr
  emd.qval <- apply(qvals, 1, min)

  emd.qval <- as.matrix(emd.qval)
  colnames(emd.qval) <- "q-value"

  # adjust q-values of 0 to 1/(nperm+1)
  emd.qval[emd.qval == 0] <- 1/(nperm+1)

  if (verbose)
    message("done.")

  emd <- cbind(emd, fc, emd.qval)

  EMDomics(data, samplesA, samplesB, emd, emd.perm)

}


#' @export
#' @title Calculate EMD score for a single gene
#' @details The data in \code{vec} is divided into "group A" and "group B" by the
#' identifiers given in \code{samplesA} and \code{samplesB}. The \code{\link{hist}}
#' function is used to generate histograms for the two resulting groups, and the
#' densities are retrieved and passed to \code{\link[emdist]{emd2d}} to compute the
#' EMD score.
#' @param vec A named vector containing data (e.g. expression data) for a single
#' gene.
#' @param samplesA A vector of sample names identifying samples in \code{vec}
#' that belong to "group A".
#' @param samplesB A vector of sample names identifying samples in \code{vec}
#' that belong to "group B".
#' @param binSize The bin size to be used when generating histograms for
#' "group A" and "group B".
#' @return The emd score is returned.
#' @examples
#' # 100 samples
#' vec <- rnorm(100)
#' names(vec) <- paste("sample", 1:100, sep="")
#'
#' # "group A" = first 50, "group B" = second 50
#' groupA <- names(vec)[1:50]
#' groupB <- names(vec)[51:100]
#'
#' calculate_emd_gene(vec, groupA, groupB)
#' @seealso \code{\link[emdist]{emd2d}}
calculate_emd_gene <- function(vec, samplesA, samplesB, binSize=0.2) {

  dataA <- vec[samplesA]
  dataB <- vec[samplesB]

  bins <- seq(floor(min(c(dataA, dataB))),
              ceiling(max(c(dataA, dataB))),
              by=binSize )

  histA <- hist(dataA, breaks=bins, plot=FALSE)
  histB <- hist(dataB, breaks=bins, plot=FALSE)

  densA <- as.matrix(histA$density)
  densB <- as.matrix(histB$density)

  emdist::emd2d(densA, densB)

}


#' @export
#' @title Create an EMDomics object
#' @description This is the constructor for objects of class 'EMDomics'. It
#' is used in \code{\link{calculate_emd}} to construct the return value.
#' @param data A matrix containing genomics data (e.g. gene expression levels).
#' The rownames should contain gene identifiers, while the column names should
#' contain sample identifiers.
#' @param samplesA A vector of sample names identifying samples in \code{data}
#' that belong to "group A". The names must corresponding to column names
#' in \code{data}.
#' @param samplesB A vector of sample names identifying samples in \code{data}
#' that belong to "group B". The names must corresponding to column names
#' in \code{data}.
#' @param emd A matrix containing a row for each gene in \code{data}, and with
#' the following columns:
#' \itemize{
#' \item \code{emd} The calculated emd score.
#' \item \code{fc} The log2 fold change of "group A" samples relative to "group B"
#' samples.
#' \item \code{q-value} The calculated q-value.
#' }
#' The row names should specify the gene identifiers for each row.
#' @param emd.perm A matrix containing a row for each gene in \code{data}, and
#' with a column containing emd scores for each random permutation calculated
#' via \code{\link{calculate_emd}}.
#' @return The function combines it's arguments in a list, which is assigned class
#' 'EMDomics'. The resulting object is returned.
#' @seealso \code{\link{calculate_emd}}
EMDomics <- function(data, samplesA, samplesB, emd, emd.perm) {

  structure(list("data"=data, "samplesA"=samplesA, "samplesB"=samplesB,
                 "emd"=emd, "emd.perm"=emd.perm),
            class = "EMDomics")

}
