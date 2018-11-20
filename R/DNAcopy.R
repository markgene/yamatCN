# Wrapper around DNAcopy


#' Run DNAcopy analysis.
#'
#' @param x A \code{minfi} object such as
#'   \code{\link[minfi]{RGChannelSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}} or
#'   \code{\link[minfi]{GenomicRatioSet-class}}.
#' @param lrr A matrix of log2 transformed ratio of query and reference copy
#'   number (sum of methylation and unmethylation signals).
#' @param seed An integer scalar to set the seed. Default to 1.
#' @param verbose A logical scalar. Default to TRUE.
#' @param ... Arugments passed to \code{\link[DNAcopy]{segment}}.
#' @return A list of two element:
#'   \itemize{
#'     \item \code{DNAcopy} a \code{\link[DNAcopy]{DNAcopy}} object.
#'     \item \code{CNA} a code{\link[DNAcopy]{CNA}} object (not smoothed).
#'   }
#' @details The function provides a template of run DNAcopy analysis. There are
#'   many parameters in smoothing and segmentation. The function uses the
#'   default settings for smoothing.
#' @export
dnacopy_analysis <- function(x, lrr, seed = 1, verbose = TRUE, ...) {
  output <- list()
  # Creating DNAcopy CNA object.
  if (verbose) {
    message("Creating DNAcopy CNA object...")
    tictoc::tic()
  }
  cna_obj <- create_dnacopy_cna(x, lrr)
  output$CNA <- cna_obj
  if (verbose) tictoc::toc()
  # Smooth
  if (verbose) {
    message("Smoothing...")
    tictoc::tic()
  }
  cna_obj <-
    DNAcopy::smooth.CNA(
      cna_obj,
      smooth.region = 10,
      outlier.SD.scale = 4,
      smooth.SD.scale = 2,
      trim = 0.025
    )
  if (verbose) tictoc::toc()
  # Segmentation
  if (verbose) {
    message("Segmentation...")
    tictoc::tic()
  }
  segment_verbose <- 1
  set.seed(seed)
  output$DNAcopy <- DNAcopy::segment(cna_obj, verbose = segment_verbose, ...)
  if (verbose) tictoc::toc()
  output
}


#' Run DNAcopy analysis: DNAcopy's default settings.
#'
#' Use the default settings of \code{\link[DNAcopy]{segment}}.
#'
#' @param x A \code{minfi} object such as
#'   \code{\link[minfi]{RGChannelSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}} or
#'   \code{\link[minfi]{GenomicRatioSet-class}}.
#' @param lrr A matrix of log2 transformed ratio of query and reference copy
#'   number (sum of methylation and unmethylation signals).
#' @param seed An integer scalar to set the seed. Default to 1.
#' @param verbose A logical scalar. Default to TRUE.
#' @return A \code{\link[DNAcopy]{DNAcopy}} object.
#' @export
dnacopy_analysis.default <- function(x, lrr, seed = 1, verbose = TRUE) {
  dnacopy_analysis(
    x,
    lrr,
    seed = seed,
    verbose = verbose,
    # weights = NULL,
    alpha = 0.01,
    nperm = 10000,
    p.method = "hybrid",
    min.width = 2,
    kmax = 25,
    nmin = 200,
    eta = 0.05,
    sbdry = NULL,
    trim = 0.025,
    undo.splits = "none",
    undo.prune = 0.05,
    undo.SD = 3
  )
}


#' Run DNAcopy analysis: Conumee's settings.
#'
#' Use Conumee's settings of \code{\link[DNAcopy]{segment}}. Notice: Conumee
#' uses the bins of probes as markers rather than probes themselves.
#'
#' @param x A \code{minfi} object such as
#'   \code{\link[minfi]{RGChannelSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}} or
#'   \code{\link[minfi]{GenomicRatioSet-class}}.
#' @param lrr A matrix of log2 transformed ratio of query and reference copy
#'   number (sum of methylation and unmethylation signals).
#' @param seed An integer scalar to set the seed. Default to 1.
#' @param verbose A logical scalar. Default to TRUE.
#' @return A \code{\link[DNAcopy]{DNAcopy}} object.
#' @export
dnacopy_analysis.conumee <- function(x, lrr, seed = 1, verbose = TRUE) {
  dnacopy_analysis(
    x,
    lrr,
    seed = seed,
    verbose = verbose,
    # weights = NULL,
    alpha = 0.001,
    nperm = 50000,
    p.method = "hybrid",
    min.width = 5,
    kmax = 25,
    nmin = 200,
    eta = 0.05,
    sbdry = NULL,
    trim = 0.025,
    undo.splits = "sdundo",
    undo.prune = 0.05,
    undo.SD = 2.2
  )
}


#' Run DNAcopy analysis: Yamat's settings.
#'
#' Use Yamat's settings of \code{\link[DNAcopy]{segment}}.
#'
#' @param x A \code{minfi} object such as
#'   \code{\link[minfi]{RGChannelSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}} or
#'   \code{\link[minfi]{GenomicRatioSet-class}}.
#' @param lrr A matrix of log2 transformed ratio of query and reference copy
#'   number (sum of methylation and unmethylation signals).
#' @param seed An integer scalar to set the seed. Default to 1.
#' @param verbose A logical scalar. Default to TRUE.
#' @return A \code{\link[DNAcopy]{DNAcopy}} object.
#' @export
dnacopy_analysis.yamat <- function(x, lrr, seed = 1, verbose = TRUE) {
  dnacopy_analysis(
    x,
    lrr,
    seed = seed,
    verbose = verbose,
    # weights = NULL,
    alpha = 0.001,
    nperm = 50000,
    p.method = "hybrid",
    min.width = 5,
    kmax = 25,
    nmin = 200,
    eta = 0.05,
    sbdry = NULL,
    trim = 0.025,
    undo.splits = "sdundo",
    undo.prune = 0.05,
    undo.SD = 3
  )
}


#' Create \code{\link[DNAcopy]{CNA}} object.
#'
#' @param x A \code{minfi} object such as
#'   \code{\link[minfi]{RGChannelSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}} or
#'   \code{\link[minfi]{GenomicRatioSet-class}}.
#' @param lrr A matrix of log2 transformed ratio of query and reference copy
#'   number (sum of methylation and unmethylation signals).
#' @param verbose A logical scalar. Default to TRUE.
#' @return A \code{\link[DNAcopy]{CNA}} object.
#' @export
create_dnacopy_cna <- function(x, lrr, verbose = TRUE) {
  if (missing(x))
    stop("Require argument x")
  if (missing(lrr))
    stop("Require argument lrr")
  # Get loci information
  if (verbose) {
    message("Getting loci information...")
    tictoc::tic()
  }
  anno <- minfi::getAnnotation(x, lociNames = rownames(lrr))
  if (verbose) tictoc::toc()
  # Create CNA object.
  if (verbose) {
    message("Creating CNA object...")
    tictoc::tic()
  }
  cna_obj <- DNAcopy::CNA(
    genomdat = lrr[rownames(anno),],
    chrom = anno$chr,
    maploc = anno$pos,
    data.type = "logratio",
    sampleid = colnames(lrr),
    presorted = FALSE
  )
  if (verbose) tictoc::toc()
  cna_obj
}


#' Summarize \code{DNAcopy} segments.
#'
#' @param x A \code{\link[DNAcopy]{DNAcopy}} object.
#' @param ... Arugments passed to \code{\link[DNAcopy]{segment.p}}.
#' @return A \code{data.frame}.
#' @export
summarize_dnacopy_segments <- function(x, ...) {
  summary_df <- DNAcopy::segments.summary(x)
  pval_df <- DNAcopy::segments.p(x, ...)
  cbind(summary_df, pval_df[, c("bstat", "pval")])
}


