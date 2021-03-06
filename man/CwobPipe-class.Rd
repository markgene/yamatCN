% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipe-classes.R
\docType{class}
\name{CwobPipe-class}
\alias{CwobPipe-class}
\alias{CwobPipe}
\title{CwobPipe pipeline}
\description{
CwobPipe pipeline
}
\section{Slots}{

\describe{
\item{\code{dat}}{An environment which has the following objects:
\itemize{
  \item \code{minfi}: an object of
    \code{\link[minfi]{GenomicRatioSet-class}} or
    \code{\link[minfi]{GenomicMethylSet-class}}.
  \item \code{dnacopy} an list of two items, i) \code{DNAcopy} an object of
    \code{\link[DNAcopy]{DNAcopy}} class, and ii) \code{CNA} an object of
    \code{\link[DNAcopy]{CNA}} class
  \item \code{segments} a \code{data.frame} of segments.
  \item \code{lrr} A \code{matrix} of LRR.
  \item \code{lrr_shift} a numeric scalar. LRR is shifted to minimize the
    median absolute deviation from all bins to zero to determine the
    copy-number neutral state.
}}

\item{\code{meta}}{A list contain meta information.}
}}

