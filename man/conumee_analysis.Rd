% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conumee-pipe.R
\name{conumee_analysis}
\alias{conumee_analysis}
\title{Main Conumee CNV analysis}
\usage{
conumee_analysis(query, ref, anno, intercept = TRUE, seed = 1, ...)
}
\arguments{
\item{query}{A \code{\link[conumee]{CNV.data-class}} object of query sample
(single sample).}

\item{ref}{A \code{\link[conumee]{CNV.data-class}} object of reference set.}

\item{anno}{A \code{\link[conumee]{CNV.anno-class}} object. Use
\code{\link[conumee]{CNV.create_anno}} or \code{\link{CNV.create_anno2}}
to create.}

\item{intercept}{A logical scalar setting the \code{intercept} argument of
\code{\link[conumee]{CNV.fit}}. Should intercept be considered? Defaults
to TRUE.}

\item{seed}{An integer scalar to set the seed. Default to 1.}

\item{...}{Any arguments passed to \code{\link[conumee]{CNV.segment}}.}
}
\value{
An object of \code{\link[conumee]{CNV.analysis-class}}.
}
\description{
Main Conumee CNV analysis
}
