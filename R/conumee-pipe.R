# Pipeline Based-upon Conumee Package.

#' CNV pipeline based-upon Conumee package.
#'
#' @param ref Reference samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param qry Query samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param report_dir A character scalar of reporting directory.
#' @param norm_method A character scalar of method, including raw, illumina,
#'   swan, quantile, noob, funnorm, yamat, dkfz, quantile. Default to "swan".
#' @param batch factor or vector indicating batches. Default to \code{NULL}, -
#'   do not remove batch effect.
#' @param batch2 optional factor or vector indicating a second series of
#'   batches. Default to \code{NULL}, - do not remove batch effect.
#' @param param A list for generating sample CNV report. Default to \code{list()}.
#' @param overwrite A logical scalar. Default to FALSE.
#' @param verbose A logical scalar. Default to TRUE.
#' @return TBA.
#' @details The \code{ref} and \code{qry} should include at least the following
#'   columns \code{Sample_Name}, \code{Gender}, \code{Batch}, \code{Sample_Prep}.
#'   Batch effects are removed on methylation and unmethylation signals (log2
#'   transformed) separately.
#' @export
cn_pipe_conumee <-
  function(ref,
           qry,
           report_base,
           norm_method = c("swan", "dkfz", "illumina", "raw", "yamat", "noob", "funnorm", "quantile"),
           batch = NULL,
           batch2 = NULL,
           param = list(),
           overwrite = FALSE,
           verbose = TRUE) {
  }


