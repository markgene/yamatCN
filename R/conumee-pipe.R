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
           report_dir,
           norm_method = c("swan", "illumina", "raw", "quantile", "noob", "funnorm", "yamat", "dkfz"),
           batch = NULL,
           batch2 = NULL,
           param = list(),
           overwrite = FALSE,
           verbose = TRUE) {
    # Check arguments.
    norm_method <- match.arg(norm_method)
    .check_cn_pipe_conumee(
      ref = ref,
      qry = qry,
      report_dir = report_dir,
      norm_method = norm_method,
      batch = batch,
      batch2 = batch2,
      param = param
    )
  }


#' Check arguments of \code{cn_pipe_conumee}.
#'
#' @param ref Reference samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param qry Query samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param report_dir A character scalar of reporting directory.
#' @param norm_method A character scalar of method, including raw, illumina,
#'   swan, quantile, noob, funnorm, yamat, dkfz, quantile.
#' @param batch factor or vector indicating batches.
#' @param batch2 optional factor or vector indicating a second series of
#'   batches.
#' @param param A list for generating sample CNV report.
#' @return TRUE if pass all checks (invisible).
#' @details The \code{ref} and \code{qry} should include at least the following
#'   columns \code{Sample_Name}, \code{Gender}, \code{Batch}, \code{Sample_Prep}.
#'   Batch effects are removed on methylation and unmethylation signals (log2
#'   transformed) separately.
#' @noRd
.check_cn_pipe_conumee <-
  function(ref,
           qry,
           report_dir,
           norm_method,
           batch,
           batch2,
           param) {
    if (missing(ref))
      stop("ref is required!")
    if (missing(qry))
      stop("qry is required!")
    if (missing(report_dir))
      stop("report_dir is required!")
    norm_methods <-
      c("yamat",
        "dkfz",
        "illumina",
        "raw",
        "swan",
        "noob",
        "funnorm",
        "quantile")
    if (!all(norm_method %in% norm_methods))
      stop("Invalid norm_method!")
    # If batch and batch2 are in phenotype information?
    # Sample preparation is coded in the same way?
    qry_df <- minfi::pData(qry)
    ref_df <- minfi::pData(ref)
    common_colnames <- colnames(qry_df)[colnames(qry_df) %in% colnames(ref_df)]
    if (!is.null(batch)) {
      if (!batch %in% common_colnames)
        stop("batch is not present in both query and reference.")
      qry_batch <- unique(minfi::pData(qry)[, batch])
      ref_batch <- unique(minfi::pData(ref)[, batch])
      if (!all(qry_batch %in% ref_batch)) {
        stop("Query has batch which are not found in reference.")
      }
    }
    if (!is.null(batch2)) {
      if (!batch2 %in% common_colnames)
        stop("batch2 is not present in both query and reference.")
      qry_batch2 <- unique(minfi::pData(qry)[, batch2])
      ref_batch2 <- unique(minfi::pData(ref)[, batch2])
      if (!all(qry_batch2 %in% ref_batch2)) {
        stop("Query has batch2 which are not found in reference.")
      }
    }
    invisible(TRUE)
  }

