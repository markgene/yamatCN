# Pipeline described in the paper Feng, G. A Statistical Method to Estimate DNA
# Copy Number from Illumina High-Density Methylation Arrays. Systems Biomedicine
# (2013).

#' MethylCNV pipeline.
#'
#' The method was described in the paper Feng, G. A Statistical Method to
#' Estimate DNA Copy Number from Illumina High-Density Methylation Arrays.
#' Systems Biomedicine (2013).
#'
#' @param ref Reference samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param qry Query samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param report_dir A character scalar of reporting.
#' @param norm_method A character scalar of method, including raw, illumina,
#'   swan, quantile, noob, funnorm, yamat, dkfz, methylcnv. Default to
#'   "methylcnv".
#' @param batch factor or vector indicating batches. Default to \code{NULL}, -
#'   do not remove batch effect.
#' @param batch2 optional factor or vector indicating a second series of
#'   batches. Default to \code{NULL}, - do not remove batch effect.
#' @param pipe_file A character scalar of the filename storing the pipeline
#'   result sample by sample. Default to "methylCNV.Rda". Turn off the saving by
#'   set it to NULL.
#' @param dnacopy_seed An integer scalar to set the seed. Default to 1.
#' @param overwrite A logical scalar. Default to FALSE.
#' @param verbose A logical scalar. Default to TRUE.
#' @return An object of \code{\link{CMethylCNVPipe}} class.
#' @details The \code{ref} and \code{qry} should include at least the following
#'   columns \code{Sample_Name}, \code{Gender}, \code{Batch}, \code{Sample_Prep}.
#'   Batch effects are removed on methylation and unmethylation signals (log2
#'   transformed) separately.
#' @note The method was named methylCNV according to Feng's paper. It is supposed
#'   to be implemented in the Lumi package. However, I did not find it by
#'   searching the name or CNV in its manual.
#' @export
methylcnv_pipe <-
  function(ref,
           qry,
           report_dir,
           norm_method = c("methylcnv", "swan", "illumina", "raw", "quantile", "noob", "funnorm", "yamat", "dkfz"),
           batch = NULL,
           batch2 = NULL,
           pipe_file = "methylCNV.Rda",
           dnacopy_seed = 1,
           overwrite = FALSE,
           verbose = TRUE) {
    # Check arguments.
    .check_args_pipe(report_dir = report_dir)
    # Preprocess.
    if (verbose) {
      message("Preprocessing...")
      tictoc::tic()
    }
    x <-
      preprocess(
        ref = ref,
        qry = qry,
        norm_method = norm_method,
        batch = batch,
        batch2 = batch2,
        overwrite = overwrite,
        verbose = verbose
      )
    if (verbose) tictoc::toc()
    # Calculate LRR
    if (verbose) {
      message("Calculating LRR...")
      tictoc::tic()
    }
    ref_idx <- .reference_indices(x, label = "ref")
    qry_idx <- .query_indices(x, label = "query")
    lrr <-
      cal_log2_ratio(
        query = minfi::getCN(x[, qry_idx]),
        reference = minfi::getCN(x[, ref_idx]),
        method = "mean"
      )
    if (verbose) tictoc::toc()
    # DNAcopy
    if (verbose) {
      message("Segmentation with DNAcopy...")
      tictoc::tic()
    }
    dnacopy_obj <-
      dnacopy_analysis.default(
        x = x[, qry_idx],
        lrr = lrr,
        seed = dnacopy_seed,
        verbose = verbose
      )
    dnacopy_segment_summary <- summarize_dnacopy_segments(dnacopy_obj$DNAcopy)
    if (verbose) tictoc::toc()
    # Create MethylCNVPipe object
    dat <- new.env(parent = parent.frame())
    dat$dnacopy <- dnacopy_obj
    dat$minfi <- x
    dat$segments <- dnacopy_segment_summary
    dat$lrr <- lrr
    pipe_obj <- MethylCNVPipe(dat = dat)
    # Save the result
    if (is.null(pipe_file)) {
      if (verbose)
        message("Skip saving the pipeline.")
    } else {
      save_pipe(
        pipe_obj,
        outdir = report_dir,
        filename = pipe_file,
        overwrite = overwrite,
        verbose = verbose
      )
    }
    pipe_obj
  }
