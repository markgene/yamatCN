# Tune Pipeline Based-upon Conumee Package.

#' Run Conumee pipeline with different preprocessing methods.
#'
#' Run Conumee pipeline with different preprocessing methods, i.e. the
#' combination of normalization methods and batch effect correction. The
#' function serves as a template for writing your own function to run
#' the pipeline with different parameters, rather than a function of
#' general purpose.
#'
#' @param ref Reference samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param qry Query samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param report_dir A character scalar of reporting.
#' @param norm_methods A character vector of methods, including raw, illumina,
#'   swan, quantile, noob, funnorm, yamat, dkfz, quantile, methylcnv. Default
#'   to all of them.
#' @param batch factor or vector indicating batches. Default to \code{NULL}, -
#'   do not remove batch effect.
#' @param batch2 optional factor or vector indicating a second series of
#'   batches. Default to \code{NULL}, - do not remove batch effect.
#' @param pipe_file A character scalar of the filename storing the pipeline
#'   result sample by sample. Turn off the saving by set it to \code{NULL}.
#'   Default to \code{NULL}.
#' @param dnacopy_seed An integer scalar to set the seed. Default to 1.
#' @param overwrite A logical scalar. Default to FALSE.
#' @param verbose A logical scalar. Default to TRUE.
#' @param ... Any arguments passed to \code{\link[conumee]{CNV.segment}}.
#' @return An object of \code{\link{ConumeePipe}} class.
#' @details The \code{ref} and \code{qry} should include at least the following
#'   columns \code{Sample_Name}, \code{Gender}, \code{Batch}, \code{Sample_Prep}.
#'   Batch effects are removed on methylation and unmethylation signals (log2
#'   transformed) separately.
#'   The function outputs the result in a directory. Each subdirectory is a
#'   normalzation method. A deeper subdirectory is named in a pattern of
#'   "\code{batch}-\code{batch2}". Then, it contains the subdirectories of
#'   sample names, in which there is a segment table named "segments.csv".
#' @export
tune_conumee_pipe <-
  function(ref,
           qry,
           report_dir,
           norm_methods = c("swan", "illumina", "raw", "quantile", "noob", "funnorm", "yamat", "dkfz", "methylcnv"),
           batch = NULL,
           batch2 = NULL,
           pipe_file = NULL,
           dnacopy_seed = 1,
           overwrite = FALSE,
           verbose = TRUE,
           ...) {
    # Check arguments.
    .check_args_pipe(report_dir = report_dir)
    # Batch grid
    batch_grid <- .batch_grid(batch, batch2)
    # Segment grid
    segment_grid <- .conumee_segment_grid()
    # Search grid
    lapply(
      norm_methods,
      function(m) {
        by(batch_grid, seq(nrow(batch_grid)), function(batch_row) {
          b1 <- batch_row$batch
          b2 <- batch_row$batch2
          if (verbose) {
            message("Run method=", m, ", batch=", b1, ", batch2=", b2)
          }
          if (is.na(b1)) b1 <- NULL
          if (is.na(b2)) b2 <- NULL
          # Preprocess.
          x <-
            preprocess(
              ref = ref,
              qry = qry,
              norm_method = m,
              batch = b1,
              batch2 = b2,
              overwrite = overwrite,
              verbose = verbose
            )
          by(segment_grid, seq(nrow(segment_grid)), function(segment_args) {
            outdir <- file.path(report_dir, m, batch_row$name, segment_args$name)
            # Run Conumee-based CNV analysis
            conumee_backbone(
              x = x,
              report_dir = outdir,
              pipe_file = pipe_file,
              dnacopy_seed = dnacopy_seed,
              overwrite = overwrite,
              verbose = verbose,
              alpha = segment_args$alpha,
              min.width = segment_args$min.width,
              undo.SD = segment_args$undo.SD,
              undo.splits = "sdundo",
              ...
            ) %>%
              write_segments(., outdir)
          })
        })
      }) -> tmp
  }


.batch_grid <- function(batch, batch2) {
  if (is.null(batch)) {
    batch <- NA
  } else {
    batch <- c(batch, NA)
  }
  if (is.null(batch2)) {
    batch2 <- NA
  } else {
    batch2 <- c(batch2, NA)
  }
  expand.grid(batch = batch, batch2 = batch2) %>%
    dplyr::filter(!(is.na(batch) & !is.na(batch2))) %>%
    dplyr::mutate(name = paste(batch, batch2, sep = "-"))
}


.conumee_segment_grid <-
  function(alpha = c(0.001, 0.005, 0.01),
           min.width = c(2, 3, 4, 5),
           undo.SD = seq(from = 2, to = 3, by = 0.1)) {
    expand.grid(
      alpha = alpha,
      min.width = min.width,
      undo.SD = undo.SD
    ) %>%
      dplyr::mutate(name = paste0("alpha", alpha, "-minw", min.width, "-SD", undo.SD))
  }
