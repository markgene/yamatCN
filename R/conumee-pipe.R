# Pipeline Based-upon Conumee Package.

#' CNV pipeline based-upon Conumee package.
#'
#' @param ref Reference samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param qry Query samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param report_dir A character scalar of reporting.
#' @param norm_method A character scalar of method, including raw, illumina,
#'   swan, quantile, noob, funnorm, yamat, dkfz, quantile. Default to "swan".
#' @param batch factor or vector indicating batches. Default to \code{NULL}, -
#'   do not remove batch effect.
#' @param batch2 optional factor or vector indicating a second series of
#'   batches. Default to \code{NULL}, - do not remove batch effect.
#' @param pipe_file A character scalar of the filename storing the pipeline
#'   result sample by sample. Default to "conumee.Rda". Turn off the saving by
#'   set it to NULL.
#' @param overwrite A logical scalar. Default to FALSE.
#' @param verbose A logical scalar. Default to TRUE.
#' @param ... Any arguments passed to \code{\link[conumee]{CNV.segment}}.
#' @return An object of \code{\link{ConumeePipe}} class.
#' @details The \code{ref} and \code{qry} should include at least the following
#'   columns \code{Sample_Name}, \code{Gender}, \code{Batch}, \code{Sample_Prep}.
#'   Batch effects are removed on methylation and unmethylation signals (log2
#'   transformed) separately.
#' @export
conumee_pipe <-
  function(ref,
           qry,
           report_dir,
           norm_method = c("swan", "illumina", "raw", "quantile", "noob", "funnorm", "yamat", "dkfz"),
           batch = NULL,
           batch2 = NULL,
           pipe_file = "conumee.Rda",
           overwrite = FALSE,
           verbose = TRUE) {
    # Check arguments.
    .check_args_conumee_pipe(report_dir = report_dir)
    # Preprocess.
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
    # Run Conumee-based CNV analysis
    # Create Conumee annotation.
    if (verbose) {
      message("Preparing CNV analysis: annotation...")
      tictoc::tic()
    }
    cnv_anno <- CNV.create_anno.yamat(x, detail_regions_preset = "conumee")
    if (verbose) tictoc::toc()
    # Load data to Conumee.
    if (verbose) {
      message("Preparing CNV analysis: loading data...")
      tictoc::tic()
    }
    cnv_dat <- .CNV_load(x)
    if (verbose) tictoc::toc()
    # Main Conumee CNV analysis.
    if (verbose) {
      message("Running CNV analysis: main analysis...")
      tictoc::tic()
    }
    ref_samples <- .reference_sample_names(x, label = "ref")
    qry_samples <- .query_sample_names(x, label = "query")
    cnv_results <- lapply(
      qry_samples,
      function(nm) {
        i <- which(qry_samples == nm)
        if (verbose) {
          message("Run CNV for sample ", "(", i, "/", length(qry_samples), "): ", nm)
          tictoc::tic()
        }
        cnv_res <-
          conumee_analysis(cnv_dat[nm], cnv_dat[ref_samples], anno = cnv_anno)
        if (verbose) tictoc::toc()
        cnv_res
      })
    if (verbose) tictoc::toc()
    # Create ConumeePipe object
    dat <- new.env(parent = parent.frame())
    dat$conumee_results <- cnv_results
    dat$minfi_obj <- x
    dat$batch <- batch
    dat$batch2 <- batch2
    cp_obj <- ConumeePipe(dat = dat)
    # Save the result
    if (is.null(pipe_file)) {
      if (verbose)
        message("Skip saving the pipeline result sample by sample.")
    } else {
      save_pipe(cp_obj, outdir = report_dir, filename = pipe_file)
    }
    cp_obj
  }


#' Check arguments of \code{cn_pipe_conumee}.
#'
#' @param report_dir A character scalar of reporting.
#' @return TRUE if pass all checks (invisible).
#' @noRd
.check_args_conumee_pipe <-
  function(report_dir) {
    if (missing(report_dir))
      stop("report_dir is required!")
    invisible(TRUE)
  }


#' A wrapper of \code{\link[conumee]{CNV.load}()}.
#'
#' @param x An object of \code{\link[minfi]{GenomicRatioSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}}.
#' @return An object of \code{\link[conumee]{CNV.data-class}} object.
#' @noRd
.CNV_load <- function(x) {
  if (inherits(x, "GenomicRatioSet")) {
    cnv_dat <- conumee::CNV.load(x)
  } else if (inherits(x, "GenomicMethylSet")) {
    cnv_dat <-
      minfi::ratioConvert(x) %>%
      conumee::CNV.load()
  } else {
    stop(
      "Invalid class. Require GenomicRatioSet or GenomicMethylSet, but ",
      as.character(class(RGsetEx), " is provided.")
    )
  }
  cnv_dat
}


#' Main Conumee CNV analysis
#'
#' @param query A \code{\link[conumee]{CNV.data-class}} object of query sample
#'   (single sample).
#' @param ref A \code{\link[conumee]{CNV.data-class}} object of reference set.
#' @param anno A \code{\link[conumee]{CNV.anno-class}} object. Use
#'   \code{\link[conumee]{CNV.create_anno}} or \code{\link{CNV.create_anno2}}
#'   to create.
#' @param intercept A logical scalar setting the \code{intercept} argument of
#'    \code{\link[conumee]{CNV.fit}}. Should intercept be considered? Defaults
#'    to TRUE.
#' @param ... Any arguments passed to \code{\link[conumee]{CNV.segment}}.
#' @return An object of \code{\link[conumee]{CNV.analysis-class}}.
#' @export
conumee_analysis <- function(query, ref, anno, intercept = TRUE, ...) {
  conumee::CNV.fit(
    query = query,
    ref = ref,
    anno = anno) %>%
    conumee::CNV.bin() %>%
    conumee::CNV.detail() %>%
    conumee::CNV.segment(...)
}



#' Query sample names.
#'
#' @param x An object of \code{\link[minfi]{RGChannelSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}}.
#' @param label A character vector labeling query samples
#' @return A character vector of query sample names.
#' @noRd
.query_sample_names <- function(x, label = "query") {
  query_sample_names <-
    minfi::pData(x) %>%
    as.data.frame() %>%
    dplyr::filter(ref_query %in% label) %>%
    dplyr::select(Sample_Name) %>%
    unlist()
  names(query_sample_names) <- query_sample_names
  query_sample_names
}


#' Reference sample names.
#'
#' @param x An object of \code{\link[minfi]{RGChannelSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}}.
#' @param label A character vector labeling reference samples
#' @return A character vector of reference sample names.
#' @noRd
.reference_sample_names <- function(x, label = "ref") {
  ref_sample_names <-
    minfi::pData(x) %>%
    as.data.frame() %>%
    dplyr::filter(ref_query %in% label) %>%
    dplyr::select(Sample_Name) %>%
    unlist()
  names(ref_sample_names) <- ref_sample_names
  ref_sample_names
}
