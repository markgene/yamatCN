# Pipeline Based-upon Conumee Package.

#' CNV pipeline based-upon Conumee package.
#'
#' @param ref Reference samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param qry Query samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param report_dir A character scalar of reporting.
#' @param norm_method A character scalar of method, including raw, illumina,
#'   swan, quantile, noob, funnorm, yamat, dkfz, quantile, methylcnv. Default
#'   to "swan".
#' @param batch factor or vector indicating batches. Default to \code{NULL}, -
#'   do not remove batch effect.
#' @param batch2 optional factor or vector indicating a second series of
#'   batches. Default to \code{NULL}, - do not remove batch effect.
#' @param batch_effect_method A character scalar of methods, including
#'   "mum" (methylation and unmethylation separately), "cn" (copy number),
#'   "beta" (beta-value). Default to "cn".
#' @param pipe_file A character scalar of the filename storing the pipeline
#'   result sample by sample. Default to "conumee.Rda". Turn off the saving by
#'   set it to NULL.
#' @param dnacopy_seed An integer scalar to set the seed. Default to 1.
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
           norm_method = c("swan", "illumina", "raw", "quantile", "noob", "funnorm", "yamat", "dkfz", "methylcnv"),
           batch = NULL,
           batch2 = NULL,
           batch_effect_method = c("cn", "mum", "beta"),
           pipe_file = "conumee.Rda",
           dnacopy_seed = 1,
           overwrite = FALSE,
           verbose = TRUE,
           ...) {
    # Check arguments.
    .check_args_pipe(report_dir = report_dir)
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
    conumee_backbone(
      x = x,
      report_dir = report_dir,
      pipe_file = pipe_file,
      dnacopy_seed = dnacopy_seed,
      overwrite = overwrite,
      verbose = verbose,
      ...
    )
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
      conumee::CNV.load(., names = minfi::sampleNames(x))
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
#' @param seed An integer scalar to set the seed. Default to 1.
#' @param ... Any arguments passed to \code{\link[conumee]{CNV.segment}}.
#' @return An object of \code{\link[conumee]{CNV.analysis-class}}.
#' @export
conumee_analysis <- function(query, ref, anno, intercept = TRUE, seed = 1, ...) {
  conumee::CNV.fit(
    query = query,
    ref = ref,
    anno = anno) %>%
    conumee::CNV.bin() %>%
    conumee::CNV.detail() -> cnv_anl
  set.seed(seed)
  conumee::CNV.segment(cnv_anl, ...)
}


#' CNV pipeline backbone.
#'
#' CNV pipeline backbone is the steps after preprocessing, including loading
#' data into Conumee, analysis and write the result.
#'
#' @param x An object of \code{\link[minfi]{GenomicRatioSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}}.
#' @param report_dir A character scalar of reporting.
#' @param pipe_file A character scalar of the filename storing the pipeline
#'   result sample by sample. Default to "conumee.Rda". Turn off the saving by
#'   set it to NULL.
#' @param dnacopy_seed An integer scalar to set the seed. Default to 1.
#' @param overwrite A logical scalar. Default to FALSE.
#' @param verbose A logical scalar. Default to TRUE.
#' @param ... Any arguments passed to \code{\link[conumee]{CNV.segment}}.
#' @return An object of \code{\link{ConumeePipe}} class.
#' @details The \code{ref} and \code{qry} should include at least the following
#'   columns \code{Sample_Name}, \code{Gender}, \code{Batch}, \code{Sample_Prep}.
#'   Batch effects are removed on methylation and unmethylation signals (log2
#'   transformed) separately.
#' @export
conumee_backbone <- function(x,
                             report_dir,
                             pipe_file = "conumee.Rda",
                             dnacopy_seed = 1,
                             overwrite = FALSE,
                             verbose = TRUE,
                             ...) {
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
  ref_idx <- .reference_indices(x, label = "ref")
  qry_idx <- .query_indices(x, label = "query")
  cnv_results <- lapply(
    qry_idx,
    function(idx) {
      i <- which(qry_idx == idx)
      nm <- minfi::sampleNames(x)[idx]
      if (verbose) {
        message("Run CNV for sample ", "(", i, "/", length(qry_idx), "): ", nm)
        tictoc::tic()
      }
      cnv_res <-
        conumee_analysis(cnv_dat[i], cnv_dat[ref_idx], anno = cnv_anno, seed = dnacopy_seed, ...)
      if (verbose) tictoc::toc()
      cnv_res
    })
  names(cnv_results) <- minfi::sampleNames(x)[qry_idx]
  if (verbose) tictoc::toc()
  # Create ConumeePipe object
  dat <- new.env(parent = parent.frame())
  dat$conumee_results <- cnv_results
  dat$minfi <- x
  cp_obj <- ConumeePipe(dat = dat)
  # Save the result
  if (is.null(pipe_file)) {
    if (verbose)
      message("Skip saving the pipeline result sample by sample.")
  } else {
    save_pipe(
      cp_obj,
      outdir = report_dir,
      filename = pipe_file,
      overwrite = overwrite,
      verbose = verbose
    )
  }
  cp_obj
}

