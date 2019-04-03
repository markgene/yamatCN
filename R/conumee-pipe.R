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
  conumee::CNV.segment(cnv_anl, ...) -> cnv_anl
  cnv_anl@seg$summary %>%
    dplyr::mutate(seg.mean.shifted = seg.mean - cnv_anl@bin$shift) %>%
    dplyr::mutate(cn.shifted = .to_absolute(seg.mean.shifted)) %>%
    .add_cytoband_to_df(.) %>%
    .add_detail_regions_to_df(.) -> segment_df
  cnv_anl@seg$summary <- segment_df
  cnv_anl
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


#' Convert a summary \code{data.frame} to IGV CBS format.
#' @noRd
.to_igv_segment <- function(summary_df, sample_id) {
  summary_df %>%
    dplyr::mutate(Sample = sample_id) %>%
    dplyr::mutate(Chromosome = stringr::str_remove_all(chrom, "^chr")) %>%
    dplyr::rename(Start = loc.start, End = loc.end, Num_Probes = num.mark, Segment_Mean = seg.mean) %>%
    dplyr::select(Sample, Chromosome, Start, End, Num_Probes, Segment_Mean)
}


#' Convert to Shiny segments.
#' @noRd
.to_shiny_segment <- function(summary_df, sample_id, gender = c("Female", "Male")) {
  gender <- match.arg(gender)
  summary_df %>%
    dplyr::arrange(desc(abs(seg.mean.shifted))) %>%
    dplyr::mutate(
      index = seq_len(nrow(summary_df)),
      caseID = sample_id,
      controlID = sample_id,
      chr = stringr::str_remove_all(chrom, "^chr"),
      CN = round(cn.shifted, digits = 3),
      gender = gender
    ) %>%
    dplyr::rename(
      start = loc.start,
      end = loc.end
    ) %>%
    dplyr::select(index, caseID, controlID, chr, start, end, CN, gender, cytoband, detail_region)
}


#' Add cytoband to segment summary \code{data.frame}.
#' @noRd
.add_cytoband_to_df <- function(df) {
  cytobands_file <- system.file("extdata", "cytoBands.txt", package = "yamatCN")
  cytobands_df <- read.table(cytobands_file, stringsAsFactors = FALSE)
  colnames(cytobands_df) <- c("chrom", "start", "end", "name", "gieStain")
  cytobands_gr <- GenomicRanges::makeGRangesFromDataFrame(cytobands_df, keep.extra.columns = TRUE)
  gr <- GenomicRanges::makeGRangesFromDataFrame(df, start.field = "loc.start", end.field = "loc.end")
  ovlp_df <- GenomicRanges::findOverlaps(gr, cytobands_gr) %>%
    as.data.frame()
  ovlp_df$cytoband <- mcols(cytobands_gr)$name[ovlp_df$subjectHits]
  ovlp_summary <- ovlp_df %>%
    dplyr::group_by(queryHits) %>%
    dplyr::summarise(first_band = dplyr::first(cytoband),
                     last_band  = dplyr::last(cytoband)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cytoband = ifelse(
      first_band == last_band,
      first_band,
      paste(first_band, last_band, sep = "-")
    )) %>%
    dplyr::select(queryHits, cytoband)
  df$cytoband <- ""
  df$cytoband[ovlp_summary$queryHits] <- ovlp_summary$cytoband
  df
}

#' Add detailed regions to segment summary \code{data.frame}.
#' @noRd
.add_detail_regions_to_df <- function(df) {
  dr_gr <- detail_regions()
  gr <- GenomicRanges::makeGRangesFromDataFrame(df, start.field = "loc.start", end.field = "loc.end")
  ovlp_df <- GenomicRanges::findOverlaps(gr, dr_gr) %>%
    as.data.frame()
  ovlp_df$detail_region <- mcols(dr_gr)$name[ovlp_df$subjectHits]
  ovlp_summary <- ovlp_df %>%
    dplyr::group_by(queryHits) %>%
    dplyr::summarise(dr = paste(detail_region, sep = ", ", collapse = ", ")) %>%
    dplyr::ungroup() %>%
    dplyr::select(queryHits, dr)
  df$detail_region <- ""
  df$detail_region[ovlp_summary$queryHits] <- ovlp_summary$dr
  df
}
