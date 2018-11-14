# Preprocessing (normalization and batch effect correction) for Copy-number
# Analysis


#' Preprocess
#'
#' Preprocessing include normalization and batch effect removal.
#'
#' @param ref Reference samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param qry Query samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param norm_method A character scalar of method, including raw, illumina,
#'   swan, quantile, noob, funnorm, yamat, dkfz, quantile, methylcnv. Default
#'   to "swan".
#' @param batch factor or vector indicating batches. Default to \code{NULL}, -
#'   do not remove batch effect.
#' @param batch2 optional factor or vector indicating a second series of
#'   batches. Default to \code{NULL}, - do not remove batch effect.
#' @param overwrite A logical scalar. Default to FALSE.
#' @param verbose A logical scalar. Default to TRUE.
#' @return An object of \code{\link[minfi]{GenomicMethylSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}}. If \code{norm_method} is quantile,
#'   an object of \code{\link[minfi]{GenomicRatioSet-class}} is returned.
#' @details The \code{ref} and \code{qry} should include at least the following
#'   columns \code{Sample_Name}, \code{Gender}, \code{Batch}, \code{Sample_Prep}.
#'   Batch effects are removed on methylation and unmethylation signals (log2
#'   transformed) separately.
#' @export
preprocess <-
  function(ref,
           qry,
           norm_method = c("swan", "illumina", "raw", "quantile", "noob", "funnorm", "yamat", "dkfz", "methylcnv"),
           batch = NULL,
           batch2 = NULL,
           overwrite = FALSE,
           verbose = TRUE) {
    # Check arguments.
    norm_method <- match.arg(norm_method)
    .check_args_preprocess(
      ref = ref,
      qry = qry,
      norm_method = norm_method,
      batch = batch,
      batch2 = batch2
    )
    # Combine and normalize
    if (verbose) {
      message("Combining ref and query and then normalizing...")
      tictoc::tic()
    }
    x <- .combine(ref = ref, qry = qry, batch = batch, batch2 = batch2) %>%
      yamat::normalize(rgset = ., norm_method = norm_method)
    if (verbose) tictoc::toc()
    # Remove batch effects
    if (!is.null(batch)) {
      if (verbose) {
        message("Removing batch effects...")
        tictoc::tic()
      }
      x <- yamat::remove_batch_effect(x, batch = batch, batch2 = batch2)
      if (verbose) tictoc::toc()
    }
    x
  }


#' Check arguments of \code{preprocess}.
#'
#' @param ref Reference samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param qry Query samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param norm_method A character scalar of method, including raw, illumina,
#'   swan, quantile, noob, funnorm, yamat, dkfz, quantile, methylcnv.
#' @param batch factor or vector indicating batches.
#' @param batch2 optional factor or vector indicating a second series of
#'   batches.
#' @return TRUE if pass all checks (invisible).
#' @details The \code{ref} and \code{qry} should include at least the following
#'   columns \code{Sample_Name}, \code{Gender}, \code{Batch}, \code{Sample_Prep}.
#'   Batch effects are removed on methylation and unmethylation signals (log2
#'   transformed) separately.
#' @noRd
.check_args_preprocess <-
  function(ref,
           qry,
           norm_method,
           batch,
           batch2) {
    if (missing(ref))
      stop("ref is required!")
    if (missing(qry))
      stop("qry is required!")
    norm_methods <-
      c("yamat",
        "dkfz",
        "illumina",
        "raw",
        "swan",
        "noob",
        "funnorm",
        "quantile",
        "methylcnv")
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


#' Combine reference and query samples.
#'
#' Process the phenotype information and combine reference and query samples.
#'
#' @param ref Reference samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param qry Query samples stored in an object of
#'   \code{\link[minfi]{RGChannelSet-class}}.
#' @param batch factor or vector indicating batches.
#' @param batch2 optional factor or vector indicating a second series of
#'   batches.
#' @return An object of \code{\link[minfi]{RGChannelSet-class}} combining ref
#'   and qry. A column \code{ref_query} is added to indicate if a sample is
#'   query or reference.
#' @details The function first adds a column \code{ref_query} is added to
#'   indicate if a sample is query or reference. Then, it distinguish NAs in the
#'   columns batch and batch2 of query and reference samples. If it is a query
#'   sample, assign it to "qNA"; if a reference sample, assign it to "rNA". At
#'   last, only the same columns of phenotype information are kept, and then
#'   calls \code{\link[minfi]{combineArrays}} to combine query and reference.
#' @noRd
.combine <- function(ref, qry, batch, batch2) {
  qry_df <- minfi::pData(qry)
  ref_df <- minfi::pData(ref)
  # Add ref_query column
  qry_df$ref_query <- "query"
  ref_df$ref_query <- "ref"
  common_colnames <- colnames(qry_df)[colnames(qry_df) %in% colnames(ref_df)]
  # Distinguish NAs from query and reference if batch and batch2 are set.
  if (!is.null(batch)) {
    if (!batch %in% common_colnames)
      stop("batch is not present in both query and reference.")
    qry_batch <- qry_df[, batch]
    ref_batch <- ref_df[, batch]
    qry_df[, batch] <- ifelse(is.na(qry_batch), "qNA", paste0("q", qry_batch))
    ref_df[, batch] <- ifelse(is.na(ref_batch), "rNA", paste0("q", ref_batch))
  }
  if (!is.null(batch2)) {
    if (!batch2 %in% common_colnames)
      stop("batch2 is not present in both query and reference.")
    qry_batch2 <- qry_df[, batch2]
    ref_batch2 <- ref_df[, batch2]
    qry_df[, batch2] <- ifelse(is.na(qry_batch2), "qNA", paste0("q", qry_batch2))
    ref_df[, batch2] <- ifelse(is.na(ref_batch2), "rNA", paste0("q", ref_batch2))
  }
  # Phenotype data should have the same columns.
  common_colnames <- colnames(qry_df)[colnames(qry_df) %in% colnames(ref_df)]
  minfi::pData(qry) <- qry_df[, common_colnames]
  minfi::pData(ref) <- ref_df[, common_colnames]
  # Combine
  minfi::combineArrays(qry, ref)
}
