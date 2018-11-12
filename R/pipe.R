# Misc function of copy-number pipelines.

#' Check arguments of pipelines.
#'
#' @param report_dir A character scalar of reporting.
#' @return TRUE if pass all checks (invisible).
#' @noRd
.check_args_pipe <-
  function(report_dir) {
    if (missing(report_dir))
      stop("report_dir is required!")
    invisible(TRUE)
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
