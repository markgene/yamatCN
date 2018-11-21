# Report pipeline methods.

#' Report an object of pipeline classes.
#'
#' @param x An object.
#' @param outdir A character scalar of output directory.
#' @param overwrite A logical scalar. Default to FALSE.
#' @param verbose A logical scalar. Default to TRUE.
#' @param ... Any other arguments.
#' @return A character vector of file names.
#' @export
setGeneric(
  name = "report_pipe",
  def = function(x,
                 outdir,
                 filename,
                 overwrite = FALSE,
                 verbose = TRUE,
                 ...) {
    standardGeneric("report_pipe")
  }
)


#' @describeIn report_pipe x is \code{CwobPipe}.
setMethod(
  f = "report_pipe",
  signature = c(x = "CwobPipe"),
  definition = function(x,
                        outdir,
                        overwrite = FALSE,
                        verbose = TRUE) {

  }
)


#' To \code{minfi} sample name.
#'
#' The function translates sample IDs of \code{DNAcopy} to those of the input
#' \code{minfi} object. If a sample ID begins with a digital number, the
#' \code{DNAcopy} package will add "X" to the beginning.
#'
#' @param x An object of \code{\link{CwobPipe}} class.
#' @param sample_id A character scalar of sample ID used in \code{DNAcopy} result.
#' @return A character scalar of sample name used in \code{minfi} object.
#' @noRd
.to_minfi_sample_name <- function(x, sample_id) {
  sample_ids <- minfi::sampleNames(x@dat$minfi)
  if (!sample_id %in% sample_ids) {
    x_sample_ids <- paste0("X", sample_ids)
    if (sample_id %in% x_sample_ids) {
      sample_id <- sample_ids[sample_id == x_sample_ids]
    } else {
      stop("Fail to find sample ID in DNAcopy result in minfi object.")
    }
  }
  sample_id
}


#' To \code{DNAcopy} sample IDs.
#'
#' The function translates sample IDs used in \code{minfi} object to those
#' in \code{DNAcopy} result. If a sample ID begins with a digital number, the
#' \code{DNAcopy} package will add "X" to the beginning.
#'
#' @param minfi_obj A minfi object such as an object of
#'     \code{\link[minfi]{GenomicRatioSet-class}} or
#'     \code{\link[minfi]{GenomicMethylSet-class}}.
#' @return A character vector of sample IDs used in \code{DNAcopy}.
#' @noRd
.to_DNAcopy_sample_ids <- function(minfi_obj) {
  sample_ids <- minfi::sampleNames(minfi_obj)
  ifelse(stringr::str_detect(sample_ids, "^[0-9]"),
         paste0("X", sample_ids),
         sample_ids)
}
