# S4 classes to represent pipelines.

#' A S4 class represent CNV pipeline.
#'
#' @slot dat An environment which contains data.
#' @slot meta A list contain meta information.
#' @export
Pipe <- setClass(
  "Pipe",
  slots = c(dat = "environment", meta = "list"),
  contains = "VIRTUAL"
)

setValidity(
  "Pipe",
  method = function(object) {
    if (!is(dat, "environment"))
      paste("dat slot must be an object of environment class.")
    if (!is(meta, "list"))
      paste("meta slot must be an object of list class.")
  }
)


setMethod(
  f = "initialize",
  signature = "Pipe",
  function(.Object,
           dat = new.env(parent = parent.frame()),
           meta = list(),
           ...) {
    .Object <- callNextMethod(.Object, ...)
    .Object@dat <- dat
    if (!length(meta)) {
      meta$package_version <- packageVersion("yamatCN")
      meta$created_time <- Sys.time()
    }
    .Object@meta <- meta
    .Object
  }
)


setMethod(
  f = "show",
  signature = "Pipe",
  function(object) {
    cat(paste(class(object), "object\n"))
    dat_names <- paste(names(object@dat), sep = ", ", collapse = ", ")
    cat(paste0("Data:            ", dat_names, "\n"))
    cat(paste0("Package version: ", object@meta$package_version, "\n"))
    cat(paste0("Created time:    ", object@meta$created_time, "\n"))
  }
)


#' Conumee CNV pipeline
#'
#' @slot dat An environment which has the following objects:
#' \itemize{
#'   \item \code{minfi_obj}: an object of
#'     \code{\link[minfi]{GenomicRatioSet-class}} or
#'     \code{\link[minfi]{GenomicMethylSet-class}}.
#'   \item \code{conumee_results}: a list of \code{\link[conumee]{CNV.analysis-class}}.
#' }
#' @slot meta A list contain meta information.
#' @export
ConumeePipe <- setClass(
  "ConumeePipe",
  slots = c(dat = "environment", meta = "list"),
  contains = "Pipe"
)


#' MethylCNV pipeline
#'
#' @slot dat An environment which has the following objects:
#' \itemize{
#'   \item \code{minfi_obj}: an object of
#'     \code{\link[minfi]{GenomicRatioSet-class}} or
#'     \code{\link[minfi]{GenomicMethylSet-class}}.
#'   \item \code{dnacopy_obj} an object of \code{\link[DNAcopy]{DNAcopy}} class.
#'   \item \code{segments} a \code{data.frame} of segments.
#' }
#' @slot meta A list contain meta information.
#' @export
MethylCNVPipe <- setClass(
  "MethylCNVPipe",
  slots = c(dat = "environment", meta = "list"),
  contains = "Pipe"
)


#' CwobPipe pipeline
#'
#' @slot dat An environment which has the following objects:
#' \itemize{
#'   \item \code{minfi_obj}: an object of
#'     \code{\link[minfi]{GenomicRatioSet-class}} or
#'     \code{\link[minfi]{GenomicMethylSet-class}}.
#'   \item \code{dnacopy_obj} an object of \code{\link[DNAcopy]{DNAcopy}} class.
#'   \item \code{segments} a \code{data.frame} of segments.
#'   \item \code{lrr_shift} a numeric scalar. LRR is shifted to minimize the
#'     median absolute deviation from all bins to zero to determine the
#'     copy-number neutral state.
#' }
#' @slot meta A list contain meta information.
#' @export
CwobPipe <- setClass(
  "CwobPipe",
  slots = c(dat = "environment", meta = "list"),
  contains = "Pipe"
)

