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
#'   \item \code{batch}: a character scalar of batch.
#'   \item \code{batch2}: a character scalar of batch2.
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
#' }
#' @slot meta A list contain meta information.
#' @export
MethylCNVPipe <- setClass(
  "MethylCNVPipe",
  slots = c(dat = "environment", meta = "list"),
  contains = "Pipe"
)
