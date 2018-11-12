# An S4 class to represent Conumee pipeline.

#' Conumee CNV pipeline
#'
#' @slot dat An environment which has the following objects. \code{minfi_obj},
#'   an object of \code{\link[minfi]{GenomicRatioSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}}
#' @export
ConumeePipe <- setClass(
  "ConumeePipe",
  slots = c(dat = "environment")
)


setValidity(
  "ConumeePipe",
  method = function(object) {
    if (!is(dat, "environment"))
      paste("dat slot must be an object of environment class.")
  }
)


setMethod(
  f = "initialize",
  signature = "ConumeePipe",
  function(.Object,
           dat = new.env(parent = parent.frame()),
           ...) {
    .Object <- callNextMethod(.Object, ...)
    .Object@dat <- dat
    .Object
  }
)


setMethod(
  f = "show",
  signature = "ConumeePipe",
  function(object) {
    cat("ConumeePipe object\n")
  }
)
