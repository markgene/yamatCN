# Write segments methods.

#' Write segments by sample.
#'
#' @param x An object.
#' @param outdir A character scalar of output directory.
#' @param filename A character scalar of output file name.
#' @param overwrite A logical scalar. Default to FALSE.
#' @param verbose A logical scalar. Default to TRUE.
#' @param ... Any other arguments.
#' @return A character vector of file names.
#' @export
setGeneric(
  name = "write_segments",
  def = function(x,
                 outdir,
                 filename,
                 overwrite = FALSE,
                 verbose = TRUE,
                 ...) {
    standardGeneric("write_segments")
  }
)


#' @describeIn write_segments x is \code{ConumeePipe}.
setMethod(
  f = "write_segments",
  signature = c(x = "ConumeePipe"),
  definition = function(x,
                        outdir,
                        filename = "conumee-segments.tab",
                        overwrite = FALSE,
                        verbose = TRUE) {

  }
)



