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

