# Save pipeline methods.

#' Save an object of pipeline classes.
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
  name = "save_pipe",
  def = function(x,
                 outdir,
                 filename,
                 overwrite = FALSE,
                 verbose = TRUE,
                 ...) {
    standardGeneric("save_pipe")
  }
)


#' @describeIn save_pipe x is \code{ConumeePipe}.
setMethod(
  f = "save_pipe",
  signature = c(x = "ConumeePipe"),
  definition = function(x,
                        outdir,
                        filename = "conumee.Rda",
                        overwrite = FALSE,
                        verbose = TRUE) {
    .save_pipe(x,
               outdir,
               filename = filename,
               overwrite = overwrite,
               verbose = verbose)
  }
)


#' @describeIn save_pipe x is \code{MethylCNVPipe}.
setMethod(
  f = "save_pipe",
  signature = c(x = "MethylCNVPipe"),
  definition = function(x,
                        outdir,
                        filename = "methylcnv.Rda",
                        overwrite = FALSE,
                        verbose = TRUE) {
    .save_pipe(x,
               outdir,
               filename = filename,
               overwrite = overwrite,
               verbose = verbose)
  }
)


.save_pipe <- function(x, outdir, filename, overwrite = FALSE, verbose = TRUE) {
  if (missing(x))
    stop("Require argument x")
  if (missing(outdir))
    stop("Require argument outdir")
  if (missing(filename))
    stop("Require argument filename")
  if (!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE)
  filename <- file.path(outdir, filename)
  if (overwrite | !file.exists(filename)) {
    save(x, file = filename)
  } else {
    if (verbose)
      message("File exists. Skip.")
  }
}
