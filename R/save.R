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
    sample_dirs <- yamat::init_report(x@dat$minfi_obj, outdir)
    qry_names <-
      .query_sample_names(x@dat$minfi_obj, label = "query")
    sapply(seq(length(qry_names)),
           function(i) {
             cnv_anl_obj <- x@dat$conumee_results[[qry_names[i]]]
             minfi_obj <- x@dat$minfi_obj[, i]
             filename <- file.path(sample_dirs[i], filename)
             if (overwrite | !file.exists(filename)) {
               save(cnv_anl_obj, minfi_obj, file = filename)
             } else {
               if (verbose)
                 message("File exists. Skip.")
             }
             filename
           })
  }
)
