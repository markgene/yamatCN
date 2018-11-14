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


#' Output Conumee analysis result as table of one sample
#'
#' @param cnv_anl A \code{CNV.analysis} object (single sample)
#' @param outdir A character scalar of output directory
#' @param what A character scalar. This should be (an unambiguous
#' abbreviation of) one of 'probes', 'bins', 'detail' or 'segments'.
#' Defaults to 'segments'.
#' @param outfile A character scalar of output file name.
#' @param prefix A character scalar of output file prefix. Default is
#' sample name. The output file name is the prefix + "-" + \code{what}
#' argument + ".csv". For example, if sample name is "MA001" and
#' \code{what} is "probes", the output file name will be "MA001-probes.csv".
#' @param overwrite A logical scalar. Default to FALSE.
#' @return A character scalar of output file path.
#' @noRd
.write_conumee_result <-
  function(cnv_anl,
           outdir,
           what = "segments",
           prefix = NULL,
           outfile = NULL,
           overwrite = FALSE) {
    if (is.null(outfile)) {
    if (is.null(prefix)) {
      outfile <- paste0(names(cnv_anl), "-", what, ".csv")
    } else {
      outfile <- paste0(prefix, what, ".csv")
    }
  }
  out_fp <- file.path(outdir, outfile)
  if (file.exists(out_fp) & !overwrite) {
    message("Output file ", out_fp,  " has existed. Skip writing.")
  } else {
    df <- conumee::CNV.write(cnv_anl, what = what)
    write.csv(df, file = out_fp, row.names = FALSE, quote = FALSE)
  }
  out_fp
}


