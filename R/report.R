# Report pipeline methods.

#' Report an object of pipeline classes.
#'
#' @param x An object.
#' @param outdir A character scalar of output directory.
#' @param genome_plot_file A character scalar of the genome plot file which
#'   shows all chromosomes.
#' @param genome_plot_width An integer scalar. Default to 9.
#' @param genome_plot_height An integer scalar. Default to 15.
#' @param chr_per_row An integer scalar of the number of chromosomes per row
#'   in genome plot. Default to 4.
#' @param chr_plot_width An integer scalar of the width of individual chromosome
#'   plot. Default to 9.
#' @param chr_plot_height An integer scalar of the height of individual
#'   chromosome plot. Default to 6.
#' @param size An integer scalar of the threshold to define CNVs. Default to
#'   \code{5e6} (5 Mb).
#' @param cn_boundary An numeric vector of length two which defines the
#'   boundaries to define CNVs. If a segment has a segment mean lower than
#'   the first element or higher than the second element, it is a CNV. It is
#'   the absolute value, so 2 means no copy number change. Default to
#'   \code{c(1.8, 2.2)}.
#' @param segment_file A character scalar of segments file which is a table of
#'   segments. Default to "segments.tab".
#' @param overwrite A logical scalar. Default to FALSE.
#' @param verbose A logical scalar. Default to TRUE.
#' @param ... Any other arguments.
#' @return TBD.
#' @export
setGeneric(
  name = "report_pipe",
  def = function(x,
                 outdir,
                 genome_plot_file = "genome-plot.png",
                 genome_plot_width = 9,
                 genome_plot_height = 15,
                 chr_per_row = 4,
                 chr_plot_width = 9,
                 chr_plot_height = 6,
                 size = 5e6,
                 cn_boundary = c(1.8, 2.2),
                 segment_file = "segments.tab",
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
                        genome_plot_file = "genome-plot-cwob.png",
                        genome_plot_width = 9,
                        genome_plot_height = 15,
                        chr_per_row = 4,
                        chr_plot_width = 9,
                        chr_plot_height = 6,
                        size = 5e6,
                        cn_boundary = c(1.8, 2.2),
                        segment_file = "segments.tab",
                        overwrite = FALSE,
                        verbose = TRUE) {
    print(genome_plot_file)
    qry_idx <- .query_indices(x@dat$minfi, label = "query")
    sample_dirs <- yamat::init_report(x@dat$minfi[, qry_idx], outdir)
    sample_ids_dnacopy <- .to_DNAcopy_sample_ids(x@dat$minfi[, qry_idx])
    # Gender
    gender <- .gender(x@dat$minfi[, qry_idx])
    # Plot genome
    if (verbose) {
      message("Plotting genomes...")
    }
    genome_plot_files <- sapply(
      seq(length(sample_dirs)),
      function(i) {
        if (verbose) {
          message("Plotting genome for sample ", names(sample_dirs[i]))
        }
        plot_file <- file.path(sample_dirs[i], genome_plot_file)
        if (overwrite | !file.exists(plot_file)) {
          p <- .plot_genome_single_sample(
            x = x,
            sample_id = sample_ids_dnacopy[i],
            gender = gender[i],
            chr_per_row = chr_per_row
          )
          ggplot2::ggsave(
            filename = plot_file,
            plot = p,
            height = genome_plot_height,
            width = genome_plot_width
          )
        }
        plot_file
      }
    )
    # Plot chromosomes contains CNVs.
    if (verbose) {
      message("Plotting chromosomes...")
    }
    chr_plot_files <- lapply(
      seq(length(sample_dirs)),
      function(i) {
        if (verbose) {
          message("Plotting chromosomes for sample ", names(sample_dirs[i]))
        }
        plot_lst <- .plot_chromosomes_single_sample(
          x = x,
          sample_id = sample_ids_dnacopy[i],
          gender = gender[i],
          size = size,
          cn_boundary = cn_boundary
        )
        if (length(plot_lst)) {
          sapply(
            seq(length(plot_lst)),
            function(k) {
              plot_file <- paste0(names(plot_lst)[k], "-cwob", ".png")
              plot_file <- file.path(sample_dirs[i], plot_file)
              if (overwrite | !file.exists(plot_file)) {
                if (verbose)
                  message("Saving to ", plot_file)
                ggplot2::ggsave(
                  filename = plot_file,
                  plot = plot_lst[[k]],
                  height = chr_plot_height,
                  width = chr_plot_width
                )
              }
              plot_file
            }
          )
        }
      }
    )
    # Segment table
    if (verbose) {
      message("Writing segment table...")
    }
    segment_files <- lapply(
      seq(length(sample_dirs)),
      function(i) {
        if (verbose) {
          message("Writing segment table for sample ", names(sample_dirs[i]))
        }
        z <- x@dat$segments %>%
          dplyr::filter(ID == sample_ids_dnacopy[i]) %>%
          dplyr::select(-ID)
        if ("lrr_shift" %in% names(x@dat)) {
          z <- dplyr::mutate(z, CN = 2 ** (seg.mean + 1 - x@dat$lrr_shift))
        } else {
          z <- dplyr::mutate(z, CN = 2 ** (seg.mean + 1))
        }
        if (gender[i] == "F") {
          z <- dplyr::filter(z, chromosome != "chrY")
        }
        outfile <- file.path(sample_dirs[i], segment_file)
        write.table(x = z, file = outfile, quote = FALSE, row.names = FALSE)
      }
    )
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


#' Get gender.
#'
#' @param minfi_obj A minfi object such as an object of
#'     \code{\link[minfi]{GenomicRatioSet-class}} or
#'     \code{\link[minfi]{GenomicMethylSet-class}}.
#' @return A character vector of gender.
#' @noRd
.gender <- function(minfi_obj) {
  pheno_df <- as.data.frame(minfi::pData(minfi_obj))
  colnames(pheno_df) <- toupper(colnames(pheno_df))
  if ("GENDER" %in% colnames(pheno_df)) {
    gender <- pheno_df$GENDER
  } else if ("SEX" %in% colnames(pheno_df)) {
    gender <- pheno_df$SEX
  } else {
    minfi::getSex(minfi_obj, cutoff = -2) %>%
      dplyr::select(predictedSex) %>%
      unlist() -> gender
  }
  if (!all(gender %in% c("M", "F"))) {
    stop("Require to encode gender in single character M for male or F for female.")
  }
  gender
}
