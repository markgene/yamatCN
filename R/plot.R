# Plot CNVs.

#' Genome plot for single sample.
#'
#' @param x An object of \code{\link{CwobPipe}} or \code{\link{MethylCNVPipe}}
#'   class.
#' @param sample_id A character scalar of sample ID.
#' @param gender A character scalar of gender. It is coded as either "M" or "F".
#' @param cn_boundary An numeric vector of length two which defines the
#'   boundaries to define CNVs. If a segment has a segment mean lower than
#'   the first element or higher than the second element, it is a CNV. It is
#'   the absolute value, so 2 means no copy number change. Default to
#'   \code{c(1.8, 2.2)}.
#' @param chr_per_row An integer scalar of chromosome per row in the plot.
#'   Default to 4.
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @noRd
.plot_genome_single_sample <-
  function(x,
           sample_id,
           gender = c("M", "F"),
           cn_boundary = c(1.8, 2.2),
           chr_per_row = 4) {
    dat <- .plot_single_sample_prep(x, sample_id, gender = gender)
    cnView(
      dat$lrr,
      z = dat$z,
      chr = "all",
      genome = "hg19",
      CNscale = "absolute",
      cn_boundary = cn_boundary
    ) +
      ggplot2::facet_wrap( ~ chromosome, scales = "free_x", ncol = chr_per_row)
  }


#' Plot chromosomes for single sample.
#'
#' @param x An object of \code{\link{CwobPipe}} or \code{\link{MethylCNVPipe}}
#'   class.
#' @param sample_id A character scalar of sample ID.
#' @param gender A character scalar of gender. It is coded as either "M" or "F".
#' @param size An integer scalar of the threshold to define CNVs. Default to
#'   \code{5e6} (5 Mb).
#' @param cn_boundary An numeric vector of length two which defines the
#'   boundaries to define CNVs. If a segment has a segment mean lower than
#'   the first element or higher than the second element, it is a CNV. It is
#'   the absolute value, so 2 means no copy number change. Default to
#'   \code{c(1.8, 2.2)}.
#' @param verbose A logical scalar. Default to TRUE.
#' @return A list of \code{\link[gtable]{gtable}} class.
#' @noRd
.plot_chromosomes_single_sample <-
  function(x,
           sample_id,
           gender = c("M", "F"),
           size = 5e6,
           cn_boundary = c(1.8, 2.2),
           verbose = TRUE) {
    dat <- .plot_single_sample_prep(x, sample_id, gender = gender)
    dat$z %>%
      dplyr::filter(end - start + 1 > size) %>%
      dplyr::filter(segmean > cn_boundary[2] | segmean < cn_boundary[1]) %>%
      dplyr::select(chromosome) %>%
      dplyr::distinct() %>%
      unlist() -> sele_chr
    plot_lst <- list()
    if (length(sele_chr)) {
      lapply(sele_chr, function(chr) {
        if (verbose)
          message("Plotting chromosome ", chr)
        cnView(
          dat$lrr,
          z = dat$z,
          chr = chr,
          genome = "hg19",
          CNscale = "absolute"
        )
      }) -> plot_lst
      names(plot_lst) <- sele_chr
    }
    plot_lst
  }


#' Prepare for the data of a single sample.
#'
#' @param x An object of \code{\link{CwobPipe}} or \code{\link{MethylCNVPipe}}
#'   class.
#' @param sample_id A character scalar of sample ID.
#' @param gender A character scalar of gender. It is coded as either "M" or "F".
#' @return A list of two elements:
#'   \itemize{
#'     \item \code{lrr}: A \code{data.frame} of LRRs which has three columns,
#'       \code{chromosome}, \code{coordinate}, and \code{cn}. \code{cn} is
#'       the absolute value of copy number, i.e. 2 means two copies without any
#'       alternation.
#'     \item \code{z}: A\code{data.fram} of segment which has four columns,
#'       \code{chromosome}, \code{start}, \code{end} and \code{segmean}.
#'   }
#' @noRd
.plot_single_sample_prep <- function(x, sample_id, gender = c("M", "F")) {
  if (!inherits(x, "CwobPipe") & !inherits(x, "MethylCNVPipe")) {
    stop("Argument x should be an object of either CwobPipe or MethylCNVPipe class.")
  }
  if (!gender %in% c("M", "F")) {
    stop("Argument gender should be either M or F.")
  }
  if (!sample_id %in% unique(x@dat$segments$ID)) {
    stop("Fail to find sample ID in segment result: ", sample_id)
  }
  z <- x@dat$segments %>%
    dplyr::filter(ID == sample_id) %>%
    dplyr::rename(chromosome = chrom, start = loc.start, end = loc.end, segmean = seg.mean)
  if ("lrr_shift" %in% names(x@dat)) {
    z <- dplyr::mutate(z, segmean = 2 ** (segmean + 1 - x@dat$lrr_shift))
  } else {
    z <- dplyr::mutate(z, segmean = 2 ** (segmean + 1))
  }
  z <- dplyr::select(z, chromosome, start, end, segmean)
  if (gender == "F") {
    z <- dplyr::filter(z, chromosome != "chrY")
  }
  # LRR
  cna_df <- as.data.frame(x@dat$dnacopy$CNA)
  if (!sample_id %in% colnames(cna_df)) {
    stop("Fail to find sample ID in CNA data frame: ", sample_id)
  }
  cna_df <- cna_df[, c("chrom", "maploc", sample_id)]
  colnames(cna_df) <- c("chromosome", "coordinate", "cn")
  cna_df %>%
    dplyr::mutate(cn = 2 ** (cn + 1)) %>%
    dplyr::mutate(chromosome = forcats::fct_relevel(chromosome, paste0("chr", c(1:22, "X", "Y")))) %>%
    dplyr::select(chromosome, coordinate, cn) -> cna_df
  if (gender == "F") {
    cna_df <- dplyr::filter(cna_df, chromosome != "chrY")
  }
  list(lrr = cna_df, z = z)
}
