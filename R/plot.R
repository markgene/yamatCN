# Plot CNVs.

#' Genome plot for single sample.
#' @noRd
.plot_genome_single_sample <- function(x, sample_id, gender = c("M", "F"), chr_per_row = 4) {
  dat <- .plot_single_sample_prep(x, sample_id, gender = gender)
  cnView(dat$lrr, z = dat$z, chr = "all", genome = "hg19", CNscale = "absolute") +
    ggplot2::facet_wrap(~chromosome, scales = "free_x", ncol = chr_per_row)
}


#' Plot chromosomes for single sample.
#' @noRd
.plot_chromosomes_single_sample <-
  function(x,
           sample_id,
           gender = c("M", "F"),
           size = 5e6,
           cn_boundary = c(1.5, 2.5),
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
    # list(plots = plot_lst, chr = sele_chr)
  }


#' Prepare for the data of a single sample.
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
    dplyr::rename(chromosome = chrom, start = loc.start, end = loc.end, segmean = seg.mean) %>%
    #dplyr::mutate(segmean = 2 ** (segmean + 1)) %>%
    dplyr::mutate(segmean = 2 ** (segmean + 1 - x@dat$lrr_shift)) %>%
    dplyr::select(chromosome, start, end, segmean)
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
