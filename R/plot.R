# Plot CNVs.

.plot_single_sample_prep <- function(x, sample_id, gender = c("M", "F")) {
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
