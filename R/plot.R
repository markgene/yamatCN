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
#'   \code{c(1.319508, 2.828427)}.
#' @param chr_per_row An integer scalar of chromosome per row in the plot.
#'   Default to 4.
#' @return A \code{\link{ggplot2}{ggplot}} object.
#' @noRd
.plot_genome_single_sample <-
  function(x,
           sample_id,
           gender = c("M", "F"),
           cn_boundary = c(1.319508, 2.828427),
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
#'   \code{c(1.319508, 2.828427)}.
#' @param verbose A logical scalar. Default to TRUE.
#' @return A list of \code{\link[gtable]{gtable}} class.
#' @noRd
.plot_chromosomes_single_sample <-
  function(x,
           sample_id,
           gender = c("M", "F"),
           size = 5e6,
           cn_boundary = c(1.319508, 2.828427),
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
          CNscale = "absolute",
          cn_boundary = cn_boundary
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
    z <- dplyr::mutate(z, segmean = .to_absolute(segmean - x@dat$lrr_shift))
  } else {
    z <- dplyr::mutate(z, segmean = .to_absolute(segmean))
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
  if ("lrr_shift" %in% names(x@dat)) {
    cna_df <- dplyr::mutate(cna_df, cn = .to_absolute(cn - x@dat$lrr_shift))
  } else {
    cna_df <- dplyr::mutate(cna_df, cn = .to_absolute(cn))
  }
  cna_df %>%
    dplyr::mutate(chromosome = forcats::fct_relevel(chromosome, paste0("chr", c(1:22, "X", "Y")))) %>%
    dplyr::select(chromosome, coordinate, cn) -> cna_df
  if (gender == "F") {
    cna_df <- dplyr::filter(cna_df, chromosome != "chrY")
  }
  list(lrr = cna_df, z = z)
}


#' Prepare for the data of a single sample for \code{\link[conumee]{CNV.analysis}}
#' object.
#'
#' @param x A \code{\link[conumee]{CNV.analysis}} object.
#' @param gender A character scalar of gender. It is coded as either "M" or "F".
#' @param scaled Character string specifying if copy number calls supplied are
#' relative (i.e.copy neutral == 0) or absolute (i.e. copy neutral ==2). One of
#' "relative" or "absolute"
#' @return A list of two elements:
#'   \itemize{
#'     \item \code{lrr}: A \code{data.frame} of LRRs which has three columns,
#'       \code{chromosome}, \code{coordinate}, and \code{cn}. \code{cn} is
#'       the absolute value of copy number, i.e. 2 means two copies without any
#'       alternation.
#'     \item \code{loci}: A \code{data.frame} of LRRs of loci which has three
#'       columns, \code{chromosome}, \code{coordinate}, and \code{cn}. \code{cn}
#'       is the absolute value of copy number, i.e. 2 means two copies without
#'       any alternation.
#'     \item \code{z}: A\code{data.fram} of segment which has four columns,
#'       \code{chromosome}, \code{start}, \code{end} and \code{segmean}.
#'   }
#' @noRd
.plot_prep <- function(x, gender = c("M", "F"), scaled = c("relative", "absolute")) {
  if (!inherits(x, "CNV.analysis")) {
    stop("Argument x should be an object of CNV.analysis class.")
  }
  if (!gender %in% c("M", "F")) {
    stop("Argument gender should be either M or F.")
  }
  scaled <- match.arg(scaled)
  z <- x@seg$summary %>%
    dplyr::rename(chromosome = chrom, start = loc.start, end = loc.end, segmean = seg.mean)
  if (scaled == "absolute") {
    z <- dplyr::mutate(z, segmean = .to_absolute(segmean - x@bin$shift))
  } else if (scaled == "relative") {
    z <- dplyr::mutate(z, segmean = segmean - x@bin$shift)
  }
  if (gender == "F") {
    z <- dplyr::filter(z, chromosome != "chrY")
  }
  # LRR of bins
  if (scaled == "absolute") {
    cn <- .to_absolute(x@bin$ratio - x@bin$shift)
  } else if (scaled == "relative") {
    cn <- x@bin$ratio - x@bin$shift
  }
  cna_df <-
    data.frame(
      chromosome = GenomicRanges::seqnames(x@anno@bins),
      coordinate = GenomicRanges::values(x@anno@bins)$midpoint,
      cn = cn
    )
  if (gender == "F") {
    cna_df <- dplyr::filter(cna_df, chromosome != "chrY")
  }
  # LRR of loci
  loci_names <- names(x@anno@probes)
  if (scaled == "absolute") {
    loci_cn <- .to_absolute(x@fit$ratio[loci_names] - x@bin$shift)
  } else if (scaled == "relative") {
    loci_cn <- x@fit$ratio[loci_names] - x@bin$shift
  }
  loci_df <-
    data.frame(
      chromosome = GenomicRanges::seqnames(x@anno@probes),
      coordinate = GenomicRanges::start(x@anno@probes),
      cn = loci_cn
    )
  if (gender == "F") {
    loci_df <- dplyr::filter(loci_df, chromosome != "chrY")
  }
  list(lrr = cna_df, z = z, loci = loci_df)
}


#' Construct copy-number single sample plot
#'
#' Given a data frame construct a plot to display raw copy number calls for a
#' single sample. I borrow the function from \code{\link[GenVisR]{cnView}} and
#' fix the bugs.
#'
#' @name cnView
#' @param x Object of class data frame with rows representing copy number calls
#'   from a single sample. The data frame must contain columns with the following
#'   names "chromosome", "coordinate", "cn", and optionally "p_value"
#'   (see details).
#' @param y Object of class data frame with rows representing cytogenetic bands
#'   for a chromosome. The data frame must contain columns with the following
#'   names "chrom", "chromStart", "chromEnd", "name", "gieStain" for plotting the
#'   ideogram (optional: see details).
#' @param z Object of class data frame with row representing copy number segment
#'   calls. The data frame must contain columns with the following names
#'   "chromosome", "start", "end", "segmean" (optional: see details)
#' @param genome Character string specifying a valid UCSC genome (see details).
#' @param chr Character string specifying which chromosome to plot one of
#'   "chr..." or "all".
#' @param CNscale Character string specifying if copy number calls supplied are
#'   relative (i.e.copy neutral == 0) or absolute (i.e. copy neutral ==2). One
#'   of "relative" or "absolute".
#' @param ideogram_txtAngle Integer specifying the angle of cytogenetic labels
#'   on the ideogram subplot.
#' @param ideogram_txtSize Integer specifying the size of cytogenetic labels on
#'   the ideogram subplot.
#' @param plotLayer Valid ggplot2 layer to be added to the copy number plot.
#' @param ideogramLayer Valid ggplot2 layer to be added to the ideogram
#'   sub-plot.
#' @param out Character vector specifying the the object to output, one of
#'   "data", "grob", or "plot", defaults to "plot" (see returns).
#' @param cn_boundary An numeric vector of length two which defines the
#'   boundaries to define CNVs. If a segment has a segment mean lower than
#'   the first element or higher than the second element, it is a CNV. It is
#'   the absolute value, so 2 means no copy number change. Default to
#'   \code{c(1.319508, 2.828427)}.
#' @param gr A genomic range of \code{\link[GenomicRanges]{GRanges}} class.
#'   If it is provided, neglect \code{chr} and plot the region defined. Default
#'   to \code{NULL}.
#' @details cnView is able to plot in two modes specified via the `chr`
#' parameter, these modes are single chromosome view in which an ideogram is
#' displayed and genome view where chromosomes are faceted. For the single
#' chromosome view cytogenetic band information is required giving the
#' coordinate, stain, and name of each band. As a convenience cnView stores this
#' information for the following genomes "hg19", "hg38", "mm9", "mm10", and
#' "rn5". If the genome assembly supplied to the `genome` parameter is not one
#' of the 5 afore mentioned genome assemblies cnView will attempt to query the
#' UCSC MySQL database to retrieve this information. Alternatively the user can
#' manually supply this information as a data frame to the `y` parameter, input
#' to the `y` parameter take precedence of input to `genome`.
#'
#' cnView is also able to represent p-values for copy-number calls if they are
#' supplied via the "p_value" column in the argument supplied to x. The presence
#' of this column in x will set a transparency value to copy-number calls with
#' calls of less significance becoming more transparent.
#'
#' If it is available cnView can plot copy-number segment calls on top of raw
#' calls supplied to parameter `x` via the parameter `z`.
#' @noRd
cnView <- function(x,
                   y = NULL,
                   z = NULL,
                   genome = 'hg19',
                   chr = 'chr1',
                   CNscale = "absolute",
                   ideogram_txtAngle = 45,
                   ideogram_txtSize = 5,
                   plotLayer = NULL,
                   ideogramLayer = NULL,
                   out = "plot",
                   cn_boundary = c(1.319508, 2.828427),
                   gr = NULL
)
{
  # Perform a basic quality check
  input <- GenVisR:::cnView_qual(x, y, z, genome, CNscale=CNscale)
  x <- input[[1]]
  y <- input[[2]]
  z <- input[[3]]

  # Obtain Cytogenetic Band information
  # use y input or query UCSC for the data if it's not preloaded
  preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
  if(is.null(y) && any(genome == preloaded))
  {
    message("genome specified is preloaded, retrieving data...")
    cytobands <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == genome,]
    cytobands <- cytobands[,-which(colnames(cytobands) == "genome")]
  } else if(is.null(y)) {
    # Obtain data for UCSC genome and extract relevant columns
    memo <- paste0("attempting to query UCSC mySQL database for chromosome",
                   " positions and cytogenetic information")
    message(memo)
    cytobands <- suppressWarnings(GenVisR:::multi_cytobandRet(genome=genome))
  } else {
    memo <- paste0("Detected argument supplied to y.. using y for",
                   "position and cytogenetic information")
    message(memo)
    cytobands <- y
  }

  # Create Dummy data and add to x for proper plot dimensions
  fakeStart <- stats::aggregate(data=cytobands, FUN=min, chromStart~chrom)
  colnames(fakeStart) <- c("chromosome", "coordinate")
  fakeEnd <- stats::aggregate(data=cytobands, FUN=max, chromEnd~chrom)
  colnames(fakeEnd) <- c("chromosome", "coordinate")
  dummyData <- rbind(fakeStart, fakeEnd)
  dummyData$chromosome <- as.factor(dummyData$chromosome)
  dummyData <- GenVisR:::multi_subsetChr(dummyData, chr)

  # Plot all chromosomes at once if specified
  if (is.null(gr)) {
    if(chr == 'all') {
      # plot the graphic
      p1 <-
        cnView_buildMain(
          x,
          z = z,
          dummyData,
          chr = chr,
          CNscale = CNscale,
          cn_boundary = cn_boundary
        )
    } else {
      # plot chromosome
      chromosome_plot <-
        GenVisR:::ideoView(
          cytobands,
          chromosome = chr,
          txtAngle = ideogram_txtAngle,
          txtSize = ideogram_txtSize,
          plotLayer = ideogramLayer
        )
      # if requested plot only selected chromosome
      x <- GenVisR:::multi_subsetChr(x, chr)
      if(!is.null(z)) {
        z <- GenVisR:::multi_subsetChr(z, chr)
      }
      # build the plot
      CN_plot <-
        cnView_buildMain(
          x,
          dummyData,
          z = z,
          chr = chr,
          CNscale = CNscale,
          layers = plotLayer,
          cn_boundary = cn_boundary
        )
    }
  } else {
    # Use the first one if there are many GRanges.
    gr <- gr[1]
    gr_df <- GenomicRanges::mcols(gr)
    if (!"name" %in% colnames(gr_df)) {
      stop("Require name column in meta data frame of gr.")
    }
    if (!"thick" %in% colnames(gr_df)) {
      stop("Require thick column in meta data frame of gr.")
    } else {
      if (!inherits(gr_df$thick, "IRanges")) {
        stop("Thick column in meta data frame of gr should be IRanges class.")
      }
    }
    p1 <-
      cnView_buildMain(
        x,
        z = z,
        dummyData,
        chr = as.character(GenomicRanges::seqnames(gr)),
        CNscale = CNscale,
        cn_boundary = cn_boundary
      ) +
      ggplot2::lims(x = c(
        GenomicRanges::start(gr_df$thick),
        GenomicRanges::end(gr_df$thick)
      )) +
      ggplot2::geom_segment(
        x = GenomicRanges::start(gr),
        y = -0.05,
        xend = GenomicRanges::end(gr),
        yend = -0.05,
        size = 3,
        colour = "red",
        alpha = 0.9
      ) +
      ggplot2::xlab(as.character(gr_df$name)) +
      ggplot2::theme(legend.position = "top")
  }

  # Decide what to output
  dataOut <- list(main=x, dummyData=dummyData, segments=z, cytobands=cytobands)
  if(!exists("p1", inherits=FALSE))
  {
    p1 <- GenVisR:::multi_align(chromosome_plot, CN_plot)
    output <- GenVisR:::multi_selectOut(data=dataOut, plot=p1, draw=FALSE, out=out)
  } else {
    output <- GenVisR:::multi_selectOut(data=dataOut, plot=p1, draw=FALSE, out=out)
  }

  return(output)
}


#' construct CN plot
#'
#' given a CN data frame plot points in ggplot. I borrow the function from
#' \code{\link[GenVisR]{cnView_buildMain}} and color segment by segment mean.
#'
#' @name cnView_buildMain
#' @param x a data frame with columns chromosome, coordinate, cn, p_value
#' @param y a data frame with columns chromosome, coordinate for plotting
#' boundaries
#' @param z a data frame with columns chromsome, start, end, segmean specifying
#' segments called from copy number (optional)
#' @param chr a character string specifying chromosome
#' @param CNscale Character string specifying if copy number calls supplied are
#' relative (i.e.copy neutral == 0) or absolute (i.e. copy neutral ==2). One of
#' "relative" or "absolute"
#' @param layers additional ggplot2 layers to add.
#' @param cn_boundary An numeric vector of length two which defines the
#'   boundaries to define CNVs. If a segment has a segment mean lower than
#'   the first element or higher than the second element, it is a CNV. It is
#'   the absolute value, so 2 means no copy number change. Default to
#'   \code{c(1.319508, 2.828427)}.
#' @param max_cn A numeric scalar of the maximum of CN. The CNs are sometime
#'   extremely large. The will remove the data points whose CN are greater than
#'   \code{max_cn}. Default to 6.
#' @return ggplot2 object
#' @noRd
cnView_buildMain <-
  function(x,
           y,
           z = NULL,
           chr,
           CNscale = c("absolute", "relative"),
           layers = NULL,
           cn_boundary = c(1.319508, 2.828427),
           max_cn = 6) {
    # Define various parameters of the plot
    CNscale <- match.arg(CNscale)
    dummy_data <-
      ggplot2::geom_point(
        data = y,
        mapping = ggplot2::aes_string(x = 'coordinate', y = 2),
        alpha = 0
      )
    # Replace extremely large CN with max_cn.
    x <- x[x$cn < max_cn, ]

    theme <- ggplot2::theme(
      legend.position = "top",
      # panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.text.x = ggplot2::element_text(angle=30, hjust=1))
    if (CNscale == "relative") {
      # cn fill colors
      shade_cn <- ggplot2::scale_color_gradient2(
        "Copy Number",
        midpoint = 0,
        low = '#009ACD',
        mid = 'gray65',
        high = '#C82536',
        space = 'Lab'
      )
      # y label
      ylabel <- ggplot2::ylab('Copy Number Difference')
    } else if (CNscale == "absolute") {
      # cn fill colors
      shade_cn <- ggplot2::scale_color_gradient2(
        "Copy Number",
        midpoint = 2,
        low = '#009ACD',
        mid = 'gray65',
        high = '#C82536',
        space = 'Lab'
      )
      # y label
      ylabel <- ggplot2::ylab('Absolute Copy Number')
    } else {
      memo <-
        paste0(
          "Did not recognize input to CNscale... defaulting to",
          "absolute scale, please specify \"relative\"",
          "if copy neutral calls == 0"
        )
      warning(memo)
    }

    # provide x label
    xlabel <- ggplot2::xlab('Coordinate')

    # allow an extra layer in the plot
    if(!is.null(layers)) {
      layers <- layers
    } else {
      layers <- ggplot2::geom_blank()
    }

    # if x contains a p_value column set an alpha for it and plot points
    if (any('p_value' %in% colnames(x))) {
      x$transparency <- 1 - x$p_value
      cnpoints <-
        ggplot2::geom_point(
          data = x,
          mapping = ggplot2::aes_string(
            x = 'coordinate',
            y = 'cn',
            colour = 'cn',
            alpha = 'transparency'
          )
        )
      transparency <- ggplot2::scale_alpha(guide = 'none')
    } else {
      if (CNscale == "absolute") {
        cnpoints <-
          ggplot2::geom_point(
            data = x,
            mapping = ggplot2::aes(
              x = coordinate,
              y = cn,
              colour = cn,
              alpha = 1 / (1 + exp(-abs(cn - 2)))
            )
          )
      } else if (CNscale == "relative") {
        cnpoints <-
          ggplot2::geom_point(
            data = x,
            mapping = ggplot2::aes(
              x = coordinate,
              y = cn,
              colour = cn,
              alpha = 1 / (1 + exp(-abs(cn)))
            )
          )
      }
      transparency <- ggplot2::scale_alpha(guide = 'none')
    }

    # Define segments for main plot
    if (!is.null(z)) {
      lower_bound <- cn_boundary[1]
      higher_bound <- cn_boundary[2]
      z %>%
        dplyr::filter(segmean > lower_bound & segmean < higher_bound) %>%
        ggplot2::geom_segment(
          data = .,
          mapping = ggplot2::aes(
            x = start,
            xend = end,
            y = segmean,
            yend = segmean
          ),
          color = "gray80",
          size = 1
        ) -> cnseg
      z %>%
        dplyr::filter(segmean < lower_bound) %>%
        ggplot2::geom_segment(
          data = .,
          mapping = ggplot2::aes(
            x = start,
            xend = end,
            y = segmean,
            yend = segmean
          ),
          color = "#009ACD",
          size = 1
        ) -> cnseg1
      z %>%
        dplyr::filter(segmean > higher_bound) %>%
        ggplot2::geom_segment(
          data = .,
          mapping = ggplot2::aes(
            x = start,
            xend = end,
            y = segmean,
            yend = segmean
          ),
          color = "#C82536",
          size = 1
        ) -> cnseg2
    } else {
      cnseg <- cnseg1 <- cnseg2 <- ggplot2::geom_blank()
    }

    # build the plot
    tmp <- data.frame(x = 0, y = 0)
    p1 <-
      ggplot2::ggplot(data = tmp, ggplot2::aes(x = 0)) +
      cnpoints +
      shade_cn +
      ylabel +
      xlabel +
      ggplot2::theme_bw() +
      theme +
      cnseg +
      cnseg1 +
      cnseg2 +
      ggplot2::lims(y = c(0, max_cn)) +
      ggplot2::scale_y_continuous(breaks = c(-0.5, seq(max_cn))) +
      ggplot2::scale_x_continuous(breaks = NULL) +
      dummy_data +
      transparency +
      layers

    if (chr == 'all') {
      facet <- ggplot2::facet_wrap( ~ chromosome, scales = 'free')
      p1 <- p1 + facet
    }

    return(p1)
  }

