#


#' Construct copy-number single sample plot
#'
#' Given a data frame construct a plot to display raw copy number calls for a
#' single sample. I borrow the function from \code{\link[GenVisR]{cnView}} and
#' fix the bugs.
#'
#' @name cnView
#' @param x Object of class data frame with rows representing copy number calls
#' from a single sample. The data frame must contain columns with the following
#' names "chromosome", "coordinate", "cn", and optionally "p_value"
#' (see details).
#' @param y Object of class data frame with rows representing cytogenetic bands
#' for a chromosome. The data frame must contain columns with the following
#' names "chrom", "chromStart", "chromEnd", "name", "gieStain" for plotting the
#' ideogram (optional: see details).
#' @param z Object of class data frame with row representing copy number segment
#' calls. The data frame must contain columns with the following names
#' "chromosome", "start", "end", "segmean" (optional: see details)
#' @param genome Character string specifying a valid UCSC genome (see details).
#' @param chr Character string specifying which chromosome to plot one of
#' "chr..." or "all"
#' @param CNscale Character string specifying if copy number calls supplied are
#' relative (i.e.copy neutral == 0) or absolute (i.e. copy neutral ==2). One of
#' "relative" or "absolute"
#' @param ideogram_txtAngle Integer specifying the angle of cytogenetic labels
#' on the ideogram subplot.
#' @param ideogram_txtSize Integer specifying the size of cytogenetic labels on
#' the ideogram subplot.
#' @param plotLayer Valid ggplot2 layer to be added to the copy number plot.
#' @param ideogramLayer Valid ggplot2 layer to be added to the ideogram
#' sub-plot.
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @param cn_boundary An numeric vector of length two which defines the
#'   boundaries to define CNVs. If a segment has a segment mean lower than
#'   the first element or higher than the second element, it is a CNV. It is
#'   the absolute value, so 2 means no copy number change. Default to
#'   \code{c(1.8, 2.2)}.
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
cnView <- function(x, y=NULL, z=NULL, genome='hg19', chr='chr1',
                   CNscale="absolute", ideogram_txtAngle=45,
                   ideogram_txtSize=5, plotLayer=NULL, ideogramLayer=NULL,
                   out="plot", cn_boundary = c(1.8, 2.2))
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
  if(chr == 'all')
  {
    # plot the graphic
    p1 <- cnView_buildMain(x, z=z, dummyData, chr=chr, CNscale = CNscale, cn_boundary = cn_boundary)
  } else {
    # plot chromosome
    chromosome_plot <- GenVisR:::ideoView(cytobands, chromosome=chr,
                                txtAngle=ideogram_txtAngle,
                                txtSize=ideogram_txtSize,
                                plotLayer=ideogramLayer)

    # if requested plot only selected chromosome
    x <- GenVisR:::multi_subsetChr(x, chr)
    if(!is.null(z))
    {
      z <- GenVisR:::multi_subsetChr(z, chr)
    }

    # build the plot
    CN_plot <- cnView_buildMain(x, dummyData, z=z, chr=chr, CNscale=CNscale,
                                layers=plotLayer, cn_boundary = cn_boundary)
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
#'   \code{c(1.8, 2.2)}.
#' @return ggplot2 object
#' @noRd
cnView_buildMain <-
  function(x,
           y,
           z = NULL,
           chr,
           CNscale = c("absolute", "relative"),
           layers = NULL,
           cn_boundary = c(1.8, 2.2)) {
  # Define various parameters of the plot
  CNscale <- match.arg(CNscale)
  dummy_data <-
    ggplot2::geom_point(
      data = y,
      mapping = ggplot2::aes_string(x = 'coordinate', y = 2),
      alpha = 0
    )

  theme <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
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
  if (any('p_value' %in% colnames(x)))
  {
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
            alpha = 0.5 + 1 / (1 + exp(cn - 2))
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
            alpha = 0.5 + 1 / (1 + exp(cn))
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
        color = "white",
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
    dummy_data +
    transparency +
    layers

  if (chr == 'all') {
    facet <- ggplot2::facet_wrap( ~ chromosome, scales = 'free')
    p1 <- p1 + facet
  }

  return(p1)
}
