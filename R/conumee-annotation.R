#' @import minfi
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import IlluminaHumanMethylation450kmanifest
#' @import IlluminaHumanMethylationEPICanno.ilm10b2.hg19
#' @import IlluminaHumanMethylationEPICanno.ilm10b3.hg19
#' @import IlluminaHumanMethylationEPICanno.ilm10b4.hg19
#' @import IlluminaHumanMethylationEPICmanifest
#' @importFrom rtracklayer import
NULL


#' Create Conumee annotation for CNV analysis.
#'
#' Rewrite \code{\link[conumee]{CNV.create_anno}()} to be compatible with
#' ilm10b3.hg19 and ilm10b4.hg19 annotation of EPIC array. I also change the
#' default array type to EPIC.
#'
#' @details This function collects all annotations required for CNV analysis
#' using Illumina 450k or EPIC/850K arrays. I rewrite the function so that
#' it can generate the annotation for the list of probes specified by the
#' experiments.
#'
#' @param array_type A character scalar. One of "450k", "EPIC", or "overlap".
#'   Defaults to "EPIC".
#' @param array_anno A character scalar. One of "ilm10b2.hg19", "ilm10b3.hg19",
#'   "ilm10b4.hg19". Only valid when \code{array_type} is set to "EPIC". Default
#'   to "ilm10b4.hg19".
#' @param bin_minprobes numeric. Minimum number of probes per bin. Bins
#'   are interatively merged with neighboring bin until minimum number is
#'   reached.
#' @param bin_minsize numeric. Minimum size of a bin.
#' @param bin_maxsize numeric. Maximum size of a bin. Merged bins that
#'   are larger are filtered out.
#' @param chrXY logical. Should chromosome X and Y be included in the analysis?
#' @param exclude_regions GRanges object or path to bed file containing
#'   genomic regions to be excluded.
#' @param detail_regions GRanges object or path to bed file containing
#'   genomic regions to be examined in detail.
#' @return An object of \code{\link[conumee]{CNV.anno-class}}.
#' @note I have tested with EPIC ilm10b4.hg19, but not other arrays and
#'   annotations.
CNV.create_anno2 <-
  function (array_type = "EPIC",
            array_anno = "ilm10b4.hg19",
            bin_minprobes = 15 * 2,
            bin_minsize = 50000 * 2,
            bin_maxsize = 5e+06 * 2,
            chrXY = FALSE,
            exclude_regions = NULL,
            detail_regions = NULL)
  {
    object <- new("CNV.anno")
    object@date <- date()
    a1 <- formals()
    a2 <- as.list(match.call())[-1]
    getElt <- function(an) {
      if (is.element(an, names(a2))) {
        a2[[an]]
      }
      else {
        a1[[an]]
      }
    }
    object@args <- as.list(sapply(unique(names(c(
      a1, a2
    ))), getElt,
    simplify = F))
    if (array_type == "overlap" |
        array_type == "450k" |
        (array_type == "EPIC" & array_anno == "ilm10b2.hg19")) {
      return(
        conumee::CNV.create_anno(
          array_type = array_type,
          bin_minprobes = 15,
          bin_minsize = 50000,
          bin_maxsize = 5e+06,
          chrXY = FALSE,
          exclude_regions = NULL,
          detail_regions = NULL
        )
      )
    } else {
      object@args[["bin_minsize"]] <- bin_minsize
      object@args[["array_type"]] <- array_type
      object@args[["chrXY"]] <- chrXY
    }
    if (chrXY) {
      object@genome <- data.frame(chr = paste("chr", c(1:22,
                                                       "X", "Y"), sep = ""),
                                  stringsAsFactors = FALSE)
    }
    else {
      object@genome <- data.frame(chr = paste("chr", 1:22,
                                              sep = ""),
                                  stringsAsFactors = FALSE)
    }
    rownames(object@genome) <- object@genome$chr
    message("using genome annotations from UCSC")
    tbl.chromInfo <- conumee:::tbl_ucsc$chromInfo[match(object@genome$chr,
                                                        conumee:::tbl_ucsc$chromInfo$chrom), "size"]
    object@genome$size <- tbl.chromInfo
    tbl.gap <-
      conumee:::tbl_ucsc$gap[is.element(conumee:::tbl_ucsc$gap$chrom, object@genome$chr),]
    object@gap <-
      sort(GRanges(
        as.vector(tbl.gap$chrom),
        IRanges(tbl.gap$chromStart +
                  1, tbl.gap$chromEnd),
        seqinfo = Seqinfo(object@genome$chr,
                          object@genome$size)
      ))
    tbl.cytoBand <-
      conumee:::tbl_ucsc$cytoBand[is.element(conumee:::tbl_ucsc$cytoBand$chrom,
                                             object@genome$chr),]
    pq <-
      sapply(split(tbl.cytoBand$chromEnd[grepl("p", tbl.cytoBand$name)],
                   as.vector(tbl.cytoBand$chrom[grepl("p", tbl.cytoBand$name)])),
             max)
    object@genome$pq <- start(resize(subsetByOverlaps(object@gap,
                                                      GRanges(
                                                        names(pq), IRanges(pq, pq)
                                                      )), 1, fix = "center"))
    probes450k <- probesEPIC <- GRanges()
    if (array_anno == "ilm10b4.hg19") {
      message("getting EPIC annotation ilm10b4.hg19")
      probesEPIC <-
        sort(
          minfi::getLocations(
            IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19
          )
        )
    } else if (array_anno == "ilm10b3.hg19") {
      message("getting EPIC annotation ilm10b3.hg19")
      probesEPIC <-
        sort(
          minfi::getLocations(
            IlluminaHumanMethylationEPICanno.ilm10b3.hg19::IlluminaHumanMethylationEPICanno.ilm10b3.hg19
          )
        )
    }
    if (array_type == "overlap") {
      probes <- sort(subsetByOverlaps(probes450k, probesEPIC))
    }
    else {
      probes <- unique(sort(c(probes450k, probesEPIC)))
    }
    probes <-
      probes[substr(names(probes), 1, 2) == "cg" &
               is.element(as.vector(seqnames(probes)),
                          object@genome$chr)]
    object@probes <- sort(GRanges(
      as.vector(seqnames(probes)),
      ranges(probes),
      seqinfo = Seqinfo(object@genome$chr,
                        object@genome$size)
    ))
    message(" - ", length(object@probes), " probes used")
    if (!is.null(exclude_regions)) {
      message("importing regions to exclude from analysis")
      if (class(exclude_regions) == "GRanges") {
        object@exclude <- GRanges(
          as.vector(seqnames(exclude_regions)),
          ranges(exclude_regions),
          seqinfo = Seqinfo(object@genome$chr,
                            object@genome$size)
        )
        values(object@exclude) <- values(exclude_regions)
        object@exclude <- sort(object@exclude)
      }
      else {
        object@exclude <- sort(rtracklayer::import(
          exclude_regions,
          seqinfo = Seqinfo(object@genome$chr, object@genome$size)
        ))
      }
    }
    else {
      object@exclude <- GRanges(seqinfo = Seqinfo(object@genome$chr,
                                                  object@genome$size))
    }
    if (!is.null(detail_regions)) {
      message("importing regions for detailed analysis")
      if (class(detail_regions) == "GRanges") {
        object@detail <- GRanges(
          as.vector(seqnames(detail_regions)),
          ranges(detail_regions),
          seqinfo = Seqinfo(object@genome$chr,
                            object@genome$size)
        )
        if (any(grepl("name", names(values(
          detail_regions
        ))))) {
          values(object@detail)$name <- values(detail_regions)[[grep("name",
                                                                     names(values(detail_regions)))[1]]]
        }
        if (any(grepl("IRanges", sapply(values(
          detail_regions
        ),
        class)))) {
          values(object@detail)$thick <-
            values(detail_regions)[[grep("IRanges",
                                         sapply(values(detail_regions), class))[1]]]
        }
        object@detail <- sort(object@detail)
      }
      else {
        object@detail <- sort(rtracklayer::import(
          detail_regions,
          seqinfo = Seqinfo(object@genome$chr, object@genome$size)
        ))
      }
      if (!is.element("name", names(values(object@detail)))) {
        stop("detailed region bed file must contain name column.")
      }
      if (!all(table(values(object@detail)$name) == 1)) {
        stop("detailed region names must be unique.")
      }
    }
    else {
      object@detail <- GRanges(seqinfo = Seqinfo(object@genome$chr,
                                                 object@genome$size))
    }
    if (!is.element("thick", names(values(object@detail)))) {
      values(object@detail)$thick <- resize(ranges(object@detail),
                                            fix = "center", 1e+06)
    }
    message("creating bins")
    anno.tile <-
      conumee:::CNV.create_bins(
        hg19.anno = object@genome,
        bin_minsize = bin_minsize,
        hg19.gap = object@gap,
        hg19.exclude = object@exclude
      )
    message(" - ", length(anno.tile), " bins created")
    message("merging bins")
    object@bins <- conumee:::CNV.merge_bins(
      hg19.anno = object@genome,
      hg19.tile = anno.tile,
      bin_minprobes = bin_minprobes,
      hg19.probes = object@probes,
      bin_maxsize = bin_maxsize
    )
    message(" - ", length(object@bins), " bins remaining")
    return(object)
  }




#' Get detail regions: yamat
#'
#' The detail regions are used in the authors' studies, focusing on pediatric
#' brain tumors and solid tumors.
#'
#' @param thick An integer scalar defined the flanking region size in bp.
#'   Default to 500000.
#' @return A \code{GRange} object.
#' @export
detail_regions <- function(thick = 500000) {
  e <- new.env()
  data(detail_regions, package = "conumee", envir = e)
  detail_regions <- e$detail_regions
  gene_list <-
    c("MYC",
      "MYCN",
      "CDKN2A/B",
      "PTCH1",
      "PTEN",
      "TP53",
      "CDK4",
      "CDK6",
      "ERBB2",
      "EGFR",
      "MET",
      "MDM2",
      "GLI2",
      "CCND1",
      "RB1",
      "NF1")
  keep <- mcols(detail_regions)$name %in% gene_list
  gr1 <- detail_regions[keep]
  GenomicRanges::mcols(gr1) <- GenomicRanges::mcols(gr1)[, c(1, 2)]
  df <- data.frame(
    chr = c("chr4", "chr19", "chr22", "chr12", "chr5", "chr1", "chr7", "chr6", "chr8", "chr17", "chr22"),
    start = c(55095256, 54169933, 29999545, 4382902, 1253282, 204485507, 140433812, 135502453, 67474410, 58677544, 24129118),
    end =   c(55164414, 54265683, 30094589, 4414522, 1295161, 204527248, 140624564, 135540311, 67525484, 58743640, 24176705),
    strand = "+",
    name = c("PDGFRA", "C19MC", "NF2", "CCND2", "TERT", "MDM4", "BRAF", "MYB", "MYBL1", "PPM1D", "SMARCB1"),
    stringsAsFactors = FALSE
  )
  gr2 <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
  GenomicRanges::mcols(gr2)$thick <- IRanges::IRanges(
    start = df$start - thick,
    end = df$end + thick
  )
  GenomeInfoDb::seqlevels(gr1) <- paste0("chr", c(seq(22), "X", "Y"))
  GenomeInfoDb::seqlevels(gr2) <- paste0("chr", c(seq(22), "X", "Y"))
  gr <- c(gr1, gr2)
  gr
}


#' Create CNV annotaion with \code{yamat} preset.
#'
#' The function is a wrapper of \code{\link{CNV.create_anno2}} with
#' pre-defined settings of detail regions.
#'
#' @param x An object of \code{\link[minfi]{RGChannelSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}}.
#' @param detail_regions_preset A character scalar of the name of preset detail
#'   regions, including "conumee" and "yamat". Default to "yamat".
#' @return An object of \code{\link[conumee]{CNV.anno-class}}.
#' @details Exclude regions are from \code{conumee} package.
#' @export
CNV.create_anno.yamat <-
  function(x,
           detail_regions_preset = c("yamat", "conumee")) {
    # Exclude regions
    e <- new.env()
    data(exclude_regions, package = "conumee", envir = e)
    # Include regions
    if (detail_regions_preset == "yamat") {
      detail_regions <- detail_regions()
    } else if (detail_regions_preset == "conumee") {
      data(detail_regions, package = "conumee", envir = e)
      detail_regions <- e$detail_regions
    } else {
      stop("detail_regions_preset should be either yamat or conumee.")
    }
    # Array names
    array_type <- minfi::annotation(x)["array"]
    if (array_type == "IlluminaHumanMethylation450k") {
      array_type <- "450k"
    } else if (array_type == "IlluminaHumanMethylationEPIC") {
      array_type <- "EPIC"
    } else {
      array_type <- "overlap"
    }
    # Create
    CNV.create_anno2(
      array_type = array_type,
      array_anno = minfi::annotation(x)["annotation"],
      chrXY = TRUE,
      exclude_regions = e$exclude_regions,
      detail_regions = detail_regions
    )
  }
