context("Pipeline based upon Conumee package")
library(yamatCN)
library(minfiData)

skip_flag <- TRUE

test_that("detail_regions()", {
  testthat::skip_if(skip_flag, "Skip")
  testthat::expect_error(detail_regions(), NA)
})

test_that("CNV.create_anno2(): EPIC array, ilm10b4 annotation.", {
  testthat::skip_if(skip_flag, "Skip")
  if (!skip_flag)
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  testthat::expect_error(CNV.create_anno2(), NA)
})

test_that("CNV.create_anno.yamat(): 450k with RGsetEx in minfiData package.", {
  testthat::skip_if(skip_flag, "Skip")
  if (!skip_flag)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  testthat::expect_error(CNV.create_anno.yamat(RGsetEx), NA)
})
