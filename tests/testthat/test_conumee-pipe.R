context("Pipeline based upon Conumee package")
library(yamatCN)
library(minfiData)

skip_flag <- FALSE

ref <- RGsetEx[, 1:3]
qry <- RGsetEx[, 4:6]
report_dir <- tempdir()

test_that("conumee_pipe()", {
  testthat::skip_if(skip_flag, "Skip")
  if (!skip_flag)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  testthat::expect_error(
    conumee_pipe(
      ref = ref,
      qry = qry,
      report_dir = report_dir,
      norm_method = "swan",
      batch = NULL,
      batch2 = NULL
    ),
    NA
  )
})
