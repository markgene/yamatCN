context("CWOB Pipeline")
library(yamatCN)
library(minfiData)

skip_flag <- TRUE

ref <- RGsetEx[, 1:3]
qry <- RGsetEx[, 4:6]
report_dir <- tempdir()

test_that("cwob_pipe()", {
  testthat::skip_if(skip_flag, "Skip")
  testthat::expect_error(
    cwob_pipe(
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
