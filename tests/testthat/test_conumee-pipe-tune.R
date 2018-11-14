context("Test tuning Conumee pipeline")
library(yamatCN)

skip_flag <- FALSE

library(minfiData)
ref <- RGsetEx[, 1:3]
qry <- RGsetEx[, 4:6]
report_dir <- tempdir()

test_that("conumee_pipe()", {
  testthat::skip_if(skip_flag, "Skip")
  if (!skip_flag)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  testthat::expect_error(
    tune_conumee_pipe(
      ref = ref,
      qry = qry,
      report_dir = report_dir,
      norm_methods = c("swan", "methylcnv"),
      batch = NULL,
      batch2 = NULL
    ),
    NA
  )
})
