context("Pipeline based upon Conumee package")
library(yamatCN)
library(minfiData)

ref <- RGsetEx[, 1:3]
qry <- RGsetEx[, 4:6]
report_dir <- tempdir()

test_that(".check_cn_pipe_conumee()", {
  testthat::expect_error(.check_cn_pipe_conumee(
    ref = ref,
    qry = qry,
    report_dir = report_dir,
    norm_method = "swan",
    batch = NULL,
    batch2 = NULL
  ),
  NA)
})
