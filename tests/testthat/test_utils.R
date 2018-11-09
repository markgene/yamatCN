context("Utilities")
library(yamat)


test_that("pipe operator works", {
  testthat::expect_error(c(1, 2) %>% head(), NA)
})
