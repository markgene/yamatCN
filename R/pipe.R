# Misc function of copy-number pipelines.

#' Check arguments of pipelines.
#'
#' @param report_dir A character scalar of reporting.
#' @return TRUE if pass all checks (invisible).
#' @noRd
.check_args_pipe <-
  function(report_dir) {
    if (missing(report_dir))
      stop("report_dir is required!")
    invisible(TRUE)
  }


#' Query sample names.
#'
#' @param x An object of \code{\link[minfi]{RGChannelSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}}.
#' @param label A character vector labeling query samples
#' @return A character vector of query sample names.
#' @noRd
.query_sample_names <- function(x, label = "query") {
  query_sample_names <-
    minfi::pData(x) %>%
    as.data.frame() %>%
    dplyr::filter(ref_query %in% label) %>%
    dplyr::select(Sample_Name) %>%
    unlist()
  names(query_sample_names) <- query_sample_names
  query_sample_names
}


#' Reference sample names.
#'
#' @param x An object of \code{\link[minfi]{RGChannelSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}}.
#' @param label A character vector labeling reference samples
#' @return A character vector of reference sample names.
#' @noRd
.reference_sample_names <- function(x, label = "ref") {
  ref_sample_names <-
    minfi::pData(x) %>%
    as.data.frame() %>%
    dplyr::filter(ref_query %in% label) %>%
    dplyr::select(Sample_Name) %>%
    unlist()
  names(ref_sample_names) <- ref_sample_names
  ref_sample_names
}


#' Query sample indices.
#'
#' @param x An object of \code{\link[minfi]{RGChannelSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}}.
#' @param label A character vector labeling query samples.
#' @return An integer vector of query sample indices.
#' @noRd
.query_indices <- function(x, label = "query") {
  df <- minfi::pData(x) %>%
    as.data.frame()
  which(df$ref_query == label)
}


#' Reference sample names.
#'
#' @param x An object of \code{\link[minfi]{RGChannelSet-class}} or
#'   \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}}.
#' @param label A character vector labeling reference samples
#' @return An integer vector of reference sample indices.
#' @noRd
.reference_indices <- function(x, label = "ref") {
  df <- minfi::pData(x) %>%
    as.data.frame()
  which(df$ref_query == label)
}



#' Calculate log2 ratio of observed and expected values of intensities.
#'
#' Calculate log2 ratio of observed and expected values of intensities at each
#' probe.
#'
#' @param query A matrix of query sample intensity.
#' @param reference A matrix of reference sample intensity.
#' @param method A character scalar. Default to "mean".
#' @return A matrix of log2 ratio of query samples.
#' @note If an intensity is less than 1, it is set to 1.
#' @export
cal_log2_ratio <- function(query, reference, method = c("mean", "median", "lm")) {
  method <- match.arg(method)
  query[query < 1] <- 1
  reference[reference < 1] <- 1
  r <- cor(reference, query) < 0.99
  if (any(!r)) warning("query sample seems to also be in the reference set. not used for fit.")
  switch(
    method,
    mean = cal_log2_ratio.mean(query, reference),
    median = cal_log2_ratio.median(query, reference),
    lm = cal_log2_ratio.lm(query, reference))
}


#' Calculate log2 ratio of observed and expected values of intensities.
#'
#' Calculate log2 ratio of observed and expected values of intensities at each
#' probe. Expected values are calculated at the mean of intensities in reference
#' dataset.
#'
#' @param query A matrix of query sample intensity.
#' @param reference A matrix of reference sample intensity.
#' @return A matrix of log2 ratio of query samples.
#' @noRd
cal_log2_ratio.mean <- function(query, reference) {
  log2(query/rowMeans(reference))
}


#' Calculate log2 ratio of observed and expected values of intensities.
#'
#' Calculate log2 ratio of observed and expected values of intensities at each
#' probe. Expected values are calculated at the median of intensities in
#' reference dataset.
#'
#' @param query A matrix of query sample intensity.
#' @param reference A matrix of reference sample intensity.
#' @return A matrix of log2 ratio of query samples.
#' @noRd
cal_log2_ratio.median <- function(query, reference) {
  log2(query/rowMedians(reference))
}


#' Calculate log2 ratio of observed and expected values of intensities.
#'
#' Calculate log2 ratio of observed and expected values of intensities at each
#' probe. Expected values are calculated with linear regression model.
#'
#' @param query A matrix of query sample intensity.
#' @param reference A matrix of reference sample intensity.
#' @return A matrix of log2 ratio of query samples.
#' @noRd
#' @note If the predicted value is less than 1, set it to 1.
cal_log2_ratio.lm <- function(query, reference) {
  r <- cor(reference, query) < 0.99
  log_ratio_qry <- lapply(
    seq_len(ncol(query)),
    function(i) {
      q <- query[, i]
      mod <- lm(q ~ ., data = data.frame(q = q, R = reference[, r[, i]]))
      pred <- as.data.frame(predict(mod, interval = "confidence"))
      pred$fit[pred$fit < 1] <- 1
      log2(q / pred$fit)
    }
  ) %>% do.call(cbind, .)
  colnames(log_ratio_qry) <- colnames(query)
  log_ratio_qry
}


#' Get influential (i.e. outlier) loci
#'
#' @param x An object of \code{\link[minfi]{MethylSet-class}} or
#'   \code{\link[minfi]{GenomicMethylSet-class}}.
#' @param threshold A numeric scalar to define outlier. It is the times
#'   that Cook's distance greater than the mean.
#' @param verbose A logical scalar.
#' @return A character vector of loci names.
get_influential_loci <- function(x, threshold = 4, verbose = TRUE) {
  intensity <- as.data.frame(copy_number(x))
  influential_list <- lapply(
    seq(ncol(intensity)),
    function(i) {
      if (verbose) message("Process ", i, " out of ", ncol(intensity), " samples.")
      mod <- lm(y ~ ., data = data.frame(y = intensity[, i], X = intensity[, -i]))
      cooksd <- cooks.distance(mod)
      rownames(intensity)[which(cooksd > threshold * mean(cooksd, na.rm=T))]
    })
  unique(do.call(c, influential_list))
}


#' Calculate LRR shift.
#'
#' LRR are shifted to minimize the median absolute deviation from all markers to
#' zero to determine the copy-number neutral state. The function is a wrapper
#' of \code{\link[stats]{optim}}. See the manual of \code{\link[stats]{optim}}
#' for more information about the arguments.
#'
#' @param lrr A matrix of LRR.
#' @param method A character scalar of the method to be used.
#' @param lower A numeric scalar. Default to -100.
#' @param upper A numeric scalar. Default to 100.
#' @param ... Arguments passed to \code{\link[stats]{optim}}.
#' @return A numeric scalar of shift.
#' @note I borrow the code from \code{\link[conumee]{CNV.bin}} package.
#' @noRd
.cal_lrr_shift <- function(lrr, method = "Brent", lower = -100, upper = 100, ...) {
  optim(
    par = 0,
    fn = function(s)
      median(abs(lrr - s), na.rm = TRUE),
    method = method,
    lower = lower,
    upper = upper,
    ...
  )$par
}
