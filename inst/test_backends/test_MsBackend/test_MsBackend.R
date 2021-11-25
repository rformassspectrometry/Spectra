#' This test suite runs tests for the most important methods defined in
#' `MsBackend.R` on a provided `MsBackend` instance and checks if the returned
#' data/results are in the expected format. This unit tests might be run for
#' a newly developed backend extending `MsBackend` to ensure that it is
#' compliant with it.
#'
#' To run this unit tests from another package:
#'
#' be <- <code to create and initialize the backend>
#' test_suite <- system.file("test_backends", "test_MsBackend",
#'     package = "Spectra")
#' test_dir(test_suite)
#'
#' The unit tests in this suite expect a variable `be` to be defined, which
#' has to represent an **already initialized** backend instance.

test_that("backend is valid", {
    validObject(be)
})

#' backendMerge - might not be possible for all...
#' export - might not be possible for all

test_that("acquisitionNum", {
    res <- acquisitionNum(be)
    expect_true(is.integer(res))
    expect_true(length(res) == length(be))
})

test_that("peaksData", {
})

test_that("centroided", {
})
