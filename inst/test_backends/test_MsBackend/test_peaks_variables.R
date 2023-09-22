#' Tests for peaks variables.

test_that("peaksVariables returns at least mz and intensity", {
    expect_true(all(c("mz", "intensity") %in% peaksVariables(be)))
})

test_that("peaksData returns a matrix or data.frame", {
    res <- peaksData(be)
    expect_equal(length(res), length(be))
    ok <- vapply(res, function(z) is.matrix(z) | is.data.frame(z), logical(1))
    expect_true(all(ok))
})

test_that("peaksData columns and peaksVariables match", {
    res <- peaksData(be, peaksVariables(be))
    expect_equal(colnames(res[[1L]]), peaksVariables(be))
})
