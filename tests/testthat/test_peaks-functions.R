test_that(".peaks_remove works", {
    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = 1:length(int), intensity = int)
    res <- .peaks_remove(x, 1L, centroided = FALSE)
    expect_equal(res, x)
    res <- .peaks_remove(x, 1L, centroided = FALSE, t = 3)
    int_exp <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 0, 0, 0, 0,
                 1, 5, 10, 5, 1)
    expect_equal(res[, "intensity"], int_exp)

    res <- .peaks_remove(x, 1L, centroided = TRUE)
    int_exp <- c(0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 3, 10, 6, 2, 0, 0, 0, 2, 0, 0,
                 0, 5, 10, 5, 0)
    expect_equal(res[, "intensity"], int_exp)

    res <- .peaks_remove(x, 1L, centroided = TRUE, msLevel = 2L)
    expect_equal(res[, "intensity"], int)
})

test_that(".peaks_clean works", {
    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = 1:length(int), intensity = int)

    res <- .peaks_clean(x, 1L, all = FALSE)
    expect_true(is.matrix(res))
    expect_equal(res[, "intensity"], c(0, 1, 2, 3, 1, 0, 0, 1, 3, 10, 6, 2,
                                       1, 0, 1, 2, 0, 0, 1, 5, 10, 5, 1))
    res <- .peaks_clean(x, 1L, all = TRUE)
    expect_equal(res[, "intensity"], c(1, 2, 3, 1, 1, 3, 10, 6, 2, 1, 1, 2,
                                       1, 5, 10, 5, 1))

    expect_equal(.peaks_clean(x, 2L, msLevel = 1L), x)
    expect_equal(.peaks_clean(x, 1L, msLevel = 2L), x)
})

test_that(".peaks_pick works", {
    e <- matrix(NA_real_, nrow = 0, ncol = 2,
                dimnames = list(c(), c("mz", "intensity")))
    expect_warning(.peaks_pick(e, spectrumMsLevel = 1), "empty")
    expect_equal(suppressWarnings(.peaks_pick(e, spectrumMsLevel = 1)), e)

    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = 1:length(int), intensity = int)
    expect_equal(.peaks_pick(x, spectrumMsLevel = 2, msLevel = 1), x)
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, centroided = TRUE,
                             msLevel = 1), x)
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             SNR = 1),
                 cbind(mz = c(4, 12, 18, 23), intensity = c(3, 10, 2, 10)))
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             SNR = 5),
                 cbind(mz = c(12, 23), intensity = c(10, 10)))
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             SNR = 10), e)
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 10,
                             SNR = 1),
                 cbind(mz = c(12, 23), intensity = c(10, 10)))
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             method = "SuperSmoother", SNR = 2),
                 cbind(mz = c(4, 12, 23), intensity = c(3, 10, 10)))
})
