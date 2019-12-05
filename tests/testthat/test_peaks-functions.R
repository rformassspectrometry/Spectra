test_that(".peaks_remove works", {
    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = 1:length(int), intensity = int)
    res <- .peaks_remove(x, 1L, centroided = FALSE)
    expect_equal(res, x)
    res <- .peaks_remove(x, 1L, centroided = FALSE, threshold = 3)
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

test_that(".peaks_bin works", {
    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = 1:length(int), intensity = int)

    res <- .peaks_bin(x, spectrumMsLevel = 1L)
    expect_identical(res[, 2], x[, 2])
    expect_identical(res[, 1], x[, 1] + 0.5)

    res <- .peaks_bin(x, spectrumMsLevel = 1L, binSize = 2L)
    expect_equal(res[, 2], c(1, 5, 1, 0, 1, 13, 8, 1, 3, 0, 6, 15, 1))
    res <- .peaks_bin(x, spectrumMsLevel = 1L, binSize = 2L, FUN = max)
    expect_equal(res[, 2], c(1, 3, 1, 0, 1, 10, 6, 1, 2, 0, 5, 10, 1))
    res <- .peaks_bin(x, spectrumMsLevel = 1L, msLevel = 2L, binSize = 2L)
    expect_identical(res, x)
})

test_that("joinPeaks works", {
    x <- cbind(c(31.34, 50.14, 60.3, 120.9, 230, 514.13, 874.1),
               1:7)
    y <- cbind(c(12, 31.35, 70.3, 120.9 + ppm(120.9, 5),
                 230 + ppm(230, 10), 315, 514.14, 901, 1202),
               1:9)

    res <- joinPeaks(x, y, ppm = 0)
    expect_true(nrow(res$x) == nrow(res$y))
    expect_true(nrow(res$x) == nrow(x) + nrow(y))
    res <- joinPeaks(x, y, ppm = 0, type = "inner")
    expect_true(nrow(res$x) == nrow(res$y))
    expect_true(nrow(res$x) == 0)

    ## ppm 5
    res <- joinPeaks(x, y, ppm = 5)
    expect_true(nrow(res$x) == nrow(res$y))
    expect_true(nrow(res$x) == nrow(x) + nrow(y) - 1)
    expect_true(!is.na(res$x[7, 1]))
    expect_true(!is.na(res$y[7, 1]))
    res <- joinPeaks(x, y, ppm = 5, type = "inner")
    expect_true(nrow(res$x) == nrow(res$y))

    ## ppm 10
    res <- joinPeaks(x, y, ppm = 10)
    expect_true(nrow(res$x) == nrow(res$y))
    expect_true(nrow(res$x) == nrow(x) + nrow(y) - 2)
    res <- joinPeaks(x, y, ppm = 10, type = "inner")
    expect_true(nrow(res$x) == nrow(res$y))
    expect_true(nrow(res$x) == 2)

    ## tolerance 0.01
    res <- joinPeaks(x, y, tolerance = 0.01)
    expect_true(nrow(res$x) == nrow(res$y))
    expect_true(nrow(res$x) == nrow(x) + nrow(y) - 4)
    res <- joinPeaks(x, y, tolerance = 0.01, type = "left")
    expect_true(nrow(res$x) == nrow(res$y))
    expect_equal(res$x, x)
    res <- joinPeaks(x, y, tolerance = 0.01, type = "right")
    expect_true(nrow(res$x) == nrow(res$y))
    expect_equal(res$y, y)
})

test_that(".peaks_pick works", {
    e <- matrix(NA_real_, nrow = 0, ncol = 2,
                dimnames = list(c(), c("mz", "intensity")))
    expect_warning(.peaks_pick(e, spectrumMsLevel = 1), "empty")
    expect_equal(suppressWarnings(.peaks_pick(e, spectrumMsLevel = 1)), e)

    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = seq_along(int), intensity = int)
    expect_equal(.peaks_pick(x, spectrumMsLevel = 2, msLevel = 1), x)
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, centroided = TRUE,
                             msLevel = 1), x)
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             snr = 1),
                 cbind(mz = c(4, 12, 18, 23), intensity = c(3, 10, 2, 10)))
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             snr = 5),
                 cbind(mz = c(12, 23), intensity = c(10, 10)))
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             snr = 10), e)
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 10,
                             snr = 1),
                 cbind(mz = c(12, 23), intensity = c(10, 10)))
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             method = "SuperSmoother", snr = 2),
                 cbind(mz = c(4, 12, 23), intensity = c(3, 10, 10)))

    ## use refineCentroids
    expect_equal(.peaks_pick(x, spectrumMsLevel = 1, halfWindowSize = 2,
                             method = "MAD", snr = 2, k = 1L, ),
                 cbind(mz = mapply(weighted.mean,
                                   x = list(3:5, 11:13, 22:24),
                                   w = list(int[3:5], int[11:13], int[22:24])),
                       intensity = c(3, 10, 10)))
})

test_that(".peaks_smooth works", {
    e <- matrix(NA_real_, nrow = 0, ncol = 2,
                dimnames = list(c(), c("mz", "intensity")))
    expect_warning(.peaks_smooth(e, spectrumMsLevel = 1), "empty")
    expect_equal(suppressWarnings(.peaks_smooth(e, spectrumMsLevel = 1)), e)

    cf <- matrix(0.2, nrow = 5, ncol = 5)
    int <- 1:10
    x <- cbind(mz = seq_along(int), intensity = int)
    expect_equal(.peaks_smooth(x, spectrumMsLevel = 1, coef = cf),
                 cbind(mz = x[, 1L], intensity = rep(3:8, c(3, 1, 1, 1, 1, 3))))
})

