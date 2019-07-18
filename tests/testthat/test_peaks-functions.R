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

test_that(".peaks_compare_intensities works", {
    x <- cbind(c(31.34, 50.14, 60.3, 120.9, 230, 514.13, 874.1),
               1:7)

    y <- cbind(c(12, 31.35, 70.3, x[4, 1] + x[4, 1] * 5 / 1e6,
                 230 + 230 * 10 / 1e6, 315, 514.14, 901, 1202),
               1:9)

    expect_true(is.na(.peaks_compare_intensities(x, y, ppm = 0)))

    res <- .peaks_compare_intensities(x, y, tolerance = 0.01)
    expect_equal(res, cor(c(1, 4, 5, 6), c(2, 4, 5, 7)))
    res <- .peaks_compare_intensities(x, y, tolerance = 0.01,
                                      method = "spearman")
    expect_equal(res, cor(c(1, 4, 5, 6), c(2, 4, 5, 7), method = "spearman"))

    ## ppm of 5, a single matching peak
    res <- .peaks_compare_intensities(x, y, ppm = 5,
                                      FUN = function(x, y) length(x))
    expect_identical(res, 1L)

    ## ppm of 10, two matching peaks
    res <- .peaks_compare_intensities(x, y, ppm = 10,
                                      FUN = function(x, y) length(x))
    expect_identical(res, 2L)
})
