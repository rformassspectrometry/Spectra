test_that("combinePeaksData works", {
    set.seed(123)
    mzs <- seq(1, 20, 0.1)
    ints1 <- abs(rnorm(length(mzs), 10))
    ints1[11:20] <- c(15, 30, 90, 200, 500, 300, 100, 70, 40, 20) # add peak
    ints2 <- abs(rnorm(length(mzs), 10))
    ints2[11:20] <- c(15, 30, 60, 120, 300, 200, 90, 60, 30, 23)
    ints3 <- abs(rnorm(length(mzs), 10))
    ints3[11:20] <- c(13, 20, 50, 100, 200, 100, 80, 40, 30, 20)
    ## Create the peaks matrices
    p1 <- cbind(mz = mzs + rnorm(length(mzs), sd = 0.01),
                intensity = ints1)
    p2 <- cbind(mz = mzs + rnorm(length(mzs), sd = 0.01),
                intensity = ints2)
    p3 <- cbind(mz = mzs + rnorm(length(mzs), sd = 0.009),
                intensity = ints3)

    ## With tolerance = 0 and ppm = 0 we expect to get just a joined matrix.
    expect_warning(combinePeaks(list(p1, p2, p3)), "deprecated")
    res <- combinePeaksData(list(p1, p2, p3))
    p <- rbind(p1, p2, p3)
    p <- p[order(p[, 1]), ]
    expect_identical(p, res)

    ## With tolerance = 0.05 we group all "triplets"
    res <- combinePeaksData(list(p1, p2, p3), tolerance = 0.05, mzFun = min,
                        intensityFun = min)
    expect_true(nrow(res) == length(mzs))
    expect_warning(res2 <- combinePeaks(list(p1, p2, p3), tolerance = 0.05,
                                        mzFun = min, intensityFun = min),
                   "deprecated")
    expect_equal(res, res2)

    p <- do.call(rbind, lapply(1:nrow(res), function(i) {
        c(mz = min(p1[i, 1], p2[i, 1], p3[i, 1]),
          intensity = min(p1[i, 2], p2[i, 2], p3[i, 2]))
    }))
    expect_identical(res, p)

    ## Other example.
    p1 <- cbind(mz = c(23, 24, 25, 26), intensity = c(10, 20, 30, 40))
    p2 <- cbind(mz = c(21, 23 + ppm(23, 10), 23.5, 26 - ppm(26, 10)),
                intensity = c(12, 22, 32, 42))
    p3 <- cbind(mz = c(23 - ppm(23, 5), 24), intensity = c(13, 23))

    res <- combinePeaksData(list(p1, p2, p3), ppm = 10)
    expect_true(nrow(res) == 6)
    expect_true(res[1, 2] == 12)
    expect_true(res[2, 2] == mean(c(10, 22, 13)))
    expect_true(res[3, 2] == mean(c(32)))
    expect_true(res[4, 2] == mean(c(20, 23)))
    expect_true(res[5, 2] == mean(c(30)))
    expect_true(res[6, 2] == mean(c(40, 42)))
    ## m/z
    expect_true(res[1, 1] == p2[1, 1])
    expect_true(res[2, 1] == mean(c(p1[1, 1], p2[2, 1], p3[1, 1])))
    expect_true(res[3, 1] == p2[3, 1])
    expect_true(res[4, 1] == mean(c(p1[2, 1], p3[2, 1])))
    expect_true(res[5, 1] == mean(p1[3, 1]))
    expect_true(res[6, 1] == mean(c(p1[4, 1], p2[4, 1])))

    res <- combinePeaksData(list(p1, p2, p3), ppm = 10, intensityFun = median,
                        mzFun = max)
    expect_true(nrow(res) == 6)
    expect_true(res[1, 2] == 12)
    expect_true(res[2, 2] == median(c(10, 22, 13)))
    expect_true(res[3, 2] == median(c(32)))
    expect_true(res[4, 2] == median(c(20, 23)))
    expect_true(res[5, 2] == median(c(30)))
    expect_true(res[6, 2] == median(c(40, 42)))
    ## m/z
    expect_true(res[1, 1] == p2[1, 1])
    expect_true(res[2, 1] == max(c(p1[1, 1], p2[2, 1], p3[1, 1])))
    expect_true(res[3, 1] == p2[3, 1])
    expect_true(res[4, 1] == max(c(p1[2, 1], p3[2, 1])))
    expect_true(res[5, 1] == max(p1[3, 1]))
    expect_true(res[6, 1] == max(c(p1[4, 1], p2[4, 1])))

    ## peaks = "intersect"
    p1 <- cbind(mz = c(12, 45, 64, 70), intensity = c(10, 20, 30, 40))
    p2 <- cbind(mz = c(17, 45.1, 63.9, 70.2), intensity = c(11, 21, 31, 41))
    p3 <- cbind(mz = c(12.1, 44.9, 63), intensity = c(12, 22, 32))

    expect_identical(combinePeaksData(list(p1), peaks = "intersect"), p1)

    res <- combinePeaksData(list(p1, p2, p3), peaks = "intersect",
                            minProp = 0.2)
    expect_equal(res[, 1], sort(c(p1[, 1], p2[, 1], p3[, 1])))
    idx <- order(c(p1[, 1], p2[, 1], p3[, 1]))
    expect_equal(res[, 2], c(p1[, 2], p2[, 2], p3[, 2])[idx])

    res_2 <- combinePeaksData(list(p1, p2, p3), peaks = "union", main = 3,
                              tolerance = 1)
    expect_true(nrow(res_2) == 3)
    expect_equal(unname(res_2[1, 1]), mean(c(12, 12.1)))
    expect_equal(unname(res_2[2, 1]), mean(c(45, 45.1, 44.9)))
    expect_equal(unname(res_2[3, 1]), mean(c(64, 63.9, 63)))
    expect_equal(unname(res_2[1, 2]), mean(c(10, 12)))
    expect_equal(unname(res_2[2, 2]), mean(c(20, 21, 22)))
    expect_equal(unname(res_2[3, 2]), mean(c(30, 31, 32)))

    res <- combinePeaksData(list(p1, p2, p3), peaks = "intersect")
    expect_true(nrow(res) == 0)

    res <- combinePeaksData(list(p1, p2, p3), tolerance = 0.1,
                            peaks = "intersect")
    expect_equal(res[, 1], c(mean(c(12, 12.1)), mean(c(45, 45.1, 44.9)),
                             mean(c(64, 63.9))))
    expect_equal(res[, 2], c(mean(c(10, 12)), mean(c(20, 21, 22)),
                             mean(c(30, 31))))

    res <- combinePeaksData(list(p1, p2, p3), tolerance = 0.1, mzFun = min,
                            intensityFun = max, peaks = "intersect")
    expect_equal(res[, 1], c(min(c(12, 12.1)), min(c(45, 45.1, 44.9)),
                             min(c(64, 63.9))))
    expect_equal(res[, 2], c(max(c(10, 12)), max(c(20, 21, 22)),
                             max(c(30, 31))))

    expect_error(combinePeaksData(list(p1, p2, p3), main = 5), "has to be")
})

test_that(".peaks_compare works", {
    sps <- Spectra(sciex_hd5)[120:126]
    sps <- setBackend(sps, MsBackendDataFrame())

    res <- .peaks_compare(.peaksapply(sps), .peaksapply(sps))
    expect_true(ncol(res) == length(sps))
    expect_true(nrow(res) == length(sps))
    expect_equal(diag(res), rep(1, length(sps)))

    res_2 <- .peaks_compare(.peaksapply(sps), .peaksapply(sps[3]))
    expect_true(ncol(res_2) == 1)
    expect_true(nrow(res_2) == length(sps))
    expect_identical(res_2[, 1], res[, 3])

    res_2 <- .peaks_compare(.peaksapply(sps[5]), .peaksapply(sps))
    expect_true(ncol(res_2) == length(sps))
    expect_true(nrow(res_2) == 1)
    expect_identical(res_2[1, ], res[5, ])
})
