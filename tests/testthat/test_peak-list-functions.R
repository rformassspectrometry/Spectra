test_that("combinePeaks works", {
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
    res <- combinePeaks(list(p1, p2, p3))
    p <- rbind(p1, p2, p3)
    p <- p[order(p[, 1]), ]
    expect_identical(p, res)

    ## With tolerance = 0.05 we group all "triplets"
    res <- combinePeaks(list(p1, p2, p3), tolerance = 0.05, mzFun = min,
                        intensityFun = min)
    expect_true(nrow(res) == length(mzs))

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

    res <- combinePeaks(list(p1, p2, p3), ppm = 10)
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

    res <- combinePeaks(list(p1, p2, p3), ppm = 10, intensityFun = median,
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
})
