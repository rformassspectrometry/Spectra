test_that("PrecursorMzParam works", {
    expect_error(PrecursorMzParam("other"), "should be one")
    expect_error(PrecursorMzParam(filterPeaks = 12), "should be either")
    expect_error(PrecursorMzParam(msLevel = "one"), "type integer")

    res <- PrecursorMzParam()
    expect_s4_class(res, "PrecursorMzParam")
    expect_true(is.function(res@filterPeaks))
    X <- cbind(mz = 1:5, intensity = c(100, 200, 300, 400, 500))
    expect_equal(res@filterPeaks(X), X)

    res <- PrecursorMzParam("abovePrecursor")
    expect_equal(res@filterPeaks(X, 3), X[X[, 1] < 3, ])

    res <- PrecursorMzParam("belowPrecursor")
    expect_equal(res@filterPeaks(X, 3), X[X[, 1] > 3, ])
})

test_that("neutralLoss,Spectra,PrecursorMzParam works", {
    DF <- DataFrame(msLevel = c(1L, 2L, 3L, 1L, 2L, 3L),
                    precursorMz = c(NA, 40, 20, NA, 300, 200))
    DF$mz <- IRanges::NumericList(
                          c(3, 12, 14, 15, 16, 200),
                          c(13, 23, 39, 86),
                          c(5, 7, 20, 34, 50),
                          c(5, 7, 9, 20, 100),
                          c(15, 53, 299, 300),
                          c(34, 56, 100, 200, 204, 309)
                      , compress = FALSE)
    DF$intensity <- IRanges::NumericList(1:6, 1:4, 1:5, 1:5, 1:4, 1:6,
                                         compress = FALSE)

    sps <- Spectra(DF, backend = MsBackendDataFrame())

    prm <- PrecursorMzParam(filterPeaks = "none", msLevel = 2L)
    res <- neutralLoss(sps, prm)
    expect_equal(res@metadata[[1L]], prm)
    expect_equal(mz(res)[[1L]], mz(sps)[[1L]])
    expect_equal(mz(res)[[3L]], mz(sps)[[3L]])
    expect_equal(mz(res)[[4L]], mz(sps)[[4L]])
    expect_equal(mz(res)[[6L]], mz(sps)[[6L]])
    expect_equal(mz(res)[[2L]], sort(40 - mz(sps)[[2L]]))
    expect_equal(mz(res)[[5L]], sort(300 - mz(sps)[[5L]]))

    res <- neutralLoss(sps, PrecursorMzParam(filterPeaks = "abovePrecursor",
                                             msLevel = 2:3))
    expect_equal(mz(res)[[1L]], mz(sps)[[1L]])
    expect_equal(mz(res)[[4L]], mz(sps)[[4L]])
    expect_equal(mz(res)[[2L]], c(1, 17, 27))
    expect_equal(mz(res)[[3L]], c(13, 15))
    expect_equal(mz(res)[[5L]], c(1, 247, 285))
    expect_equal(mz(res)[[6L]], c(100, 144, 166))

    res <- neutralLoss(sps, PrecursorMzParam(filterPeaks = "belowPrecursor",
                                             msLevel = 3))
    expect_equal(mz(res)[[1L]], mz(sps)[[1L]])
    expect_equal(mz(res)[[2L]], mz(sps)[[2L]])
    expect_equal(mz(res)[[4L]], mz(sps)[[4L]])
    expect_equal(mz(res)[[5L]], mz(sps)[[5L]])

    expect_equal(mz(res)[[3L]], c(-30, -14))
    expect_equal(mz(res)[[6L]], c(-109, -4))

    res <- neutralLoss(sps, PrecursorMzParam(filterPeaks = "belowPrecursor",
                                             msLevel = 2))
    expect_equal(unname(mz(res)[[2L]]), c(-46))
    expect_equal(mz(res)[[5L]], numeric())

    ## With and without NA MS level.
    a <- neutralLoss(sps[2], PrecursorMzParam())
    sps_2 <- sps
    sps_2$msLevel <- NA_integer_
    b <- neutralLoss(sps_2[2], PrecursorMzParam())
    expect_equal(mz(a), mz(b))
    b <- neutralLoss(sps, PrecursorMzParam(msLevel = 10))
    expect_equal(mz(sps), mz(b))

    ## With precursor m/z being NA.
    a <- sps[2]
    a$precursorMz <- NA_integer_
    a <- neutralLoss(a, PrecursorMzParam())
    expect_equal(mz(a)[[1L]], numeric())
})
