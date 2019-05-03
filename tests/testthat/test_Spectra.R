test_that("Spectra,DataFrame works", {
    df <- DataFrame()
    res <- Spectra(df)
    expect_true(validObject(res))
    expect_true(length(res) == 0)

    df <- DataFrame(msLevel = c(1L, 2L))
    res <- Spectra(df)
    expect_equal(msLevel(res), c(1L, 2L))
    expect_true(length(res) == 2)

    df$polarity <- "NEG"
    expect_error(Spectra(df), "wrong data type: polarity")
})

test_that("Spectra,missing works", {
    res <- Spectra()
    expect_true(length(res) == 0)

    be <- backendInitialize(MsBackendDataFrame(), files = NA_character_,
                            spectraData = DataFrame(msLevel = c(1L, 2L),
                                                    fromFile = 1L))
    res <- Spectra(backend = be)
    expect_true(length(res) == 2)
    expect_equal(msLevel(res), c(1L, 2L))
})

test_that("acquisitionNum,Spectra works", {
    sps <- Spectra()
    res <- acquisitionNum(sps)
    expect_equal(res, integer())
    df <- DataFrame(msLevel = c(1L, 1L))
    sps <- Spectra(df)
    res <- acquisitionNum(sps)
    expect_equal(res, c(NA_integer_, NA_integer_))
    df$acquisitionNum <- 1:2
    sps <- Spectra(df)
    res <- acquisitionNum(sps)
    expect_equal(res, 1:2)
})

test_that("fromFile,Spectra works", {
    sps <- Spectra()
    expect_equal(fromFile(sps), integer())

    df <- DataFrame(msLevel = c(1L, 2L))
    sps <- Spectra(df)
    expect_equal(fromFile(sps), c(1L, 1L))
})

test_that("length,Spectra works", {
    sps <- Spectra()
    expect_equal(length(sps), 0)

    df <- DataFrame(msLevel = c(1L, 2L))
    sps <- Spectra(df)
    expect_true(length(sps) == 2)
})

test_that("intensity,Spectra works", {
    sps <- Spectra()
    res <- intensity(sps)
    expect_true(is(res, "SimpleList"))
    expect_true(length(res) == 0)

    df <- DataFrame(msLevel = c(1L, 1L),
                    rtime = c(1.2, 1.4),
                    centroided = TRUE)
    df$intensity <- list(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                         c(9, 6, 0, 0, 3, 9, 4, 0, 0, 0))
    df$mz <- list(c(1:10), c(1:10))
    sps <- Spectra(df)
    res <- intensity(sps)
    expect_equal(res, SimpleList(
                          c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                          c(9, 6, 0, 0, 3, 9, 4, 0, 0, 0)))
    sps <- clean(sps, all = TRUE)
    res <- intensity(sps)
    expect_equal(res, SimpleList(c(1, 6, 3, 9, 1), c(9, 6, 3, 9, 4)))
})

test_that("msLevel,Spectra works", {
    sps <- Spectra()
    expect_equal(msLevel(sps), integer())

    df <- DataFrame(msLevel = c(1L, 2L))
    sps <- Spectra(df)
    expect_equal(msLevel(sps), c(1L, 2L))
})

test_that("mz,Spectra works", {
    sps <- Spectra()
    res <- mz(sps)
    expect_true(is(res, "SimpleList"))
    expect_true(length(res) == 0)

    df <- DataFrame(msLevel = c(1L, 1L),
                    rtime = c(1.2, 1.4),
                    centroided = TRUE)
    df$intensity <- list(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                         c(9, 6, 0, 0, 3, 9, 4, 0, 0, 0))
    df$mz <- list(c(1:10), c(1:10))
    sps <- Spectra(df)
    res <- mz(sps)
    expect_equal(res, SimpleList(1:10, 1:10))
    sps <- clean(sps, all = TRUE)
    res <- mz(sps)
    expect_equal(res, SimpleList(c(3, 4, 5, 8, 9), c(1, 2, 5, 6, 7)))
})

test_that("removePeaks,Spectra works", {
    sps <- Spectra()
    res <- removePeaks(sps, t = 10)
    expect_true(length(res@processingQueue) == 1)
    expect_equal(res@processingQueue[[1]],
                 ProcessingStep(.remove_peaks,
                                list(t = 10, msLevel. = integer())))
})

test_that("spectraData,Spectra works", {
    sps <- Spectra()
    res <- spectraData(sps)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)
    expect_true(all(c("msLevel", "rtime", "mz", "intensity") %in% colnames(res)))

    df <- DataFrame(msLevel = c(1L, 1L), rtime = c(1.2, 1.3), centroided = TRUE)
    df$mz <- list(1:10, 1:10)
    df$intensity <- list(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                         c(9, 6, 0, 0, 3, 0, 0, 0, 3, 2))
    sps <- Spectra(df)
    res <- spectraData(sps)
    expect_true(is(res, "DataFrame"))
    expect_equal(res$msLevel, c(1L, 1L))
    expect_equal(res$rtime, c(1.2, 1.3))
    expect_equal(res$precScanNum, c(NA_integer_, NA_integer_))
    expect_equal(res$mz, SimpleList(1:10, 1:10))
    expect_equal(res$intensity, SimpleList(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                                           c(9, 6, 0, 0, 3, 0, 0, 0, 3, 2)))
    sps <- clean(sps, all = TRUE)
    res <- spectraData(sps)
    expect_equal(res$mz, SimpleList(c(3, 4, 5, 8, 9), c(1, 2, 5, 9, 10)))
    expect_equal(res$intensity, SimpleList(c(1, 6, 3, 9, 1), c(9, 6, 3, 3, 2)))
})

test_that("spectraVariables,Spectra works", {
    sps <- Spectra()
    res <- spectraVariables(sps)
    exp_col <- c("msLevel", "rtime", "acquisitionNum", "scanIndex", "mz",
                 "intensity", "fromFile", "centroided", "smoothed",
                 "polarity", "precScanNum", "precursorMz", "precursorIntensity",
                 "precursorCharge", "collisionEnergy")
    expect_true(all(exp_col %in% res))
    df <- DataFrame(other_col = "a", msLevel = c(1L, 1L))
    sps <- Spectra(df)
    res <- spectraVariables(sps)
    expect_true(all(c(exp_col, "other_col") %in% res))
})

test_that("clean,Spectra works", {
    sps <- Spectra()
    res <- clean(sps)
    expect_true(length(res@processingQueue) == 1)
    expect_equal(res@processingQueue[[1]],
                 ProcessingStep(.clean_peaks,
                                list(all = FALSE, msLevel. = integer())))

    res <- clean(sps, all = TRUE, msLevel. = 2L)
    expect_true(length(res@processingQueue) == 1)
    expect_equal(res@processingQueue[[1]],
                 ProcessingStep(.clean_peaks,
                                list(all = TRUE, msLevel. = 2L)))
})

test_that("peaks,Spectra works", {
    df <- DataFrame(msLevel = c(1L, 2L), fromFile = 1L)
    df$mz <- list(1:4, 1:5)
    df$intensity <- list(1:4, 1:5)
    be <- backendInitialize(MsBackendDataFrame(), file = NA_character_, df)
    sps <- Spectra(backend = be)
    res <- peaks(sps)
    expect_true(is(res, "SimpleList"))
    expect_equal(res[[1]][, 1], 1:4)
    expect_equal(res[[2]][, 1], 1:5)

    sps <- Spectra(backend = MsBackendDataFrame())
    res <- peaks(sps)
    expect_true(is(res, "SimpleList"))
    expect_true(length(res) == 0)
})