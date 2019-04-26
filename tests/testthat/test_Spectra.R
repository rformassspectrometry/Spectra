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

test_that("msLevel,Spectra works", {
    sps <- Spectra()
    expect_equal(msLevel(sps), integer())

    df <- DataFrame(msLevel = c(1L, 2L))
    sps <- Spectra(df)
    expect_equal(msLevel(sps), c(1L, 2L))
})

test_that("removePeaks,Spectra works", {
    sps <- Spectra()
    res <- removePeaks(sps, t = 10)
    expect_true(length(res@processingQueue) == 1)
    expect_equal(res@processingQueue[[1]],
                 ProcessingStep(.remove_peaks,
                                list(t = 10, msLevel. = integer())))
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