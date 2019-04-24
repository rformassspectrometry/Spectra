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
