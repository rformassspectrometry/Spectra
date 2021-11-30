test_that("validators and constructors for MsBackendCached work", {
    be <- new("MsBackendCached")
    expect_true(validObject(be))
    ## .valid_local_data
    df <- data.frame(a = 1:4, b = c("a", "b", "c", "d"))
    expect_equal(.valid_local_data(df, 4), NULL)
    expect_match(.valid_local_data(df, 3), "Number of rows")

    be@localData <- df
    expect_error(validObject(be), "Number of rows")
    be@nspectra <- 4L
    expect_true(validObject(be))

    be <- backendInitialize(be, nspectra = 10)
    expect_equal(nrow(be@localData), 10)

    expect_equal(length(be), 10)
    expect_equal(dataStorage(be), rep("<cache>", 10))
})

test_that("spectraVariables,MsBackendCached works", {
    be <- new("MsBackendCached")
    res <- spectraVariables(be)
    expect_equal(res, names(Spectra:::.SPECTRA_DATA_COLUMNS))

    df <- data.frame(a = 1:4, b = c("a", "b", "c", "d"))
    be@localData <- df
    be@nspectra <- 4L

    res <- spectraVariables(be)
    expect_equal(res, c(names(Spectra:::.SPECTRA_DATA_COLUMNS), "a", "b"))

    be@spectraVariables <- c("msLevel", "other_col")
    res <- spectraVariables(be)
    expect_equal(res, c(names(Spectra:::.SPECTRA_DATA_COLUMNS),
                        "a", "b", "other_col"))
})

test_that(".spectra_data MsBackendCached works", {
    be <- new("MsBackendCached")

    res <- .spectra_data(be)
    expect_true(nrow(res) == 0)
    expect_equal(sort(colnames(res)), sort(names(Spectra:::.SPECTRA_DATA_COLUMNS)))

    df <- data.frame(a = 1:4, msLevel = c(1L, 2L, 1L, 3L))
    be <- backendInitialize(be, data = df)
    be@localData <- df

    res <- .spectra_data(be)
    expect_equal(res$a, df$a)
    expect_equal(res$msLevel, df$msLevel)
    expect_true(all(is.na(res$rtime)))
    expect_true(all(lengths(res$mz) == 0))

    ## Just mz
    res <- .spectra_data(be, "mz")
    expect_true(is(res, "DataFrame"))
    expect_true(is(res$mz, "NumericList"))
    expect_true(all(lengths(res$mz) == 0))

    ## Just intensity
    res <- .spectra_data(be, "intensity")
    expect_true(is(res, "DataFrame"))
    expect_true(is(res$intensity, "NumericList"))
    expect_true(all(lengths(res$intensity) == 0))

    res <- .spectra_data(be, c("msLevel", "intensity"))
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), c("msLevel", "intensity"))
    expect_equal(res$msLevel, df$msLevel)
    expect_true(is(res$intensity, "NumericList"))
    expect_true(all(lengths(res$intensity) == 0))

    be@spectraVariables <- c("rtime", "mz", "intensity", "precursorMz")
    res <- .spectra_data(be)
    expect_true(!any(colnames(res) %in% c("rtime", "mz",
                                          "intensity", "precursorMz")))
    res <- .spectra_data(be, "mz")
    expect_true(is.null(res))
})

test_that("[,MsBackendCached works", {
    be <- MsBackendCached()

    be <- backendInitialize(be, nspectra = 10)
    res <- be[c(1, 4, 3), ]
    expect_true(length(res) == 3)
    expect_true(nrow(res@localData) == 3)

    df <- data.frame(msLevel = 1L, b = 1:6)
    be <- backendInitialize(be, data = df)
    res <- be[c(6, 1, 3)]
    expect_true(length(res) == 3)
    expect_equal(res@localData$b, c(6, 1, 3))
})

test_that("$<-,MsBackendCached works", {
    be <- MsBackendCached()

    expect_error(be$other_col <- 1:3, "value has to be")

    be <- backendInitialize(be, nspectra = 4)
    be$msLevel <- 2L
    expect_equal(spectraData(be, "msLevel")$msLevel, rep(2L, 4))

    be$rtime <- c(3.1, 5.2, 3.4, 6.2)
    expect_equal(spectraData(be)$rtime, c(3.1, 5.2, 3.4, 6.2))

    expect_error(be$rtime <- "Z", "wrong data type")
})

test_that("selectSpectraVariables,MsBackendCached works", {
    be <- MsBackendCached()
    df <- data.frame(msLevel = 1L, b = 1:6)
    be <- backendInitialize(
        be, data = df, spectraVariables = c("msLevel", "rtime", "precursorMz"))

    res <- selectSpectraVariables(be, spectraVariables = c("b", "dataOrigin"))
    expect_equal(colnames(res@localData), c("b"))
})
