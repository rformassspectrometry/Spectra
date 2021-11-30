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
})
