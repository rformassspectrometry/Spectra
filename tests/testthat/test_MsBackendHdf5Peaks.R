test_df <- DataFrame(msLevel = c(1L, 2L, 2L), scanIndex = 4:6)
test_df$mz <- list(c(1.1, 1.3, 1.5), c(4.1, 5.1), c(1.6, 1.7, 1.8, 1.9))
test_df$intensity <- list(c(45.1, 34, 12), c(234.4, 1333), c(42.1, 34.2, 65, 6))

test_be <- backendInitialize(MsBackendHdf5Peaks(), file = tempfile(),
                             spectraData = test_df)

test_that("backendInitialize,MsBackendHdf5Peaks works", {
    if (!requireNamespace("rhdf5", quietly = TRUE))
        stop("Unable to load package rhdf5")
    res <- MsBackendHdf5Peaks()
    expect_true(is(res, "MsBackendHdf5Peaks"))
    expect_true(length(res) == 0)
    expect_true(validObject(res))

    ## Start with a file NA. need fromFile, need scanIndex.
    df <- DataFrame(msLevel = c(1L, 2L, 2L), scanIndex = 4:6)
    df$mz <- list(c(1.1, 1.3, 1.5), c(4.1, 5.1), c(1.6, 1.7, 1.8, 1.9))
    df$intensity <- list(c(45.1, 34, 12), c(234.4, 1333), c(42.1, 34.2, 65, 6))

    dr <- tempdir()
    res <- backendInitialize(MsBackendHdf5Peaks(), files = NA_character_,
                             spectraData = df, hdf5path = dr)
    expect_true(is(res, "MsBackendHdf5Peaks"))
    expect_identical(res$msLevel, c(1L, 2L, 2L))
    expect_equal(res@files, normalizePath(paste0(dr, "/spectra_peaks.h5")))
    expect_true(validObject(res))

    ## Two files.
    df$fromFile <- c(2L, 1L, 2L)
    res <- backendInitialize(MsBackendHdf5Peaks(), files = c("a.h5", "b.h5"),
                             spectraData = df, hdf5path = dr)
    expect_true(is(res, "MsBackendHdf5Peaks"))
    expect_equal(res@files, normalizePath(paste0(dr, c("/a.h5", "/b.h5"))))
    expect_identical(fromFile(res), c(2L, 1L, 2L))
    expect_true(validObject(res))

    ## Errors
    expect_error(backendInitialize(MsBackendHdf5Peaks(), files = NA_character_,
                                   spectraData = df, hdf5path = dr),
                 "do already exist")
    expect_error(backendInitialize(MsBackendHdf5Peaks(), spectraData = "a"),
                 "with spectrum data")
    expect_error(backendInitialize(MsBackendHdf5Peaks(), spectraData = df,
                                   files = c("a", "b", "c", "d")),
                 "Number of files does not match")

    ## Special cases
    res <- backendInitialize(MsBackendHdf5Peaks())
    expect_true(is(res, "MsBackendHdf5Peaks"))
    expect_true(length(res) == 0)
    expect_true(validObject(res))

    ## Empty m/z and/or intensity
    spd <- df
    spd$mz <- NULL
    spd$fromFile <- 1L
    expect_true(file.remove(paste0(dr, "/spectra_peaks.h5")))
    res <- backendInitialize(MsBackendHdf5Peaks(), files = NA_character_,
                             spectraData = spd, hdf5path = dr)
    expect_true(is(res, "MsBackendHdf5Peaks"))
    expect_true(validObject(res))
    expect_true(length(res) == 3)
})

test_that("intensity,MsBackendHdf5Peaks works", {
    be <- MsBackendHdf5Peaks()
    expect_equal(intensity(be), SimpleList())

    res <- intensity(test_be)
    expect_equal(res, SimpleList(test_df$intensity))
    expect_true(is.numeric(res[[1]]))
    expect_equal(length(res), length(test_be))
})

test_that("ionCount,MsBackendHdf5Peaks works", {
    be <- MsBackendHdf5Peaks()
    expect_equal(ionCount(be), numeric())

    res <- ionCount(test_be)
    expect_true(is.numeric(res))
    expect_true(length(res) == length(test_be))
    expect_identical(res, vapply(test_be$intensity, sum, numeric(1)))
})

test_that("isCentroided,MsBackendHdf5Peaks works", {
    be <- MsBackendHdf5Peaks()
    expect_equal(isCentroided(be), logical())

    res <- isCentroided(test_be)
    expect_true(all(is.na(res)))
    expect_true(is.logical(res))
    expect_true(length(res) == length(test_be))

    res <- isCentroided(sciex_hd5)
    expect_true(all(!res))
})

test_that("isEmpty,MsBackendHdf5Peaks works", {
    be <- MsBackendHdf5Peaks()
    expect_equal(isEmpty(be), logical())

    res <- isEmpty(test_be)
    expect_true(is.logical(res))
    expect_true(length(res) == length(test_be))
    expect_true(all(!res))

    df <- test_df
    df$intensity[[2]] <- numeric()
    df$mz[[2]] <- numeric()
    be <- backendInitialize(MsBackendHdf5Peaks(), files = tempfile(),
                            spectraData = df)
    res <- isEmpty(be)
    expect_identical(res, c(FALSE, TRUE, FALSE))
})

test_that("mz,MsBackendHdf5Peaks works", {
    be <- MsBackendHdf5Peaks()
    expect_equal(mz(be), SimpleList())

    res <- mz(test_be)
    expect_true(is(res, "SimpleList"))
    expect_true(is.numeric(res[[1]]))
    expect_identical(res, SimpleList(test_df$mz))
})

test_that("peaks,MsBackendHdf5Peaks works", {
    if (!requireNamespace("rhdf5", quietly = TRUE))
        stop("Unable to load package rhdf5")

    df <- DataFrame(msLevel = c(1L, 2L, 2L), scanIndex = 4:6)
    df$mz <- list(c(1.1, 1.3, 1.5), c(4.1, 5.1), c(1.6, 1.7, 1.8, 1.9))
    df$intensity <- list(c(45.1, 34, 12), c(234.4, 1333), c(42.1, 34.2, 65, 6))

    matl <- mapply(function(mz, intensity) cbind(mz, intensity),
                   df$mz, df$intensity)
    fl <- tempfile()
    be <- backendInitialize(MsBackendHdf5Peaks(), file = fl, spectraData = df)
    expect_identical(peaks(be), matl)

    df$fromFile <- c(2L, 1L, 2L)
    be <- backendInitialize(MsBackendHdf5Peaks(), spectraData = df,
                            files = c(tempfile(), tempfile()))
    expect_identical(peaks(be), matl)

    expect_identical(peaks(sciex_mzr), peaks(sciex_hd5))
})

test_that("peaksCount,MsBackendHdf5Peaks works", {
    be <- MsBackendHdf5Peaks()
    expect_equal(peaksCount(be), integer())

    res <- peaksCount(test_be)
    expect_true(is.integer(res))
    expect_true(length(res) == length(test_be))
    expect_identical(res, lengths(test_be$mz))
})

test_that("spectraData,MsBackendHdf5Peaks works", {
    be <- MsBackendHdf5Peaks()
    res <- spectraData(be)
    expect_true(nrow(res) == 0)
    expect_identical(colnames(res), spectraVariables(be))

    res <- spectraData(test_be)
    expect_identical(res$msLevel, test_df$msLevel)
    expect_identical(res$mz, SimpleList(test_df$mz))
    expect_identical(res$intensity, SimpleList(test_df$intensity))
})
