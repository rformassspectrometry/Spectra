test_df <- DataFrame(msLevel = c(1L, 2L, 2L))
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

    ## Start with missing file, no scanIndex
    df <- DataFrame(msLevel = c(1L, 2L, 2L))
    df$mz <- list(c(1.1, 1.3, 1.5), c(4.1, 5.1), c(1.6, 1.7, 1.8, 1.9))
    df$intensity <- list(c(45.1, 34, 12), c(234.4, 1333), c(42.1, 34.2, 65, 6))
    df$dataStorage <- "myfile"

    dr <- tempdir()
    res <- backendInitialize(MsBackendHdf5Peaks(), spectraData = df,
                             hdf5path = dr)
    expect_true(is(res, "MsBackendHdf5Peaks"))
    expect_identical(res$msLevel, c(1L, 2L, 2L))
    expect_equal(dataStorageLevels(res), normalizePath(paste0(dr, "/myfile.h5")))
    expect_true(validObject(res))
    expect_identical(mz(res), as(df$mz, "NumericList"))
    expect_identical(scanIndex(res), 1:3)

    ## Two dataStorage, two file names.
    df$dataStorage <- c("2", "1", "2")
    res <- backendInitialize(MsBackendHdf5Peaks(), files = c("a.h5", "b.h5"),
                             spectraData = df, hdf5path = dr)
    expect_true(is(res, "MsBackendHdf5Peaks"))
    expect_equal(dataStorageLevels(res),
                 normalizePath(paste0(dr, c("/a.h5", "/b.h5"))))
    expect_true(validObject(res))
    expect_identical(scanIndex(res), 1:3)

    ## Two dataStorage, no files provided
    res <- backendInitialize(MsBackendHdf5Peaks(),
                             spectraData = df, hdf5path = dr)
    expect_true(is(res, "MsBackendHdf5Peaks"))
    expect_equal(dataStorageLevels(res),
                 normalizePath(paste0(dr, c("/2.h5", "/1.h5"))))
    expect_true(validObject(res))
    expect_identical(scanIndex(res), 1:3)

    ## Same but with scanIndex provide
    df$scanIndex <- c(1L, 1L, 2L)
    res <- backendInitialize(MsBackendHdf5Peaks(), files = c("c", "d"),
                             spectraData = df, hdf5path = dr)
    expect_true(is(res, "MsBackendHdf5Peaks"))
    expect_equal(dataStorageLevels(res),
                 normalizePath(paste0(dr, c("/c", "/d"))))
    expect_true(validObject(res))
    expect_identical(scanIndex(res), c(1L, 1L, 2L))

    ## Errors
    expect_error(backendInitialize(MsBackendHdf5Peaks(), files = "single",
                                   spectraData = df, hdf5path = dr),
                 "provided file names has to match")
    expect_error(backendInitialize(MsBackendHdf5Peaks(), spectraData = "a"),
                 "with spectrum data")
    expect_error(backendInitialize(MsBackendHdf5Peaks(), spectraData = df,
                                   files = c("a", "b", "c", "d")),
                 "Number of provided")
    df$dataStorage <- NULL
    expect_error(backendInitialize(MsBackendHdf5Peaks(), spectraData = df,
                                   hdf5path = dr),
                 "Column \"dataStorage\" is required")

    ## Special cases
    res <- backendInitialize(MsBackendHdf5Peaks())
    expect_true(is(res, "MsBackendHdf5Peaks"))
    expect_true(length(res) == 0)
    expect_true(validObject(res))

    ## Empty m/z and/or intensity
    spd <- df
    spd$mz <- NULL
    spd$scanIndex <- 1:3
    res <- backendInitialize(MsBackendHdf5Peaks(), files = "no-peak-data",
                             spectraData = spd, hdf5path = dr)
    expect_true(is(res, "MsBackendHdf5Peaks"))
    expect_true(validObject(res))
    expect_true(length(res) == 3)
})

test_that("intensity,MsBackendHdf5Peaks works", {
    be <- MsBackendHdf5Peaks()
    expect_identical(intensity(be), NumericList(compress = FALSE))

    res <- intensity(test_be)
    expect_identical(res, NumericList(test_df$intensity, compress = FALSE))
    expect_true(is.numeric(res[[1]]))
    expect_equal(length(res), length(test_be))
})

test_that("intensity<-,MsBackendHdf5Peaks works", {
    be <- backendInitialize(MsBackendHdf5Peaks(), files = tempfile(), test_df)

    ints <- lapply(test_df$intensity, function(z) z/2)
    intensity(be) <- ints
    expect_identical(intensity(be), NumericList(ints, compress = FALSE))
    expect_identical(be@modCount, 1L)

    expect_error(intensity(be) <- 4, "has to be a list")
    expect_error(intensity(be) <- list(4), "has to match the length")
    expect_error(intensity(be) <- list(4, 3, 5), "match the number of peaks")
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
    expect_identical(mz(be), NumericList(compress = FALSE))

    res <- mz(test_be)
    expect_true(is(res, "NumericList"))
    expect_true(is.numeric(res[[1]]))
    expect_identical(res, NumericList(test_df$mz, compress = FALSE))
})

test_that("mz<-,MsBackendHdf5Peaks works", {
    be <- backendInitialize(MsBackendHdf5Peaks(), files = tempfile(), test_df)

    mzs <- lapply(test_df$mz, function(z) z/2)
    mz(be) <- mzs
    expect_identical(mz(be), NumericList(mzs, compress = FALSE))
    expect_identical(be@modCount, 1L)

    expect_error(mz(be) <- 4, "has to be a list")
    expect_error(mz(be) <- list(4), "has to match the length")
    expect_error(mz(be) <- list(4, 3, 5), "match the number of peaks")
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

    df$dataStorage <- c("2", "1", "2")
    be <- backendInitialize(MsBackendHdf5Peaks(), spectraData = df,
                            files = c(tempfile(), tempfile()))
    expect_identical(peaks(be), matl)

    expect_identical(peaks(sciex_mzr), peaks(sciex_hd5))
})

test_that("peaks<-,MsBackendHdf5Peaks works", {
    be <- backendInitialize(MsBackendHdf5Peaks(), files = tempfile(),
                            spectraData = test_df)
    pks <- peaks(be)
    pks_2 <- list(pks[[1]][2, , drop = FALSE],
                  pks[[2]],
                  pks[[3]][1:3, ])
    peaks(be) <- pks_2
    expect_identical(be@modCount, 1L)
    expect_identical(peaks(be), pks_2)
    pks_2[[2]] <- pks_2[[2]][0, ]

    peaks(be) <- pks_2
    expect_identical(be@modCount, 2L)
    expect_identical(peaks(be), pks_2)
    expect_identical(peaksCount(be), c(1L, 0L, 3L))
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
    expect_identical(res$mz, NumericList(test_df$mz, compress = FALSE))
    expect_identical(res$intensity, NumericList(test_df$intensity,
                                                compress = FALSE))
})

test_that("spectraData<-,MsBackendHdf5Peaks works", {
    be <- backendInitialize(MsBackendHdf5Peaks(), files = tempfile(),
                            spectraData = test_df)
    df <- test_df
    df$new_col <- "a"
    df$msLevel <- 1L
    df$intensity <- NULL
    df$mz <- NULL

    spectraData(be) <- df
    expect_identical(be@modCount, 0L)
    expect_identical(msLevel(be), rep(1L, 3))
    expect_identical(intensity(be), NumericList(test_df$intensity,
                                                compress = FALSE))
    expect_identical(be$new_col, c("a", "a", "a"))
    expect_identical(scanIndex(be), 1:3)

    ## Only m/z, no intensities.
    df$mz <- test_df$mz
    spectraData(be) <- df
    expect_identical(mz(be), NumericList(test_df$mz, compress = FALSE))
    expect_true(all(is.na(unlist(intensity(be)))))
    expect_identical(be@modCount, 1L)

    ## Update the m/z and intensities.
    df$mz <- list(c(1.1, 1.2), numeric(), c(1.3, 1.4))
    df$intensity <- list(c(45.3, 345.1), numeric(), c(1234.2, 12.1))

    spectraData(be) <- df
    expect_identical(intensity(be), NumericList(df$intensity, compress = FALSE))
    expect_identical(mz(be), NumericList(df$mz, compress = FALSE))
    expect_identical(peaksCount(be), c(2L, 0L, 2L))
    expect_identical(be@modCount, 2L)

    ## Error:
    expect_error(spectraData(be) <- "a", "has to be a 'DataFrame'")
    expect_error(spectraData(be) <- df[1:2, ], "have to match the length of")
})

test_that("$<-,MsBackendHdf5Peaks works", {
    be <- backendInitialize(MsBackendHdf5Peaks(), files = tempfile(), test_df)

    be$msLevel <- 1L
    expect_identical(msLevel(be), rep(1L, 3))

    ints <- lapply(test_df$intensity, function(z) z/2)
    be$intensity <- ints
    expect_identical(intensity(be), NumericList(ints, compress = FALSE))
    expect_identical(be@modCount, 1L)

    mzs <- lapply(test_df$mz, function(z) z/4)
    be$mz <- mzs
    expect_identical(mz(be), NumericList(mzs, compress = FALSE))
    expect_identical(be@modCount, 2L)

    expect_error(be$mz <- 4, "has to be a list")
    expect_error(be$intensity <- list(4), "has to match the length")
    expect_error(be$intensity <- list(4, 3, 5), "match the number of peaks")
})

test_that("[,MsBackendHdf5Peaks works", {
    fls <- c(tempfile(), tempfile())
    be <- backendInitialize(MsBackendHdf5Peaks(), files = fls,
                            spectraData = spectraData(sciex_mzr))
    idx <- sample(seq_along(be), 30)
    res <- be[idx]
    expect_identical(peaks(res), sciex_pks[idx])
    expect_identical(rtime(res), rtime(sciex_mzr)[idx])
    expect_identical(msLevel(res), msLevel(sciex_mzr)[idx])

    idx <- dataStorage(be) == fls[2]
    res <- be[idx, ]
    expect_true(all(dataStorage(res) == fls[2]))
    expect_identical(peaks(res), sciex_pks[idx])
})

test_that("backendMerge,MsBackendHdf5Peaks works", {
    tmp <- test_df
    tmp$rtime <- as.numeric(1:3)
    be1 <- backendInitialize(MsBackendHdf5Peaks(), file = tempfile(), tmp)
    tmp$rtime <- as.numeric(4:6)
    be2 <- backendInitialize(MsBackendHdf5Peaks(), file = tempfile(), tmp)
    tmp$rtime <- as.numeric(7:9)
    tmp$new_col <- "3rd backend"
    be3 <- backendInitialize(MsBackendHdf5Peaks(), file = tempfile(), tmp)

    be2$intensity <- be2$intensity / 2
    be3$intensity <- be3$intensity / 3
    be2$intensity <- be2$intensity + 3

    expect_true(be1@modCount == 0)
    expect_true(be2@modCount == 2)
    expect_true(be3@modCount == 1)

    ## Error if duplicated dataStorage
    expect_error(backendMerge(be1, be1), "same 'dataStorage' is not supported")

    ## Test merging of backends with different modCount.
    res <- backendMerge(be1, be2, be3)
    expect_identical(dataStorageLevels(res), c(dataStorageLevels(be1),
                                               dataStorageLevels(be2),
                                               dataStorageLevels(be3)))
    expect_identical(rtime(res), as.numeric(1:9))
    expect_identical(intensity(res), c(intensity(be1), intensity(be2),
                                       intensity(be3)))
    expect_identical(mz(res), c(mz(be1), mz(be2), mz(be3)))
    expect_identical(res@modCount, c(0L, 2L, 1L))

    ## Subsetting of backend with different modCount.
    tmp <- res[9:1]
    expect_identical(dataStorageLevels(tmp), rev(dataStorageLevels(res)))
    expect_identical(tmp@modCount, rev(res@modCount))
    expect_identical(peaks(tmp)[[1]], peaks(res)[[9]])
    expect_identical(peaks(tmp)[[8]], peaks(res)[[2]])

    tmp <- res[7]
    expect_identical(tmp@modCount, 1L)
    expect_identical(dataStorageLevels(tmp), dataStorageLevels(be3))

    ## Changing the data in the original object, invalidates the joined.
    be3$intensity <- be3$intensity * 2

    expect_true(validObject(be3))
    expect_error(mz(res))
})
