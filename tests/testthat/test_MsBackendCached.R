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
    expect_equal(sort(colnames(res)), sort(names(.SPECTRA_DATA_COLUMNS)))

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

test_that("spectraData,MsBackendCached works", {
    be <- new("MsBackendCached")
    be@nspectra <- 14L
    res <- spectraData(be)
    expect_s4_class(res, "DataFrame")
    expect_equal(.valid_column_datatype(res), NULL)
    expect_true(all(names(.SPECTRA_DATA_COLUMNS) %in% colnames(res)))

    expect_error(spectraData(be, "other_col"), "not available")

    res <- spectraData(be, "smoothed")
    expect_s4_class(res, "DataFrame")
    expect_equal(colnames(res), "smoothed")
    expect_true(all(is.na(res$smoothed)))
    expect_true(is.logical(res$smoothed))
})

test_that("[,MsBackendCached works", {
    be <- MsBackendCached()

    be <- backendInitialize(be, nspectra = 10)
    res <- be[c(1, 4, 3), ]
    expect_true(length(res) == 3)
    expect_true(nrow(res@localData) == 3)
    res_2 <- extractByIndex(be, c(1, 4, 3))
    expect_equal(res, res_2)

    df <- data.frame(msLevel = 1L, b = 1:6)
    be <- backendInitialize(be, data = df)
    res <- be[c(6, 1, 3)]
    expect_true(length(res) == 3)
    expect_equal(res@localData$b, c(6, 1, 3))
    res_2 <- extractByIndex(be, c(6, 1, 3))
    expect_equal(res, res_2)

    res <- be[c(6, 1, 3, 1)]
    expect_true(length(res) == 4)
    expect_equal(res@localData$b, c(6, 1, 3, 1))
    res_2 <- extractByIndex(be, c(6, 1, 3, 1))
    expect_equal(res, res_2)

    expect_equal(extractByIndex(be), be)
})

test_that("$,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)

    expect_equal(be$msLevel, rep(NA_integer_, 5))

    be <- backendInitialize(be, data = data.frame(msLevel = 1:3))
    expect_equal(be$msLevel, 1:3)

    expect_error(be$not, "not available")
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

test_that("acquisitionNum,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(acquisitionNum(be), rep(NA_integer_, length(be)))

    be <- backendInitialize(
        MsBackendCached(), data = data.frame(acquisitionNum = 1:3))
    expect_equal(acquisitionNum(be), 1:3)
})

test_that("centroided,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(centroided(be), rep(NA, length(be)))

    centroided(be) <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
    expect_equal(centroided(be), c(TRUE, FALSE, TRUE, FALSE, FALSE))
})

test_that("collisionEnergy,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(collisionEnergy(be), rep(NA_real_, length(be)))

    collisionEnergy(be) <- c(1.2, 2, 4, 2, 5)
    expect_equal(collisionEnergy(be), c(1.2, 2, 4, 2, 5))
})

test_that("dataOrigin,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(dataOrigin(be), rep(NA_character_, length(be)))

    dataOrigin(be) <- "unknown"
    expect_equal(dataOrigin(be), rep("unknown", length(be)))
})

test_that("msLevel,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(msLevel(be), rep(NA_integer_, length(be)))

    be <- backendInitialize(be, data = data.frame(msLevel = 1:4))
    expect_equal(msLevel(be), 1:4)
})

test_that("isolationWindowLowerMz,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(isolationWindowLowerMz(be), rep(NA_real_, length(be)))

    vals <- c(1.3, 4.2, 4.2, 4.5, 6.3)
    isolationWindowLowerMz(be) <- vals
    expect_equal(isolationWindowLowerMz(be), vals)
})

test_that("isolationWindowTargetMz,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(isolationWindowTargetMz(be), rep(NA_real_, length(be)))

    vals <- c(1.3, 4.2, 4.2, 4.5, 6.3)
    isolationWindowTargetMz(be) <- vals
    expect_equal(isolationWindowTargetMz(be), vals)
})

test_that("isolationWindowUpperMz,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(isolationWindowUpperMz(be), rep(NA_real_, length(be)))

    vals <- c(1.3, 4.2, 4.2, 4.5, 6.3)
    isolationWindowUpperMz(be) <- vals
    expect_equal(isolationWindowUpperMz(be), vals)
})

test_that("polarity,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(polarity(be), rep(NA_integer_, length(be)))

    vals <- c(1L, 1L, 0L, 1L, 0L)
    polarity(be) <- vals
    expect_equal(polarity(be), vals)
})

test_that("precursorCharge,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(precursorCharge(be), rep(NA_integer_, length(be)))

    be <- backendInitialize(be, data = data.frame(precursorCharge = 1:3))
    expect_equal(precursorCharge(be), 1:3)
})

test_that("precursorIntensity,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(precursorIntensity(be), rep(NA_real_, length(be)))

    be <- backendInitialize(be, data = data.frame(
                                    precursorIntensity = c(1.2, 1.5)))
    expect_equal(precursorIntensity(be), c(1.2, 1.5))
})

test_that("precursorMz,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(precursorMz(be), rep(NA_real_, length(be)))

    be <- backendInitialize(be, data = data.frame(
                                    precursorMz = c(1.2, 1.5)))
    expect_equal(precursorMz(be), c(1.2, 1.5))
})

test_that("rtime,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(rtime(be), rep(NA_real_, length(be)))

    vals <- c(1.4, 1.6, 1.8, 3.1, 5.2)
    rtime(be) <- vals
    expect_equal(rtime(be), vals)
})

test_that("scanIndex,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(scanIndex(be), rep(NA_integer_, length(be)))

    be <- backendInitialize(be, data = data.frame(scanIndex = 4:6))
    expect_equal(scanIndex(be), 4:6)
})

test_that("smoothed,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 5)
    expect_equal(smoothed(be), rep(NA, length(be)))

    smoothed(be) <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
    expect_equal(smoothed(be), c(FALSE, TRUE, FALSE, TRUE, FALSE))
})

test_that("intensity,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 3)
    res <- intensity(be)
    expect_true(is(res, "NumericList"))
    expect_true(all(lengths(res) == 0))
})

test_that("mz,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 3)
    res <- mz(be)
    expect_true(is(res, "NumericList"))
    expect_true(all(lengths(res) == 0))
})

test_that("ionCount,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 4)
    res <- ionCount(be)
    expect_true(all(res == 0))
})

test_that("isEmpty,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 4)
    res <- isEmpty(be)
    expect_true(all(res))
})

test_that("lengths,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 4)
    res <- lengths(be)
    expect_true(all(res == 0))
})

test_that("precursorMz<-,MsBackendCached works", {
    be <- backendInitialize(MsBackendCached(), nspectra = 4)
    expect_true(all(is.na(precursorMz(be))))
    precursorMz(be) <- c(1.1, 1.2, 1.3, 1.34)
    expect_equal(precursorMz(be), c(1.1, 1.2, 1.3, 1.34))
})
