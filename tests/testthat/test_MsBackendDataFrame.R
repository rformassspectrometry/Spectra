test_df <- DataFrame(msLevel = c(1L, 2L, 2L), scanIndex = 4:6)
test_df$mz <- list(c(1.1, 1.3, 1.5), c(4.1, 5.1), c(1.6, 1.7, 1.8, 1.9))
test_df$intensity <- list(c(45.1, 34, 12), c(234.4, 1333), c(42.1, 34.2, 65, 6))

test_that("backendInitialize,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_true(validObject(be))
    be <- backendInitialize(be)
    expect_true(validObject(be))
    be <- backendInitialize(be, data = DataFrame(msLevel = 2L))
    expect_true(validObject(be))
    expect_equal(dataStorage(be), "<memory>")

    be_2 <- backendInitialize(be, data = data.frame(msLevel = 2L))
    expect_equal(be, be_2)

    expect_error(backendInitialize(be, data = 4), "has to be a")

    df <- test_df
    be <- backendInitialize(be, df)
    expect_true(validObject(be))
    expect_true(is(be@spectraData$mz, "NumericList"))
    expect_true(is(be@spectraData$intensity, "NumericList"))

    df$mz <- SimpleList(df$mz)
    df$intensity <- SimpleList(df$intensity)
    be <- backendInitialize(be, df)
    expect_true(validObject(be))
    expect_true(is(be@spectraData$mz, "NumericList"))
    expect_true(is(be@spectraData$intensity, "NumericList"))
    expect_identical(be@spectraData$dataStorage, rep("<memory>", 3))
    expect_identical(be$dataStorage, rep("<memory>", 3))

    df$mz <- NumericList(df$mz)
    df$intensity <- NumericList(df$intensity)
    be <- backendInitialize(be, df)
    expect_true(validObject(be))
    expect_true(is(be@spectraData$mz, "NumericList"))
    expect_true(is(be@spectraData$intensity, "NumericList"))
})

test_that("backendMerge,MsBackendDataFrame works", {
    df <- DataFrame(msLevel = c(1L, 2L, 2L), fromFile = 1L,
                    rtime = as.numeric(1:3))
    df2 <- DataFrame(msLevel = c(2L, 1L), fromFile = 1L,
                     rtime = c(4.1, 5.2), scanIndex = 1:2)
    df3 <- DataFrame(msLevel = c(1L, 2L), fromFile = 1L,
                     precScanNum = 1L, other_col = "z")
    be <- backendInitialize(MsBackendDataFrame(), df)
    be2 <- backendInitialize(MsBackendDataFrame(), df2)
    be3 <- backendInitialize(MsBackendDataFrame(), df3)

    expect_equal(backendMerge(be), be)
    expect_error(backendMerge(be, 4), "backends of the same type")

    res <- backendMerge(be, be2, be3)
    expect_true(is(res, "MsBackendDataFrame"))
    expect_identical(res@spectraData$dataStorage, rep("<memory>", 7))
    expect_identical(dataStorage(res), rep("<memory>", 7))
    expect_identical(msLevel(res), c(1L, 2L, 2L, 2L, 1L, 1L, 2L))
    expect_identical(rtime(res), c(1:3, 4.1, 5.2, NA, NA))
    expect_identical(res@spectraData$other_col,
                     c(rep(NA_character_, 5), "z", "z"))
    expect_true(is(be3@spectraData$precScanNum, "integer"))

    ## One backend with and one without m/z
    df2$mz <- list(c(1.1, 1.2), c(1.1, 1.2))
    df2$intensity <- list(c(12.4, 3), c(123.4, 1))
    be2 <- backendInitialize(MsBackendDataFrame(), df2)
    res <- backendMerge(be, be2, be3)
    expect_identical(lengths(mz(res)), c(0L, 0L, 0L, 2L, 2L, 0L, 0L))

    ## With different dataStorage
    be$dataStorage <- c("a", "a", "a")
    be3$dataStorage <- c("z", "b")

    res <- backendMerge(be, be2, be3)
    expect_identical(res$dataStorage,
                     c("a", "a", "a", "<memory>", "<memory>", "z", "b"))
    expect_identical(res@spectraData$dataStorage,
                     c("a", "a", "a", "<memory>", "<memory>", "z", "b"))
    expect_identical(rtime(res), c(1:3, 4.1, 5.2, NA, NA))
})

test_that("acquisitionNum, MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(acquisitionNum(be), integer())
    be <- backendInitialize(be, data = DataFrame(msLevel = c(1L, 2L)))
    expect_equal(acquisitionNum(be), c(NA_integer_, NA_integer_))
    be <- backendInitialize(be, data = DataFrame(msLevel = 1L,
                                                 acquisitionNum = 1:10))
    expect_equal(acquisitionNum(be), 1:10)
})

test_that("centroided, centroided<-, MsBackendDataFrame work", {
    be <- MsBackendDataFrame()
    expect_equal(centroided(be), logical())
    be <- backendInitialize(be, data = DataFrame(msLevel = c(1L, 2L)))
    expect_equal(centroided(be), c(NA, NA))
    expect_error(centroided(be) <- "a", "to be a 'logical'")
    expect_error(centroided(be) <- c(FALSE, TRUE, TRUE), "has to be a")
    centroided(be) <- c(TRUE, FALSE)
    expect_equal(centroided(be), c(TRUE, FALSE))
    centroided(be) <- FALSE
    expect_equal(centroided(be), c(FALSE, FALSE))
})

test_that("collisionEnergy, collisionEnergy<-,MsBackendDataFrame work", {
    be <- MsBackendDataFrame()
    expect_equal(collisionEnergy(be), numeric())
    be <- backendInitialize(be, data = DataFrame(msLevel = c(1L, 2L)))
    expect_equal(collisionEnergy(be), c(NA_real_, NA_real_))
    expect_error(collisionEnergy(be) <- "a", "to be a 'numeric'")
    expect_error(collisionEnergy(be) <- c(2.3), "has to be a")
    collisionEnergy(be) <- c(2.1, 3.2)
    expect_equal(collisionEnergy(be), c(2.1, 3.2))
})

test_that("dataOrigin,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(dataOrigin(be), character())
    be <- backendInitialize(be, data = DataFrame(msLevel = c(1L, 2L)))
    expect_identical(dataOrigin(be), rep(NA_character_, 2))
    expect_error(dataOrigin(be) <- "a", "of length 2")
    dataOrigin(be) <- c("b", "a")
    expect_identical(dataOrigin(be), c("b", "a"))
})

test_that("dataStorage,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(dataStorage(be), character())
    be <- backendInitialize(be, data = DataFrame(msLevel = c(1L, 2L)))
    expect_identical(dataStorage(be), rep("<memory>", 2))
    dataStorage(be) <- c("a", "b")
    expect_identical(dataStorage(be), c("a", "b"))
    expect_error(dataStorage(be) <- c("a", NA), "not allowed")
    expect_error(dataStorage(be) <- c("a", "b", "c"), "of length")
})

test_that("intensity,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(intensity(be), NumericList(compress = FALSE))
    be <- backendInitialize(be, data = DataFrame(msLevel = c(1L, 2L)))
    expect_equal(intensity(be), NumericList(numeric(), numeric(),
                                            compress = FALSE))
    df <- DataFrame(msLevel = c(1L, 2L))
    df$intensity <- list(1:4, c(2.1, 3.4))
    df$mz <- list(1:4, 1:2)
    be <- backendInitialize(be, data = df)
    expect_equal(intensity(be), NumericList(1:4, c(2.1, 3.4), compress = FALSE))
})

test_that("intensity<-,MsBackendDataFrame works", {
    be <- backendInitialize(MsBackendDataFrame(), data = test_df)

    new_ints <- lapply(test_df$intensity, function(z) z / 2)
    intensity(be) <- new_ints
    expect_identical(intensity(be), NumericList(new_ints, compress = FALSE))

    expect_error(intensity(be) <- 3, "has to be a list")
    expect_error(intensity(be) <- list(3, 2), "match the length")
    expect_error(intensity(be) <- list(3, 2, 4), "number of peaks")
})

test_that("ionCound,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(ionCount(be), numeric())
    df <- DataFrame(msLevel = c(1L, 2L))
    df$intensity <- list(1:4, c(2.1, 3.4))
    df$mz <- list(1:4, 1:2)
    be <- backendInitialize(be, data = df)
    expect_equal(ionCount(be), c(sum(1:4), sum(c(2.1, 3.4))))
})

test_that("isEmpty,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(isEmpty(be), logical())
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, data = df)
    expect_equal(isEmpty(be), c(TRUE, TRUE))
    df$intensity <- list(1:2, 1:5)
    df$mz <- list(1:2, 1:5)
    be <- backendInitialize(be, data = df)
    expect_equal(isEmpty(be), c(FALSE, FALSE))
})

test_that("isolationWindowLowerMz,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_identical(isolationWindowLowerMz(be), numeric())

    df <- DataFrame(msLevel = c(1L, 2L, 2L, 2L, 2L))
    be <- backendInitialize(MsBackendDataFrame(), df)

    expect_identical(isolationWindowLowerMz(be), rep(NA_real_, 5))
    isolationWindowLowerMz(be) <- c(NA_real_, 2, 2, 3, 3)
    expect_identical(isolationWindowLowerMz(be), c(NA_real_, 2, 2, 3, 3))

    df$isolationWindowLowerMz <- c(NA_real_, 4, 4, 3, 1)
    be <- backendInitialize(MsBackendDataFrame(), df)
    expect_identical(isolationWindowLowerMz(be), c(NA_real_, 4, 4, 3, 1))

    expect_error(isolationWindowLowerMz(be) <- 2.1, "of length 5")
})

test_that("isolationWindowTargetMz,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_identical(isolationWindowTargetMz(be), numeric())

    df <- DataFrame(msLevel = c(1L, 2L, 2L, 2L, 2L))
    be <- backendInitialize(MsBackendDataFrame(), df)

    expect_identical(isolationWindowTargetMz(be), rep(NA_real_, 5))
    isolationWindowTargetMz(be) <- c(NA_real_, 2, 2, 3, 3)
    expect_identical(isolationWindowTargetMz(be), c(NA_real_, 2, 2, 3, 3))

    df$isolationWindowTargetMz <- c(NA_real_, 4, 4, 3, 1)
    be <- backendInitialize(MsBackendDataFrame(), df)
    expect_identical(isolationWindowTargetMz(be), c(NA_real_, 4, 4, 3, 1))

    expect_error(isolationWindowTargetMz(be) <- 2.1, "of length 5")
})

test_that("isolationWindowUpperMz,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_identical(isolationWindowUpperMz(be), numeric())

    df <- DataFrame(msLevel = c(1L, 2L, 2L, 2L, 2L))
    be <- backendInitialize(MsBackendDataFrame(), df)

    expect_identical(isolationWindowUpperMz(be), rep(NA_real_, 5))
    isolationWindowUpperMz(be) <- c(NA_real_, 2, 2, 3, 3)
    expect_identical(isolationWindowUpperMz(be), c(NA_real_, 2, 2, 3, 3))

    df$isolationWindowUpperMz <- c(NA_real_, 4, 4, 3, 1)
    be <- backendInitialize(MsBackendDataFrame(), df)
    expect_identical(isolationWindowUpperMz(be), c(NA_real_, 4, 4, 3, 1))

    expect_error(isolationWindowUpperMz(be) <- 2.1, "of length 5")
})

test_that("length,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(length(be), 0)
    be <- new("MsBackendDataFrame", spectraData = DataFrame(a = 1:3,
                                                            dataStorage = "a"))
    expect_equal(length(be), 3)
})

test_that("msLevel,MsBackendDataFrame works", {
    be <- backendInitialize(MsBackendDataFrame(),
                            DataFrame(msLevel = c(1L, 2L, 1L)))
    expect_equal(msLevel(be), c(1, 2, 1))
    be <- backendInitialize(MsBackendDataFrame(),
                            DataFrame(scanIndex = 1:4))
    expect_equal(msLevel(be), rep(NA_integer_, 4))
})

test_that("mz,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(mz(be), NumericList(compress = FALSE))
    df <- DataFrame(msLevel = c(1L, 1L))
    be <- backendInitialize(be, data = df)
    expect_equal(mz(be), NumericList(numeric(), numeric(), compress = FALSE))
    df$intensity <- list(1:3, 4)
    df$mz <- list(1:3, c(2.1))
    be <- backendInitialize(be, data = df)
    expect_equal(mz(be), NumericList(1:3, 2.1, compress = FALSE))
})

test_that("mz<-,MsBackendDataFrame works", {
    be <- backendInitialize(MsBackendDataFrame(), data = test_df)

    new_mzs <- lapply(test_df$mz, function(z) z / 2)
    mz(be) <- new_mzs
    expect_identical(mz(be), NumericList(new_mzs, compress = FALSE))

    expect_error(mz(be) <- 3, "has to be a list")
    expect_error(mz(be) <- list(3, 2), "match the length")
    expect_error(mz(be) <- list(3, 2, 4), "number of peaks")
})

test_that("peaksData,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(peaksData(be), list())
    df <- DataFrame(msLevel = c(1L, 1L))
    be <- backendInitialize(be, data = df)
    expect_equal(peaksData(be), list(cbind(mz = numeric(),
                                           intensity = numeric()),
                                     cbind(mz = numeric(),
                                           intensity = numeric())))
    df$mz <- list(1:3, c(2.1))
    df$intensity <- list(1:3, 4)
    be <- backendInitialize(be, data = df)
    expect_equal(peaksData(be), list(cbind(mz = 1:3, intensity = 1:3),
                                     cbind(mz = 2.1, intensity = 4)))
})

test_that("peaksData<-,MsBackendDataFrame works", {
    be <- backendInitialize(MsBackendDataFrame(), data = test_df)

    pks <- lapply(peaksData(be), function(z) z / 2)
    peaksData(be) <- pks
    expect_identical(peaksData(be), pks)

    expect_error(peaksData(be) <- 3, "has to be a list")
    expect_error(peaksData(be) <- list(3, 2), "match length")
    expect_error(peaksData(be) <- list(3, 2, 4), "dimensions")
})

test_that("lengths,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(lengths(be), integer())
    df <- DataFrame(msLevel = c(1L, 1L))
    be <- backendInitialize(be, data = df)
    expect_equal(lengths(be), c(0L, 0L))
    df$mz <- list(1:3, c(2.1))
    df$intensity <- list(1:3, 4)
    be <- backendInitialize(be, data = df)
    expect_equal(lengths(be), c(3L, 1L))
})

test_that("polarity, polarity<- MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(polarity(be), integer())
    df <- DataFrame(msLevel = c(1L, 1L))
    be <- backendInitialize(be, data = df)
    expect_equal(polarity(be), c(NA_integer_, NA_integer_))
    expect_error(polarity(be) <- "a", "has to be an 'integer'")
    expect_error(polarity(be) <- c(1L, 1L, 2L), "has to be")
    polarity(be) <- 0
    expect_equal(polarity(be), c(0L, 0L))
    polarity(be) <- 1:2
    expect_equal(polarity(be), 1:2)
})

test_that("precScanNum,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(precScanNum(be), integer())
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, data = df)
    expect_equal(precScanNum(be), c(NA_integer_, NA_integer_))
    df$precScanNum <- c(0L, 1L)
    be <- backendInitialize(be, data = df)
    expect_equal(precScanNum(be), c(0L, 1L))
})

test_that("precursorCharge,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(precursorCharge(be), integer())
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, data = df)
    expect_equal(precursorCharge(be), c(NA_integer_, NA_integer_))
    df$precursorCharge <- c(-1L, 1L)
    be <- backendInitialize(be, data = df)
    expect_equal(precursorCharge(be), c(-1L, 1L))
})

test_that("precursorIntensity,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(precursorIntensity(be), numeric())
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, data = df)
    expect_equal(precursorIntensity(be), c(NA_real_, NA_real_))
    df$precursorIntensity <- c(134.4, 4322.2)
    be <- backendInitialize(be, data = df)
    expect_equal(precursorIntensity(be), c(134.4, 4322.2))
})

test_that("precursorMz,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(precursorMz(be), numeric())
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, data = df)
    expect_equal(precursorMz(be), c(NA_real_, NA_real_))
    df$precursorMz <- c(134.4, 342.2)
    be <- backendInitialize(be, data = df)
    expect_equal(precursorMz(be), c(134.4, 342.2))
})

test_that("rtime, rtime<-,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(rtime(be), numeric())
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, data = df)
    expect_equal(rtime(be), c(NA_real_, NA_real_))
    expect_error(rtime(be) <- "2", "has to be a 'numeric'")
    expect_error(rtime(be) <- 1:4, "of length 2")
    rtime(be) <- c(123, 124)
    expect_equal(rtime(be), c(123, 124))
})

test_that("scanIndex,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(scanIndex(be), integer())
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, data = df)
    expect_equal(scanIndex(be), c(NA_integer_, NA_integer_))
    df$scanIndex <- c(1L, 2L)
    be <- backendInitialize(be, data = df)
    expect_equal(scanIndex(be), c(1L, 2L))
})

test_that("smoothed, smoothed<-,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(smoothed(be), logical())
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, data = df)
    expect_equal(smoothed(be), c(NA, NA))
    expect_error(smoothed(be) <- "2", "has to be a 'logical'")
    expect_error(smoothed(be) <- c(TRUE, TRUE, FALSE), "of length 1 or 2")
    smoothed(be) <- c(TRUE, FALSE)
    expect_equal(smoothed(be), c(TRUE, FALSE))
    smoothed(be) <- c(TRUE)
    expect_equal(smoothed(be), c(TRUE, TRUE))
})

test_that("spectraNames, spectraNames<-,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_null(spectraNames(be))
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, data = df)
    expect_null(spectraNames(be))
    df <- DataFrame(msLevel = c(1L, 2L))
    rownames(df) <- c("sp_1", "sp_2")
    be <- backendInitialize(be, data = df)
    expect_equal(spectraNames(be), c("sp_1", "sp_2"))
    expect_error(spectraNames(be) <- "a", "rownames length")
    spectraNames(be) <- c("a", "b")
    expect_equal(spectraNames(be), c("a", "b"))
})

test_that("tic,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(tic(be), numeric())
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, data = df)
    expect_equal(tic(be), c(NA_real_, NA_real_))
    expect_equal(tic(be, initial = FALSE), c(0, 0))
    df$totIonCurrent <- c(5, 3)
    be <- backendInitialize(be, data = df)
    expect_equal(tic(be), c(5, 3))
    expect_equal(tic(be, initial = FALSE), c(0, 0))
    df$intensity <- list(5:7, 1:4)
    df$mz <- list(1:3, 1:4)
    be <- backendInitialize(be, data = df)
    expect_equal(tic(be, initial = FALSE), c(sum(5:7), sum(1:4)))
})

test_that("spectraVariables,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(spectraVariables(be), names(.SPECTRA_DATA_COLUMNS))
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, data = df)
    expect_equal(spectraVariables(be), names(.SPECTRA_DATA_COLUMNS))
    df$other_column <- 3
    be <- backendInitialize(be, data = df)
    expect_equal(spectraVariables(be), c(names(.SPECTRA_DATA_COLUMNS),
                                         "other_column"))
})

test_that("spectraData, spectraData<-, MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)
    expect_equal(colnames(res), names(.SPECTRA_DATA_COLUMNS))

    df <- DataFrame(scanIndex = 1:2, a = "a", b = "b")
    be <- backendInitialize(be, data = df)

    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(all(names(.SPECTRA_DATA_COLUMNS) %in% colnames(res)))
    expect_equal(res$a, c("a", "a"))
    expect_equal(res$b, c("b", "b"))

    expect_error(spectraData(be) <- data.frame(4), "not valid")

    res <- spectraData(be, "msLevel")
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), "msLevel")
    expect_equal(res$msLevel, c(NA_integer_, NA_integer_))

    res <- spectraData(be, c("mz"))
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), "mz")
    expect_equal(res$mz, NumericList(numeric(), numeric(), compress = FALSE))

    res <- spectraData(be, c("a", "intensity"))
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), c("a", "intensity"))
    expect_equal(res$intensity, NumericList(numeric(), numeric(),
                                            compress = FALSE))
    expect_equal(res$a, c("a", "a"))

    spectraData(be) <- DataFrame(mzLevel = c(3L, 4L),
                                 rtime = c(1.2, 1.4), other_col = "b")
    expect_identical(rtime(be), c(1.2, 1.4))
    expect_true(any(spectraVariables(be) == "other_col"))
    expect_identical(spectraData(be, "other_col")[, 1], c("b", "b"))

    expect_error(spectraData(be) <- DataFrame(msLevel = 1:3),
                 "with 2 rows")
})

test_that("show,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    show(be)
    df <- DataFrame(rt = c(1.2, 1.3))
    be <- backendInitialize(be, df)
})

test_that("[,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_error(be[1])
    df <- DataFrame(scanIndex = 1:2, a = "a", b = "b")
    be <- backendInitialize(be, df)
    res <- be[1]
    expect_true(validObject(res))
    expect_equal(be@spectraData[1, ], res@spectraData[1, ])
    res <- be[2]
    expect_true(validObject(res))
    expect_equal(be@spectraData[2, ], res@spectraData[1, ])
    res <- be[2:1]
    expect_true(validObject(res))
    expect_equal(be@spectraData[2:1, ], res@spectraData)

    res <- be[c(FALSE, FALSE)]
    expect_true(validObject(res))
    expect_true(length(res) == 0)
    res <- be[c(FALSE, TRUE)]
    expect_true(validObject(res))
    expect_equal(be@spectraData[2, ], res@spectraData[1, ])

    expect_error(be[TRUE], "match the length of")
    expect_error(be["a"], "does not have names")

    df <- DataFrame(scanIndex = c(1L, 2L, 1L, 2L),
                    file = c("a", "a", "b", "b"))
    be <- backendInitialize(be, df)
    dataStorage(be) <- c("1", "1", "2", "2")
    res <- be[3]
    expect_true(validObject(res))
    expect_equal(dataStorage(res), "2")
    expect_equal(res@spectraData$file, "b")

    res <- be[c(3, 1)]
    expect_true(validObject(res))
    expect_equal(dataStorage(res), c("2", "1"))
    expect_equal(res@spectraData$file, c("b", "a"))
})

test_that("selectSpectraVariables,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    res <- selectSpectraVariables(be, c("dataStorage", "msLevel"))

    df <- DataFrame(msLevel = 1:2, rtime = c(2.3, 1.2),
                    other_col = 2)
    be <- backendInitialize(MsBackendDataFrame(), df)

    res <- selectSpectraVariables(be, c("dataStorage", "other_col"))
    expect_equal(colnames(res@spectraData), c("dataStorage", "other_col"))
    expect_equal(msLevel(res), c(NA_integer_, NA_integer_))

    res <- selectSpectraVariables(be, c("dataStorage", "rtime"))
    expect_equal(colnames(res@spectraData), c("dataStorage", "rtime"))

    expect_error(selectSpectraVariables(be, "rtime"), "dataStorage is/are missing")
    expect_error(selectSpectraVariables(be, "something"),
                 "something not available")
})

test_that("$,$<-,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_identical(be$msLevel, integer())

    df <- DataFrame(msLevel = 1:2, rtime = c(2.3, 1.2),
                    other_col = 2)
    be <- backendInitialize(MsBackendDataFrame(), df)
    expect_identical(be$msLevel, 1:2)
    expect_identical(be$other_col, c(2, 2))

    expect_error(be$not_there, "not available")

    be$other_col <- 4
    expect_equal(be$other_col, c(4, 4))

    df$mz <- list(1:3, 1:4)
    df$intensity <- list(c(3, 3, 3), c(4, 4, 4, 4))
    be <- backendInitialize(MsBackendDataFrame(), df)
    be$intensity <- list(c(5, 5, 5), 1:4)
    expect_equal(be$intensity, NumericList(c(5, 5, 5), 1:4, compress = FALSE))
})

test_that("filterAcquisitionNum,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(be, filterAcquisitionNum(be, n = 4))

    df <- DataFrame(acquisitionNum = c(1L, 2L, 3L, 2L, 3L, 1L, 2L, 4L),
                    msLevel = 1L)
    be <- backendInitialize(MsBackendDataFrame(), df)
    dataStorage(be) <- c("1", "1", "1", "2", "2", "3", "3", "3")
    res <- filterAcquisitionNum(be, n = c(2L, 4L))
    expect_equal(length(res), 4)
    expect_equal(dataStorage(res), c("1", "2", "3", "3"))
    expect_equal(acquisitionNum(res), c(2L, 2L, 2L, 4L))

    res <- filterAcquisitionNum(be, n = 2L, dataStorage = "2")
    expect_equal(dataStorage(res), c("1", "1", "1", "2", "3", "3", "3"))
    expect_equal(acquisitionNum(res), c(1L, 2L, 3L, 2L, 1L, 2L, 4L))

    expect_error(filterAcquisitionNum(be, n = "a"), "integer representing")
    expect_equal(filterAcquisitionNum(be), be)
})

test_that("filterDataOrigin,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(be, filterDataOrigin(be))
    df <- DataFrame(acquisitionNum = c(1L, 2L, 3L, 2L, 3L, 1L, 2L, 4L),
                    rtime = as.numeric(1:8),
                    msLevel = 1L)
    be <- backendInitialize(MsBackendDataFrame(), df)
    dataStorage(be) <- c("1", "1", "1", "2", "2", "3", "3", "3")
    expect_true(length(filterDataOrigin(be, c(3, 1))) == 0)
    expect_equal(be, filterDataOrigin(be, NA_character_))

    dataOrigin(be) <- c("1", "1", "1", "2", "2", "3", "3", "3")
    res <- filterDataOrigin(be, c(3, 1))
    expect_equal(length(res), 6)
    expect_equal(dataOrigin(res), c("3", "3", "3", "1", "1", "1"))
    expect_equal(rtime(res), c(6, 7, 8, 1, 2, 3))
    expect_equal(unique(res$dataStorage), c("3", "1"))

    res <- filterDataOrigin(be, c("2", "1"))
    expect_equal(length(res), 5)
    expect_equal(dataOrigin(res), c("2", "2", "1", "1", "1"))
    expect_equal(rtime(res), c(4, 5, 1, 2, 3))

    res <- filterDataOrigin(be, 2)
    expect_equal(rtime(res), c(4, 5))
    expect_equal(unique(res$dataStorage), "2")

    res <- filterDataOrigin(be, c(2, 3))
    expect_equal(rtime(res), c(4, 5, 6, 7, 8))
    expect_equal(unique(res$dataStorage), c("2", "3"))

    res <- filterDataOrigin(be)
    expect_equal(res, be)

    df$dataOrigin <- "1"
    be <- backendInitialize(MsBackendDataFrame(), df)
    res <- filterDataOrigin(be, 1L)
    expect_equal(res, be)
})

test_that("filterDataStorage,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(be, filterDataStorage(be))
    df <- DataFrame(acquisitionNum = c(1L, 2L, 3L, 2L, 3L, 1L, 2L, 4L),
                    rtime = as.numeric(1:8),
                    msLevel = 1L)
    be <- backendInitialize(MsBackendDataFrame(), df)
    dataStorage(be) <- c("1", "1", "1", "2", "2", "3", "3", "3")
    res <- filterDataStorage(be, c(3, 1))
    expect_equal(length(res), 6)
    expect_equal(dataStorage(res), c("3", "3", "3", "1", "1", "1"))
    expect_equal(rtime(res), c(6, 7, 8, 1, 2, 3))
    expect_equal(unique(res$dataStorage), c("3", "1"))

    res <- filterDataStorage(be, c("2", "1"))
    expect_equal(length(res), 5)
    expect_equal(dataStorage(res), c("2", "2", "1", "1", "1"))
    expect_equal(rtime(res), c(4, 5, 1, 2, 3))

    res <- filterDataStorage(be, 2)
    expect_equal(rtime(res), c(4, 5))
    expect_equal(unique(res$dataStorage), "2")

    res <- filterDataStorage(be, c(2, 3))
    expect_equal(rtime(res), c(4, 5, 6, 7, 8))
    expect_equal(unique(res$dataStorage), c("2", "3"))

    res <- filterDataStorage(be)
    expect_equal(res, be)
})

test_that("filterEmptySpectra,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(filterEmptySpectra(be), be)

    df <- DataFrame(msLevel = c(1L, 1L, 1L), rtime = c(1.2, 1.3, 1.5))
    df$mz <- SimpleList(1:3, 1:4, integer())
    df$intensity <- SimpleList(c(5, 12, 9), c(6, 9, 12, 5), numeric())
    be <- backendInitialize(MsBackendDataFrame(), df)

    res <- filterEmptySpectra(be)
    expect_true(length(res) == 2)
    expect_equal(rtime(res), c(1.2, 1.3))

    res <- filterEmptySpectra(be[3])
    expect_true(length(res) == 0)
})

test_that("filterIsolationWindow,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(be, filterIsolationWindow(be))

    df <- DataFrame(msLevel = c(1L, 2L, 2L, 2L, 2L, 2L),
                    rtime = as.numeric(1:6))
    be <- backendInitialize(MsBackendDataFrame(), df)

    expect_error(filterIsolationWindow(be, c(1.2, 1.3)), "single m/z value")
    res <- filterIsolationWindow(be, 1.2)
    expect_true(length(res) == 0)

    df$isolationWindowLowerMz <- c(NA_real_, 1.2, 2.1, 3.1, 3, 3.4)
    df$isolationWindowUpperMz <- c(NA_real_, 2.2, 3.2, 4.1, 4, 4.4)
    be <- backendInitialize(MsBackendDataFrame(), df)

    res <- filterIsolationWindow(be, 3)
    expect_equal(rtime(res), c(3, 5))
})

test_that("filterMsLevel,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(be, filterMsLevel(be))

    df <- DataFrame(msLevel = c(1L, 2L, 1L, 2L, 3L, 2L),
                    rtime = as.numeric(1:6))
    be <- backendInitialize(MsBackendDataFrame(), df)
    res <- filterMsLevel(be, 2L)
    expect_true(all(msLevel(res) == 2))
    expect_equal(rtime(res), c(2, 4, 6))

    res <- filterMsLevel(be, c(3L, 2L))
    expect_equal(rtime(res), c(2, 4, 5, 6))

    res <- filterMsLevel(be, c(3L, 5L))
    expect_equal(rtime(res), 5)
})

test_that("filterPolarity,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(be, filterPolarity(be))

    df <- DataFrame(msLevel = c(1L, 2L, 3L, 1L, 2L, 3L),
                    rtime = as.numeric(1:6), polarity = c(1L, 1L, -1L, 0L, 1L, 0L))
    be <- backendInitialize(MsBackendDataFrame(), df)
    res <- filterPolarity(be, c(1, 2))
    expect_true(all(polarity(res) == 1L))
    expect_equal(rtime(res), c(1, 2, 5))

    res <- filterPolarity(be, c(0L, -1L))
    expect_true(all(polarity(res) %in% c(0L, -1L)))
    expect_equal(rtime(res), c(3, 4, 6))
})

test_that("filterPrecursorMz,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(be, filterPrecursorMz(be))

    df <- DataFrame(msLevel = c(1L, 2L, 2L, 2L),
                    rtime = as.numeric(1:4))
    be <- backendInitialize(MsBackendDataFrame(), df)
    res <- filterPrecursorMz(be, 1.4)
    expect_true(length(res) == 0)

    df$precursorMz <- c(NA_real_, 4.43, 4.4312, 5.4)
    be <- backendInitialize(MsBackendDataFrame(), df)

    res <- filterPrecursorMz(be, 5.4)
    expect_equal(rtime(res), 4)

    res <- filterPrecursorMz(be, 4.4311)
    expect_true(length(res) == 0)

    res <- filterPrecursorMz(be, mz = 4.4311 + ppm(c(-4.4311, 4.4311), 40))
    expect_equal(rtime(res), 3)
})

test_that("filterPrecursorScan,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(be, filterPrecursorScan(be))
    expect_true(length(filterPrecursorScan(be, 3)) == 0)

    df <- DataFrame(msLevel = c(1L, 2L, 3L, 4L, 1L, 2L, 3L),
                    rtime = as.numeric(1:7))
    be <- backendInitialize(MsBackendDataFrame(), df)
    res <- filterPrecursorScan(be, 2L)
    expect_true(length(res) == 0)

    df$acquisitionNum <- c(1L, 2L, 3L, 4L, 1L, 2L, 6L)
    be <- backendInitialize(MsBackendDataFrame(), df)
    res <- filterPrecursorScan(be, 2L)
    expect_equal(acquisitionNum(res), c(2L, 2L))
    expect_equal(rtime(res), c(2, 6))

    df$precScanNum <- c(0L, 1L, 2L, 3L, 0L, 1L, 5L)
    be <- backendInitialize(MsBackendDataFrame(), df)
    res <- filterPrecursorScan(be, 2L)
    expect_equal(rtime(res), 1:6)

    res <- filterPrecursorScan(be, 5L)
    expect_true(length(rtime(res)) == 0)

    res <- filterAcquisitionNum(be, 6L)
    expect_equal(rtime(res), 7)
})

test_that("filterRt,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(filterRt(be), be)
    expect_true(length(filterRt(be, rt = 1:2)) == 0)

    df <- DataFrame(rtime = c(1, 2, 3, 4, 1, 2, 6, 7, 9),
                    index = 1:9)
    be <- backendInitialize(MsBackendDataFrame(), df)

    res <- filterRt(be, c(4, 6))
    expect_equal(rtime(res), c(4, 6))
    expect_equal(res@spectraData$index, c(4, 7))

    res <- filterRt(be, c(4, 6), 2L)
    expect_equal(res, be)

    df$msLevel <- c(1L, 2L, 1L, 2L, 1L, 2L, 2L, 2L, 2L)
    be <- backendInitialize(MsBackendDataFrame(), df)

    res <- filterRt(be, c(2, 6), msLevel = 1L)
    expect_equal(rtime(res), c(2, 3, 4, 2, 6, 7, 9))


    df$rtime <- NULL
    be <- backendInitialize(MsBackendDataFrame(), df)

    res <- filterRt(be, c(2, 6), msLevel = 1L)
    expect_equal(res, filterMsLevel(be, 2L))

    res <- filterRt(be, c(2, 6))
    expect_true(length(res) == 0)
})

test_that("split,MsBackendDataFrame works", {
    msb <- sciex_mzr
    msbl <- split(msb, f = msb$dataStorage)
    expect_true(is(msbl[[1]], "MsBackendDataFrame"))
    expect_identical(msLevel(msb)[msb$dataStorage == msb$dataStorage[1]],
                     msLevel(msbl[[1]]))
    expect_identical(intensity(msb)[msb$dataStorage == msb$dataStorage[1]],
                     intensity(msbl[[1]]))

    msb2 <- backendMerge(msbl)
    expect_identical(msb2, msb)
})

test_that("isCentroided,MsBackendDataFrame works", {
    msb <- MsBackendDataFrame()
    msb <- backendInitialize(msb, test_df)
    expect_true(all(is.na(isCentroided(msb))))
})

test_that("dropNaSpectraVariables works with MsBackendDataFrame", {
    b <- MsBackendDataFrame()
    expect_true(length(dropNaSpectraVariables(b)) == 0)
    b <- backendInitialize(b, test_df)
    res <- dropNaSpectraVariables(b)
    expect_equal(spectraVariables(b), spectraVariables(res))
    expect_equal(mz(b), mz(res))
    expect_equal(intensity(b), intensity(res))

    b$other_col <- NA
    res <- dropNaSpectraVariables(b)
    expect_true(!any(spectraVariables(res) == "other_col"))
})
