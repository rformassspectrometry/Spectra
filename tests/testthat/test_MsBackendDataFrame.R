test_that("backendInitialize,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    be <- backendInitialize(be)
    expect_true(validObject(be))
    expect_true(length(be@files) == 0)
    be <- backendInitialize(be, files = NA_character_,
                            spectraData = DataFrame(fromFile = 1L))
    expect_true(validObject(be))
    expect_error(backendInitialize(be, spectraData = DataFrame(msLevel = 1L)),
                 "fromFile is/are")
    expect_error(backendInitialize(MsBackendDataFrame(),
                                   spectraData = DataFrame(msLevel = 1L,
                                                           fromFile = 1L)),
                 "'files' can not be empty")
})

test_that("acquisitionNum, MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(acquisitionNum(be), integer())
    be <- backendInitialize(be, files = NA_character_,
                            spectraData = DataFrame(fromFile = 1L,
                                                    msLevel = c(1L, 2L)))
    expect_equal(acquisitionNum(be), c(NA_integer_, NA_integer_))
    be <- backendInitialize(be, files = NA_character_,
                            spectraData = DataFrame(fromFile = 1L,
                                                    msLevel = 1L,
                                                    acquisitionNum = 1:10))
    expect_equal(acquisitionNum(be), 1:10)
})

test_that("centroided, centroided<-, MsBackendDataFrame work", {
    be <- MsBackendDataFrame()
    expect_equal(centroided(be), logical())
    be <- backendInitialize(be, files = NA_character_,
                            spectraData = DataFrame(fromFile = 1L,
                                                    msLevel = c(1L, 2L)))
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
    be <- backendInitialize(be, files = NA_character_,
                            spectraData = DataFrame(fromFile = 1L,
                                                    msLevel = c(1L, 2L)))
    expect_equal(collisionEnergy(be), c(NA_real_, NA_real_))
    expect_error(collisionEnergy(be) <- "a", "to be a 'numeric'")
    expect_error(collisionEnergy(be) <- c(2.3), "has to be a")
    collisionEnergy(be) <- c(2.1, 3.2)
    expect_equal(collisionEnergy(be), c(2.1, 3.2))
})

test_that("fromFile,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(fromFile(be), integer())
    be <- backendInitialize(be, files = NA_character_,
                            spectraData = DataFrame(fromFile = 1L,
                                                    msLevel = c(1L, 2L)))
    expect_equal(fromFile(be), c(1L, 1L))
})

test_that("intensity,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(intensity(be), list())
    be <- backendInitialize(be, files = NA_character_,
                            spectraData = DataFrame(fromFile = 1L,
                                                    msLevel = c(1L, 2L)))
    expect_equal(intensity(be), list(numeric(), numeric()))
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    df$intensity <- list(1:4, c(2.1, 3.4))
    df$mz <- list(1:4, 1:2)
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(intensity(be), list(1:4, c(2.1, 3.4)))
})

test_that("ionCound,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(ionCount(be), numeric())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    df$intensity <- list(1:4, c(2.1, 3.4))
    df$mz <- list(1:4, 1:2)
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(ionCount(be), c(sum(1:4), sum(c(2.1, 3.4))))
})

test_that("isEmpty,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(isEmpty(be), logical())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(isEmpty(be), c(TRUE, TRUE))
    df$intensity <- list(1:2, 1:5)
    df$mz <- list(1:2, 1:5)
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(isEmpty(be), c(FALSE, FALSE))
})

test_that("length,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(length(be), 0)
    be <- new("MsBackendDataFrame", spectraData = DataFrame(a = 1:3, fromFile = 1L),
              files = NA_character_, modCount = 0L)
    expect_equal(length(be), 3)
})

test_that("msLevel,MsBackendDataFrame works", {
    be <- backendInitialize(MsBackendDataFrame(), NA_character_,
                            DataFrame(msLevel = c(1L, 2L, 1L),
                                      fromFile = 1L))
    expect_equal(msLevel(be), c(1, 2, 1))
    be <- backendInitialize(MsBackendDataFrame(), NA_character_,
                            DataFrame(scanIndex = 1:4,
                                      fromFile = 1L))
    expect_equal(msLevel(be), rep(NA_integer_, 4))
})

test_that("mz,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(mz(be), list())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 1L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(mz(be), list(numeric(), numeric()))
    df$intensity <- list(1:3, 4)
    df$mz <- list(1:3, c(2.1))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(mz(be), list(1:3, 2.1))
})

test_that("peaks,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(peaks(be), list())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 1L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(peaks(be), list(cbind(mz = numeric(), intensity = numeric()),
                                 cbind(mz = numeric(), intensity = numeric())))
    df$mz <- list(1:3, c(2.1))
    df$intensity <- list(1:3, 4)
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(peaks(be), list(cbind(mz = 1:3, intensity = 1:3),
                                 cbind(mz = 2.1, intensity = 4)))
})

test_that("peaksCount,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(peaksCount(be), integer())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 1L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(peaksCount(be), c(0L, 0L))
    df$mz <- list(1:3, c(2.1))
    df$intensity <- list(1:3, 4)
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(peaksCount(be), c(3L, 1L))
})

test_that("polarity, polarity<- MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(polarity(be), integer())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 1L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
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
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(precScanNum(be), c(NA_integer_, NA_integer_))
    df$precScanNum <- c(0L, 1L)
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(precScanNum(be), c(0L, 1L))
})

test_that("precursorCharge,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(precursorCharge(be), integer())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(precursorCharge(be), c(NA_integer_, NA_integer_))
    df$precursorCharge <- c(-1L, 1L)
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(precursorCharge(be), c(-1L, 1L))
})

test_that("precursorIntensity,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(precursorIntensity(be), numeric())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(precursorIntensity(be), c(NA_real_, NA_real_))
    df$precursorIntensity <- c(134.4, 4322.2)
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(precursorIntensity(be), c(134.4, 4322.2))
})

test_that("precursorMz,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(precursorMz(be), numeric())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(precursorMz(be), c(NA_real_, NA_real_))
    df$precursorMz <- c(134.4, 342.2)
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(precursorMz(be), c(134.4, 342.2))
})

test_that("rtime, rtime<-,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(rtime(be), numeric())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(rtime(be), c(NA_real_, NA_real_))
    expect_error(rtime(be) <- "2", "has to be a 'numeric'")
    expect_error(rtime(be) <- 1:4, "of length 2")
    rtime(be) <- c(123, 124)
    expect_equal(rtime(be), c(123, 124))
})

test_that("scanIndex,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(scanIndex(be), integer())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(scanIndex(be), c(NA_integer_, NA_integer_))
    df$scanIndex <- c(1L, 2L)
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(scanIndex(be), c(1L, 2L))
})

test_that("smoothed, smoothed<-,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(smoothed(be), logical())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
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
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_null(spectraNames(be))
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    rownames(df) <- c("sp_1", "sp_2")
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(spectraNames(be), c("sp_1", "sp_2"))
    expect_error(spectraNames(be) <- "a", "rownames length")
    spectraNames(be) <- c("a", "b")
    expect_equal(spectraNames(be), c("a", "b"))
})

test_that("tic,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(tic(be), numeric())
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(tic(be), c(NA_real_, NA_real_))
    expect_equal(tic(be, initial = FALSE), c(0, 0))
    df$totIonCurrent <- c(5, 3)
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(tic(be), c(5, 3))
    expect_equal(tic(be, initial = FALSE), c(0, 0))
    df$intensity <- list(5:7, 1:4)
    df$mz <- list(1:3, 1:4)
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(tic(be, initial = FALSE), c(sum(5:7), sum(1:4)))
})

test_that("spectraVariables,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    expect_equal(spectraVariables(be), names(Spectra:::.SPECTRA_DATA_COLUMNS))
    df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(spectraVariables(be), names(Spectra:::.SPECTRA_DATA_COLUMNS))
    df$other_column <- 3
    be <- backendInitialize(be, files = NA_character_, spectraData = df)
    expect_equal(spectraVariables(be), c(names(Spectra:::.SPECTRA_DATA_COLUMNS),
                                         "other_column"))
})

test_that("spectraData, spectraData<-, MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)
    expect_equal(colnames(res), names(Spectra:::.SPECTRA_DATA_COLUMNS))

    df <- DataFrame(fromFile = c(1L, 1L), scanIndex = 1:2, a = "a", b = "b")
    be <- backendInitialize(be, files = NA_character_, spectraData = df)

    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(all(names(Spectra:::.SPECTRA_DATA_COLUMNS) %in% colnames(res)))
    expect_equal(res$a, c("a", "a"))
    expect_equal(res$b, c("b", "b"))

    res <- spectraData(be, "msLevel")
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), "msLevel")
    expect_equal(res$msLevel, c(NA_integer_, NA_integer_))

    res <- spectraData(be, c("mz"))
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), "mz")
    expect_equal(res$mz, list(numeric(), numeric()))

    res <- spectraData(be, c("a", "intensity"))
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), c("a", "intensity"))
    expect_equal(res$intensity, list(numeric(), numeric()))
    expect_equal(res$a, c("a", "a"))
})

test_that("show,MsBackendDataFrame works", {
    be <- MsBackendDataFrame()
    show(be)
    df <- DataFrame(fromFile = c(1L, 1L), rt = c(1.2, 1.3))
    be <- backendInitialize(be, files = NA_character_, df)
})
