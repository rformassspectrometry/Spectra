fl <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
sciex_mzr <- backendInitialize(MsBackendMzR(), files = fl)
sciex_pks <- peaks(sciex_mzr)
fl <- dir(system.file("proteomics", package = "msdata"), full.names = TRUE)
tmt_mzr <- backendInitialize(MsBackendMzR(), files = fl[5])

test_that("initializeBackend,MsBackendMzR works", {
    fl <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
    expect_error(backendInitialize(MsBackendMzR()), "Parameter 'files'")
    expect_error(backendInitialize(MsBackendMzR(), files = c(fl, fl)),
                 "Duplicated")
    be <- backendInitialize(MsBackendMzR(), files = fl)
    expect_true(validObject(be))
    expect_true(is(be, "MsBackendMzR"))
    expect_equal(be@files, fl)
    expect_equal(be@modCount, c(0L, 0L))
    expect_equal(nrow(be@spectraData), 1862)
    expect_equal(be@spectraData$scanIndex, c(1:931, 1:931))
    expect_equal(be@spectraData$fromFile, Rle(rep(1:2, each = 931)))
    expect_true(isReadOnly(be))
})

test_that("acquisitionNum, MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(acquisitionNum(be), integer())
    expect_true(is(sciex_mzr@spectraData$acquisitionNum, "integer"))
    expect_equal(acquisitionNum(sciex_mzr), c(1:931, 1:931))
})

test_that("centroided, centroided<-, MsBackendMzR work", {
    be <- MsBackendMzR()
    expect_equal(centroided(be), logical())
    expect_true(is(centroided(sciex_mzr), "logical"))
    expect_true(is(sciex_mzr@spectraData$centroided, "Rle"))
    expect_true(!all(centroided(sciex_mzr)))

    expect_error(centroided(sciex_mzr) <- "a", "to be a 'logical'")
    expect_error(centroided(sciex_mzr) <- c(FALSE, TRUE, TRUE), "has to be a")
    centroided(sciex_mzr) <- TRUE
    expect_true(is(centroided(sciex_mzr), "logical"))
    expect_true(is(sciex_mzr@spectraData$centroided, "Rle"))
    expect_true(all(centroided(sciex_mzr)))

    centroided(sciex_mzr) <- rep(FALSE, length(sciex_mzr))
    expect_true(is(centroided(sciex_mzr), "logical"))
    expect_true(is(sciex_mzr@spectraData$centroided, "Rle"))
    expect_true(!all(centroided(sciex_mzr)))
})

test_that("collisionEnergy, collisionEnergy<-,MsBackendMzR work", {
    be <- MsBackendMzR()
    expect_equal(collisionEnergy(be), numeric())

    expect_true(is(collisionEnergy(sciex_mzr), "numeric"))
    expect_true(is(sciex_mzr@spectraData$collisionEnergy, "Rle"))
    expect_true(all(collisionEnergy(sciex_mzr) == 0))

    expect_error(collisionEnergy(sciex_mzr) <- "a", "to be a 'numeric'")
    expect_error(collisionEnergy(sciex_mzr) <- c(2.3), "has to be a")

    rn <- rnorm(length(sciex_mzr))
    collisionEnergy(sciex_mzr) <- rn
    expect_equal(collisionEnergy(sciex_mzr), rn)
})

test_that("fromFile,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(fromFile(be), integer())

    expect_true(is(fromFile(sciex_mzr), "integer"))
    expect_true(is(sciex_mzr@spectraData$fromFile, "Rle"))
    expect_equal(fromFile(sciex_mzr), rep(1:2, each = 931))
})

test_that("intensity,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(intensity(be), list())

    res <- intensity(sciex_mzr)
    expect_true(is(res, "list"))
    expect_true(is.numeric(res[[1]]))
    expect_equal(length(res), length(sciex_mzr))
})

test_that("ionCount,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(ionCount(be), numeric())

    res <- ionCount(sciex_mzr)
    expect_true(is.numeric(res))
    expect_true(length(res) == length(sciex_mzr))
    expect_equal(res, vapply(sciex_pks, function(z) sum(z[, 2]), numeric(1)))
})

test_that("isCentroided,MsBackendMzR works", {
    res <- MsBackendMzR()
    expect_equal(isCentroided(be), logical())

    res <- isCentroided(sciex_mzr)
    expect_true(is.logical(res))
    expect_true(length(res) == length(sciex_mzr))
    expect_true(all(!res))
})

test_that("isEmpty,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(isEmpty(be), logical())

    res <- isEmpty(sciex_mzr)
    expect_true(is.logical(res))
    expect_true(length(res) == length(sciex_mzr))
    expect_true(all(!res))
})

test_that("msLevel,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(msLevel(be), integer())

    expect_true(is(msLevel(sciex_mzr), "integer"))
    expect_true(is(sciex_mzr@spectraData$msLevel, "Rle"))
    expect_true(all(msLevel(sciex_mzr) == 1L))

    expect_true(sum(msLevel(tmt_mzr) == 2) == 451)
})

test_that("mz,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(mz(be), list())

    res <- mz(sciex_mzr)
    expect_true(is.list(res))
    expect_true(is.numeric(res[[1]]))
    expect_true(!any(vapply(res, is.unsorted, logical(1))))
    expect_equal(length(res), length(sciex_mzr))
})

test_that("peaks,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(peaks(be), list())

    res <- peaks(sciex_mzr)
    expect_true(is(res, "list"))
    expect_equal(length(res), length(object))
    expect_true(is(res[[1]], "matrix"))
    expect_equal(colnames(res[[1]]), c("mz", "intensity"))

    tmp_one <- backendInitialize(MsBackendMzR(), sciex_mzr@files[1])
    res_one <- peaks(tmp_one)
    expect_equal(res[1:length(res_one)], res_one)
})

test_that("peaksCount,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(peaksCount(be), integer())

    res <- peaksCount(sciex_mzr)
    expect_true(is.integer(res))
    expect_true(length(res) == length(sciex_mzr))
})

test_that("polarity, polarity<- MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(polarity(be), integer())

    expect_true(is(polarity(sciex_mzr), "integer"))
    expect_true(is(sciex_mzr@spectraData$polarity, "Rle"))
    expect_true(all(polarity(sciex_mzr) == 1L))

    expect_error(polarity(sciex_mzr) <- "a", "has to be an 'integer'")
    expect_error(polarity(sciex_mzr) <- c(1L, 1L, 2L), "has to be")

    polarity(sciex_mzr) <- 0
    expect_true(all(polarity(sciex_mzr) == 0L))
    polarity(sciex_mzr) <- 1:length(sciex_mzr)
    expect_equal(polarity(sciex_mzr), 1:length(sciex_mzr))
    expect_equal(sciex_mzr@spectraData$polarity, 1:length(sciex_mzr))
})

test_that("precScanNum,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(precScanNum(be), integer())

    expect_true(is(precScanNum(sciex_mzr), "integer"))
    expect_true(is(sciex_mzr@spectraData$precScanNum, "Rle"))
    expect_true(all(precScanNum(sciex_mzr) == 0L))

    expect_true(is(tmt_mzr@spectraData$precScanNum, "integer"))
    expect_true(length(unique(precScanNum(tmt_mzr))) > 1)
})

test_that("precursorCharge,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(precursorCharge(be), integer())

    expect_true(is(precursorCharge(sciex_mzr), "integer"))
    expect_true(is(sciex_mzr@spectraData$precursorCharge, "Rle"))
    expect_true(all(precursorCharge(sciex_mzr) == 0L))

    expect_true(is(precursorCharge(tmt_mzr), "integer"))
    expect_true(is(tmt_mzr@spectraData$precursorCharge, "integer"))
    expect_true(length(unique(precursorCharge(tmt_mzr))) > 1)
})

test_that("precursorIntensity,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(precursorIntensity(be), numeric())

    expect_true(is(precursorIntensity(sciex_mzr), "numeric"))
    expect_true(is(sciex_mzr@spectraData$precursorIntensity, "Rle"))
    expect_true(all(precursorIntensity(sciex_mzr) == 0L))

    expect_true(is(precursorIntensity(tmt_mzr), "numeric"))
    expect_true(is(tmt_mzr@spectraData$precursorIntensity, "numeric"))
    expect_true(length(unique(precursorIntensity(tmt_mzr))) > 1)
})

test_that("precursorMz,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(precursorMz(be), numeric())

    expect_true(is(precursorMz(sciex_mzr), "numeric"))
    expect_true(is(sciex_mzr@spectraData$precursorMz, "Rle"))
    expect_true(all(precursorMz(sciex_mzr) == 0L))

    expect_true(is(precursorMz(tmt_mzr), "numeric"))
    expect_true(is(tmt_mzr@spectraData$precursorMz, "numeric"))
    expect_true(length(unique(precursorMz(tmt_mzr))) > 1)
})

test_that("rtime, rtime<-,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(rtime(be), numeric())

    expect_true(is(rtime(sciex_mzr), "numeric"))
    expect_true(is(sciex_mzr@spectraData$rtime, "numeric"))
    expect_true(length(unique(rtime(sciex_mzr))) > 1)

    expect_error(rtime(sciex_mzr) <- "a", "'numeric' of length")
    expect_error(rtime(sciex_mzr) <- 0.2, "'numeric' of length")

    rts <- rtime(sciex_mzr)
    rtime(sciex_mzr) <- rep(0.1, length(sciex_mzr))
    expect_true(is(rtime(sciex_mzr), "numeric"))
    expect_true(is(sciex_mzr@spectraData$rtime, "Rle"))

    rtime(sciex_mzr) <- rts
    expect_equal(rtime(sciex_mzr), rts)
})

test_that("scanIndex,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(scanIndex(be), integer())

    expect_true(is(scanIndex(sciex_mzr), "integer"))
    expect_true(is(sciex_mzr@spectraData$scanIndex, "integer"))
    expect_true(length(unique(scanIndex(sciex_mzr))) > 1)
})

test_that("smoothed, smoothed<-,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(smoothed(be), logical())

    expect_true(is(smoothed(sciex_mzr), "logical"))
    expect_true(all(is.na(smoothed(sciex_mzr))))

    expect_error(smoothed(sciex_mzr) <- "2", "has to be a 'logical'")
    expect_error(smoothed(sciex_mzr) <- c(TRUE, TRUE, FALSE), "of length 1")

    smoothed(sciex_mzr) <- TRUE
    expect_true(is(smoothed(sciex_mzr), "logical"))
    expect_true(is(sciex_mzr@spectraData$smoothed, "Rle"))
    expect_true(all(smoothed(sciex_mzr)))

    smoothed(sciex_mzr) <- rep(FALSE, length(sciex_mzr))
    expect_true(is(smoothed(sciex_mzr), "logical"))
    expect_true(is(sciex_mzr@spectraData$smoothed, "Rle"))
    expect_true(all(!smoothed(sciex_mzr)))
})

## test_that("spectraNames, spectraNames<-,MsBackendMzR works", {
##     be <- MsBackendMzR()
##     expect_null(spectraNames(be))
##     df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
##     be <- backendInitialize(be, files = NA_character_, spectraData = df)
##     expect_null(spectraNames(be))
##     df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
##     rownames(df) <- c("sp_1", "sp_2")
##     be <- backendInitialize(be, files = NA_character_, spectraData = df)
##     expect_equal(spectraNames(be), c("sp_1", "sp_2"))
##     expect_error(spectraNames(be) <- "a", "rownames length")
##     spectraNames(be) <- c("a", "b")
##     expect_equal(spectraNames(be), c("a", "b"))
## })

## test_that("tic,MsBackendMzR works", {
##     be <- MsBackendMzR()
##     expect_equal(tic(be), numeric())
##     df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
##     be <- backendInitialize(be, files = NA_character_, spectraData = df)
##     expect_equal(tic(be), c(NA_real_, NA_real_))
##     expect_equal(tic(be, initial = FALSE), c(0, 0))
##     df$totIonCurrent <- c(5, 3)
##     be <- backendInitialize(be, files = NA_character_, spectraData = df)
##     expect_equal(tic(be), c(5, 3))
##     expect_equal(tic(be, initial = FALSE), c(0, 0))
##     df$intensity <- list(5:7, 1:4)
##     df$mz <- list(1:3, 1:4)
##     be <- backendInitialize(be, files = NA_character_, spectraData = df)
##     expect_equal(tic(be, initial = FALSE), c(sum(5:7), sum(1:4)))
## })

## test_that("spectraVariables,MsBackendMzR works", {
##     be <- MsBackendMzR()
##     expect_equal(spectraVariables(be), names(.SPECTRA_DATA_COLUMNS))
##     df <- DataFrame(fromFile = 1L, msLevel = c(1L, 2L))
##     be <- backendInitialize(be, files = NA_character_, spectraData = df)
##     expect_true(all(names(.SPECTRA_DATA_COLUMNS) %in% spectraVariables(be)))
##     df$other_column <- 3
##     be <- backendInitialize(be, files = NA_character_, spectraData = df)
##     expect_true(all(c(names(.SPECTRA_DATA_COLUMNS), "other_column") %in%
##                     spectraVariables(be)))
## })

## test_that("spectraData, spectraData<-, MsBackendMzR works", {
##     be <- MsBackendMzR()
##     res <- spectraData(be)
##     expect_true(is(res, "DataFrame"))
##     expect_true(nrow(res) == 0)
##     expect_true(all(names(.SPECTRA_DATA_COLUMNS) %in% colnames(res)))

##     df <- DataFrame(fromFile = c(1L, 1L), scanIndex = 1:2, a = "a", b = "b")
##     be <- backendInitialize(be, files = NA_character_, spectraData = df)

##     res <- spectraData(be)
##     expect_true(is(res, "DataFrame"))
##     expect_true(all(names(.SPECTRA_DATA_COLUMNS) %in% colnames(res)))
##     expect_equal(res$a, c("a", "a"))
##     expect_equal(res$b, c("b", "b"))

##     res <- spectraData(be, "msLevel")
##     expect_true(is(res, "DataFrame"))
##     expect_equal(colnames(res), "msLevel")
##     expect_equal(res$msLevel, c(NA_integer_, NA_integer_))

##     res <- spectraData(be, c("mz"))
##     expect_true(is(res, "DataFrame"))
##     expect_equal(colnames(res), "mz")
##     expect_equal(res$mz, list(numeric(), numeric()))

##     res <- spectraData(be, c("a", "intensity"))
##     expect_true(is(res, "DataFrame"))
##     expect_equal(colnames(res), c("a", "intensity"))
##     expect_equal(res$intensity, list(numeric(), numeric()))
##     expect_equal(res$a, c("a", "a"))

##     ## spectraData<-
##     df <- DataFrame(fromFile = 1L, msLevel = 1L, rtime = 1.2)
##     expect_error(spectraData(be) <- df, "with 2 rows")
##     df <- DataFrame(fromFile = c(1L, 1L), msLevel = 2L, rtime = c(1.2, 1.3))
##     res <- spectraData(be) <- df
##     expect_equal(rtime(be), c(1.2, 1.3))
## })

## test_that("show,MsBackendMzR works", {
##     be <- MsBackendMzR()
##     show(be)
##     df <- DataFrame(fromFile = c(1L, 1L), rt = c(1.2, 1.3))
##     be <- backendInitialize(be, files = NA_character_, df)
##     show(be)
## })
