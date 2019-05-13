
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
    expect_equal(intensity(be), SimpleList())

    res <- intensity(sciex_mzr)
    expect_true(is(res, "SimpleList"))
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
    be <- MsBackendMzR()
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
    expect_equal(mz(be), SimpleList())

    res <- mz(sciex_mzr)
    expect_true(is(res, "SimpleList"))
    expect_true(is.numeric(res[[1]]))
    expect_true(!any(vapply(res, is.unsorted, logical(1))))
    expect_equal(length(res), length(sciex_mzr))
})

test_that("peaks,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(peaks(be), list())

    res <- peaks(sciex_mzr)
    expect_true(is(res, "list"))
    expect_equal(length(res), length(sciex_mzr))
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

test_that("spectraNames, spectraNames<-,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_null(spectraNames(be))

    tmp <- sciex_mzr
    res <- spectraNames(tmp)
    expect_equal(res, NULL)
    expect_error(spectraNames(tmp) <- "a", "rownames length")
    spectraNames(tmp) <- paste0("sp_", 1:length(tmp))
    expect_equal(rownames(tmp@spectraData), paste0("sp_", 1:length(tmp)))
    expect_equal(spectraNames(tmp), paste0("sp_", 1:length(tmp)))
})

test_that("tic,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(tic(be), numeric())

    res_initial <- tic(sciex_mzr)
    res_file <- tic(sciex_mzr, initial = FALSE)
    expect_true(is.numeric(res_initial))
    expect_true(is.numeric(res_file))
    expect_true(length(res_initial) == length(sciex_mzr))
    expect_true(length(res_file) == length(sciex_mzr))

    expect_equal(res_initial, res_file)
})

test_that("spectraVariables,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(spectraVariables(be), names(.SPECTRA_DATA_COLUMNS))

    res <- spectraVariables(sciex_mzr)
    expect_true(all(res %in% c(colnames(sciex_mzr@spectraData),
                               names(.SPECTRA_DATA_COLUMNS))))
})

test_that("spectraData, spectraData<-, MsBackendMzR works", {
    be <- MsBackendMzR()
    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)
    expect_true(all(names(.SPECTRA_DATA_COLUMNS) %in% colnames(res)))

    tmp <- sciex_mzr
    res <- spectraData(tmp)
    expect_true(all(names(.SPECTRA_DATA_COLUMNS) %in% colnames(res)))
    expect_true(all(colnames(tmp@spectraData) %in% colnames(res)))
    expect_true(is.logical(res$smoothed))
    expect_equal(nrow(res), length(tmp))

    spectraData(tmp)$new_col <- 1
    expect_true(any(colnames(tmp@spectraData) == "new_col"))
    expect_true(is(tmp@spectraData$new_col, "Rle"))
    expect_true(any(spectraVariables(tmp) == "new_col"))

    res <- spectraData(tmp, columns = c("msLevel", "new_col", "rtime"))
    expect_equal(colnames(res), c("msLevel", "new_col", "rtime"))
    expect_true(is.integer(res$msLevel))
    expect_true(is.numeric(res$new_col))
    expect_true(is.numeric(res$rtime))

    res <- spectraData(tmp, columns = c("msLevel", "new_col", "smoothed"))
    expect_equal(colnames(res), c("msLevel", "new_col", "smoothed"))
    expect_true(is.integer(res$msLevel))
    expect_true(is.numeric(res$new_col))
    expect_true(is.logical(res$smoothed))
    expect_true(all(is.na(res$smoothed)))

    spd <- spectraData(tmp, columns = c("msLevel", "rtime", "fromFile"))
    expect_error(spectraData(tmp) <- spd, "scanIndex")
    spd <- spectraData(tmp, columns = c("msLevel", "rtime", "fromFile",
                                        "scanIndex"))
    spectraData(tmp) <- spd
    expect_true(all(is.na(centroided(tmp))))
    expect_true(all(is.na(polarity(tmp))))
    expect_equal(mz(tmp), mz(sciex_mzr))
})

test_that("show,MsBackendMzR works", {
    be <- MsBackendMzR()
    show(be)

    show(sciex_mzr)
})

test_that("[,MsBackendMzR works", {
    tmp <- sciex_mzr[13:25]
    expect_true(validObject(tmp))
    expect_equal(length(tmp), 13)
    expect_equal(tmp@spectraData$scanIndex, 13:25)
    expect_true(all(is.na(smoothed(tmp))))

    ints <- intensity(tmp)
    spd <- spectraData(tmp)
    expect_equal(ints, spd$intensity)

    tmp <- sciex_mzr[1000]
    expect_equal(length(tmp), 1)
    spd <- spectraData(tmp)
    expect_equal(spd$mz, mz(tmp))
})

test_that("selectSpectraVariables,MsBackendMzR works", {
    be <- sciex_mzr

    res <- selectSpectraVariables(be, c("fromFile", "msLevel", "rtime",
                                        "scanIndex"))
    expect_equal(colnames(res@spectraData), c("fromFile", "msLevel", "rtime",
                                              "scanIndex"))
    expect_error(selectSpectraVariables(be, c("fromFile", "msLevel")),
                 "scanIndex is/are missing")
})

test_that("$,$<-,MsBackendDataFrame works", {
    tmp <- sciex_mzr
    tmp$new_col <- 5
    expect_true(any(spectraVariables(tmp) == "new_col"))
    expect_true(all(tmp$new_col == 5))
    expect_equal(rtime(tmp), rtime(sciex_mzr))
    expect_true(is.numeric(tmp$new_col))
    expect_true(is(tmp@spectraData$new_col, "Rle"))

    expect_error(tmp$mz <- SimpleList(1:4, 1:6), "not support replacing mz")
    expect_error(tmp$new_col <- c(2, 4), "either 1 or")
})
