test_that("initializeBackend,MsBackendMzR works", {
    fl <- normalizePath(
        dir(system.file("sciex", package = "msdata"), full.names = TRUE))
    expect_error(backendInitialize(MsBackendMzR()), "Parameter 'files'")
    expect_error(backendInitialize(MsBackendMzR(), files = 4),
                 "expected to be a character")
    be <- backendInitialize(MsBackendMzR(), files = fl)
    expect_true(validObject(be))
    expect_true(is(be, "MsBackendMzR"))
    expect_equal(unique(be$dataStorage), fl)
    expect_equal(nrow(be@spectraData), 1862)
    expect_equal(be@spectraData$scanIndex, c(1:931, 1:931))
    expect_equal(be@spectraData$dataStorage, rep(fl, each = 931))
    expect_true(isReadOnly(be))
})

test_that("backendMerge,MsBackendDataFrame works for MsBackendMzR too", {
    splt <- split(sciex_mzr, dataStorage(sciex_mzr))
    expect_equal(peaksData(splt[[1]]), sciex_pks[1:931])
    expect_equal(peaksData(splt[[2]]), sciex_pks[932:1862])
    res <- backendMerge(splt)
    expect_equal(res, sciex_mzr)

    res <- backendMerge(splt[2:1])
    dstor <- unique(dataStorage(res))
    expect_equal(unique(res$dataStorage), unique(sciex_mzr$dataStorage)[2:1])
    expect_equal(rtime(res)[dataStorage(res) == dstor[2]],
                 rtime(sciex_mzr)[dataStorage(sciex_mzr) == dstor[2]])
    expect_equal(rtime(res)[dataStorage(res) == dstor[1]],
                 rtime(sciex_mzr)[dataStorage(sciex_mzr) == dstor[1]])

    splt[[2]]@spectraData$some_col <- "a"
    res <- backendMerge(c(splt, tmt_mzr))
    expect_equal(unique(res$dataStorage), c(unique(sciex_mzr$dataStorage),
                                            unique(tmt_mzr$dataStorage)))
    expect_true(all(res@spectraData$some_col[dataStorage(res) == dstor[1]] == "a"))
    expect_true(all(is.na(res@spectraData$some_col[dataStorage(res) == dstor[2]])))
})

test_that("acquisitionNum,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(acquisitionNum(be), integer())
    expect_true(is(sciex_mzr@spectraData$acquisitionNum, "integer"))
    expect_equal(acquisitionNum(sciex_mzr), c(1:931, 1:931))
})

test_that("centroided, centroided<-, MsBackendMzR work", {
    be <- MsBackendMzR()
    expect_equal(centroided(be), logical())
    expect_true(is(centroided(sciex_mzr), "logical"))
    expect_true(is(sciex_mzr@spectraData$centroided, "logical"))
    expect_true(!all(centroided(sciex_mzr)))

    expect_error(centroided(sciex_mzr) <- "a", "to be a 'logical'")
    expect_error(centroided(sciex_mzr) <- c(FALSE, TRUE, TRUE), "has to be a")
    centroided(sciex_mzr) <- TRUE
    expect_true(is(centroided(sciex_mzr), "logical"))
    expect_true(is(sciex_mzr@spectraData$centroided, "logical"))
    expect_true(all(centroided(sciex_mzr)))

    centroided(sciex_mzr) <- rep(FALSE, length(sciex_mzr))
    expect_true(is(centroided(sciex_mzr), "logical"))
    expect_true(is(sciex_mzr@spectraData$centroided, "logical"))
    expect_true(!all(centroided(sciex_mzr)))
})

test_that("collisionEnergy, collisionEnergy<-,MsBackendMzR work", {
    be <- MsBackendMzR()
    expect_equal(collisionEnergy(be), numeric())

    expect_true(is(collisionEnergy(sciex_mzr), "numeric"))
    expect_true(!any(colnames(sciex_mzr@spectraData) == "collisionEnergy"))
    expect_true(all(is.na(collisionEnergy(sciex_mzr))))

    expect_error(collisionEnergy(sciex_mzr) <- "a", "to be a 'numeric'")
    expect_error(collisionEnergy(sciex_mzr) <- c(2.3), "has to be a")

    rn <- rnorm(length(sciex_mzr))
    collisionEnergy(sciex_mzr) <- rn
    expect_equal(collisionEnergy(sciex_mzr), rn)
})

test_that("dataStorage,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(dataStorage(be), character())

    expect_true(is(dataStorage(sciex_mzr), "character"))
    expect_true(is(sciex_mzr@spectraData$dataStorage, "character"))
    expect_equal(dataStorage(sciex_mzr), rep(sciex_file, each = 931))
})

test_that("intensity,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(intensity(be), NumericList(compress = FALSE))

    res <- intensity(sciex_mzr)
    expect_true(is(res, "NumericList"))
    expect_true(is.numeric(res[[1]]))
    expect_equal(length(res), length(sciex_mzr))
})

test_that("intensity<-,MsBackendMzR works", {
    be <- MsBackendMzR()

    expect_error(intensity(be) <- list, "does not support replacing intensity")
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

test_that("isolationWindowLowerMz,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_identical(isolationWindowLowerMz(be), numeric())

    be <- tmt_mzr
    expect_true(is.numeric(isolationWindowLowerMz(be)))
    expect_true(all(is.na(isolationWindowLowerMz(be)[msLevel(be) == 1])))
    expect_true(all(!is.na(isolationWindowLowerMz(be)[msLevel(be) == 2])))

    isolationWindowLowerMz(be) <- rep(2, length(be))
    expect_true(is(be@spectraData$isolationWindowLowerMz, "numeric"))
    expect_true(all(isolationWindowLowerMz(be) == 2))

    expect_error(isolationWindowLowerMz(be) <- 2, "of length 509")

    be <- sciex_mzr
    expect_true(all(is.na(isolationWindowLowerMz(be))))
})

test_that("isolationWindowTargetMz,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_identical(isolationWindowTargetMz(be), numeric())

    be <- tmt_mzr
    expect_true(is.numeric(isolationWindowTargetMz(be)))
    expect_true(all(is.na(isolationWindowTargetMz(be)[msLevel(be) == 1])))
    expect_true(all(!is.na(isolationWindowTargetMz(be)[msLevel(be) == 2])))
    expect_true(all(isolationWindowTargetMz(be)[msLevel(be) == 2] >
                    isolationWindowLowerMz(be)[msLevel(be) == 2]))

    isolationWindowTargetMz(be) <- rep(2, length(be))
    expect_true(is(be@spectraData$isolationWindowTargetMz, "numeric"))
    expect_true(all(isolationWindowTargetMz(be) == 2))

    expect_error(isolationWindowTargetMz(be) <- 2, "of length 509")

    be <- sciex_mzr
    expect_true(all(is.na(isolationWindowTargetMz(be))))
})

test_that("isolationWindowUpperMz,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_identical(isolationWindowUpperMz(be), numeric())

    be <- tmt_mzr
    expect_true(is.numeric(isolationWindowUpperMz(be)))
    expect_true(all(is.na(isolationWindowUpperMz(be)[msLevel(be) == 1])))
    expect_true(all(!is.na(isolationWindowUpperMz(be)[msLevel(be) == 2])))
    expect_true(all(isolationWindowUpperMz(be)[msLevel(be) == 2] >
                    isolationWindowTargetMz(be)[msLevel(be) == 2]))

    isolationWindowUpperMz(be) <- rep(2, length(be))
    expect_true(is(be@spectraData$isolationWindowUpperMz, "numeric"))
    expect_true(all(isolationWindowUpperMz(be) == 2))

    expect_error(isolationWindowUpperMz(be) <- 2, "of length 509")

    be <- sciex_mzr
    expect_true(all(is.na(isolationWindowUpperMz(be))))
})

test_that("msLevel,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(msLevel(be), integer())

    expect_true(is(msLevel(sciex_mzr), "integer"))
    expect_true(is(sciex_mzr@spectraData$msLevel, "integer"))
    expect_true(all(msLevel(sciex_mzr) == 1L))

    expect_true(sum(msLevel(tmt_mzr) == 2) == 451)
})

test_that("mz,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(mz(be), NumericList(compress = FALSE))

    res <- mz(sciex_mzr)
    expect_true(is(res, "NumericList"))
    expect_true(is.numeric(res[[1]]))
    expect_true(!any(vapply(res, is.unsorted, logical(1))))
    expect_equal(length(res), length(sciex_mzr))
})

test_that("mz<-,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_error(mz(be) <- list(), "does not support replacing")
})

test_that("peaksData,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(peaksData(be), list())

    res <- peaksData(sciex_mzr)
    expect_true(is(res, "list"))
    expect_equal(length(res), length(sciex_mzr))
    expect_true(is(res[[1]], "matrix"))
    expect_equal(colnames(res[[1]]), c("mz", "intensity"))

    tmp_one <- backendInitialize(MsBackendMzR(), sciex_file[1])
    res_one <- peaksData(tmp_one)
    expect_equal(res[1:length(res_one)], res_one)

    ## Arbitrary ordering.
    idx <- sample(1:length(sciex_mzr))
    be <- sciex_mzr[idx]
    pks <- peaksData(be)
    expect_identical(pks, sciex_pks[idx])
})

test_that("lengths,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(lengths(be), integer())

    res <- lengths(sciex_mzr)
    expect_true(is.integer(res))
    expect_true(length(res) == length(sciex_mzr))
})

test_that("polarity, polarity<- MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(polarity(be), integer())

    expect_true(is(polarity(sciex_mzr), "integer"))
    expect_true(is(sciex_mzr@spectraData$polarity, "integer"))
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
    expect_false(any(colnames(sciex_mzr@spectraData) == "precScanNum"))
    expect_true(all(is.na(precScanNum(sciex_mzr))))

    expect_true(is(tmt_mzr@spectraData$precScanNum, "integer"))
    expect_true(length(unique(precScanNum(tmt_mzr))) > 1)
})

test_that("precursorCharge,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(precursorCharge(be), integer())

    expect_true(is(precursorCharge(sciex_mzr), "integer"))
    expect_false(any(colnames(sciex_mzr@spectraData) == "precursorCharge"))
    expect_true(all(is.na(precursorCharge(sciex_mzr))))

    expect_true(is(precursorCharge(tmt_mzr), "integer"))
    expect_true(is(tmt_mzr@spectraData$precursorCharge, "integer"))
    expect_true(length(unique(precursorCharge(tmt_mzr))) > 1)
})

test_that("precursorIntensity,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(precursorIntensity(be), numeric())

    expect_true(is(precursorIntensity(sciex_mzr), "numeric"))
    expect_false(any(colnames(sciex_mzr@spectraData) == "precursorIntensity"))
    expect_true(all(is.na(precursorIntensity(sciex_mzr))))

    expect_true(is(precursorIntensity(tmt_mzr), "numeric"))
    expect_true(is(tmt_mzr@spectraData$precursorIntensity, "numeric"))
    expect_true(length(unique(precursorIntensity(tmt_mzr))) > 1)
})

test_that("precursorMz,MsBackendMzR works", {
    be <- MsBackendMzR()
    expect_equal(precursorMz(be), numeric())

    expect_true(is(precursorMz(sciex_mzr), "numeric"))
    expect_false(any(colnames(sciex_mzr@spectraData) == "precursorMz"))
    expect_true(all(is.na(precursorMz(sciex_mzr))))

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
    expect_true(is(sciex_mzr@spectraData$rtime, "numeric"))

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
    expect_true(is(sciex_mzr@spectraData$smoothed, "logical"))
    expect_true(all(smoothed(sciex_mzr)))

    smoothed(sciex_mzr) <- rep(FALSE, length(sciex_mzr))
    expect_true(is(smoothed(sciex_mzr), "logical"))
    expect_true(is(sciex_mzr@spectraData$smoothed, "logical"))
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
    res <- .spectra_data_mzR(tmp)
    expect_true(all(names(.SPECTRA_DATA_COLUMNS) %in% colnames(res)))
    expect_true(all(colnames(tmp@spectraData) %in% colnames(res)))
    expect_true(is.logical(res$smoothed))
    expect_equal(nrow(res), length(tmp))

    spectraData(tmp)$new_col <- 1
    expect_true(any(colnames(tmp@spectraData) == "new_col"))
    expect_true(is(tmp@spectraData$new_col, "numeric"))
    expect_true(any(spectraVariables(tmp) == "new_col"))

    res <- .spectra_data_mzR(tmp, columns = c("msLevel", "new_col", "rtime"))
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

    ## Empty object.
    tmp_sub <- tmp[integer()]
    expect_identical(tmp_sub$new_col, numeric())
    spd <- spectraData(tmp_sub, columns = c("msLevel", "rtime", "new_col"))
    expect_true(nrow(spd) == 0)
    expect_identical(spd$msLevel, integer())
    expect_identical(spd$rtime, numeric())
    expect_identical(spd$new_col, numeric())
    expect_identical(tmp_sub$mz, IRanges::NumericList(compress = FALSE))

    spd <- spectraData(tmp, columns = c("msLevel", "rtime", "dataStorage"))
    expect_error(spectraData(tmp) <- spd, "scanIndex")
    spd <- spectraData(tmp, columns = c("msLevel", "rtime", "dataStorage",
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
    expect_true(validObject(tmp))
    expect_equal(length(tmp), 1)
    spd <- spectraData(tmp)
    expect_equal(spd$mz, mz(tmp))
})

test_that("selectSpectraVariables,MsBackendMzR works", {
    be <- sciex_mzr

    res <- selectSpectraVariables(be, c("dataStorage", "msLevel", "rtime",
                                        "scanIndex"))
    expect_equal(colnames(res@spectraData), c("dataStorage", "msLevel", "rtime",
                                              "scanIndex"))
    expect_error(selectSpectraVariables(be, c("dataStorage", "msLevel")),
                 "scanIndex is/are missing")
})

test_that("$,$<-,MsBackendMzR works", {
    tmp <- sciex_mzr
    tmp$new_col <- 5
    expect_true(any(spectraVariables(tmp) == "new_col"))
    expect_true(all(tmp$new_col == 5))
    expect_equal(rtime(tmp), rtime(sciex_mzr))
    expect_true(is.numeric(tmp$new_col))
    expect_true(is(tmp@spectraData$new_col, "numeric"))

    expect_error(tmp$mz <- NumericList(1:4, 1:6, compress = FALSE),
                 "not support replacing mz")
    expect_error(tmp$new_col <- c(2, 4), "either 1 or")
})

test_that("export,MsBackendMzR works", {
    df <- DataFrame(msLevel = c(1L, 2L, 2L), scanIndex = 4:6)
    df$mz <- list(c(1.1, 1.3, 1.5), c(4.1, 5.1), c(1.6, 1.7, 1.8, 1.9))
    df$intensity <- list(c(45.1, 34, 12), c(234.4, 1333), c(42.1, 34.2, 65, 6))
    sps <- Spectra(df)

    fl <- tempfile()
    expect_error(export(sps, file = fl), "is required")
    expect_error(export(sps, MsBackendDataFrame(), fl), "does not support")
    export(sps, MsBackendMzR(), file = fl)
    res <- Spectra(backendInitialize(MsBackendMzR(), fl))
    expect_equal(mz(res), mz(sps))

    expect_error(export(sps, MsBackendMzR(), file = c(fl, fl)),
                 "of length 1")

    ## Export to multiple files
    fls <- c(tempfile(), tempfile())
    fls <- fls[c(1, 2, 1)]
    export(sps, MsBackendMzR(), fls)

    a <- Spectra(backendInitialize(MsBackendMzR(), fls[1]))
    expect_true(length(a) == 2)
    expect_equal(mz(a), mz(sps[c(1, 3)]))
    b <- Spectra(backendInitialize(MsBackendMzR(), fls[2]))
    expect_true(length(b) == 1)
    expect_equal(mz(b), mz(sps[2]))

    sps <- Spectra(sciex_mzr)
    sps <- filterRt(sps, c(100, 150))
    sps <- setBackend(sps, MsBackendDataFrame())
    fls <- c(tempfile(), tempfile())
    names(fls) <- unique(sps$dataOrigin)
    export(sps, MsBackendMzR(), file = fls[sps$dataOrigin], copy = TRUE)
    res <- Spectra(backendInitialize(MsBackendMzR(), fls))
    expect_equal(rtime(res), rtime(sps))
    expect_equal(mz(res), mz(sps))

    export(sps, MsBackendMzR(), file = fl)
    res <- Spectra(backendInitialize(MsBackendMzR(), fl))
    expect_equal(rtime(res), rtime(sps))
})

test_that("dropNaSpectraVariables works with MsBackendMzR", {
    res <- dropNaSpectraVariables(sciex_mzr)
    expect_equal(mz(res[1]), mz(sciex_mzr[1]))
    expect_true(length(spectraVariables(res)) <
                length(spectraVariables(sciex_mzr)))
})

test_that("supportsSetBackend,MsBackendMzR", {
    expect_false(supportsSetBackend(MsBackendMzR()))
})
