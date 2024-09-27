test_that("MsBackend methods throw errors", {
    setClass("DummyBackend",
             contains = "MsBackend")
    dm <- new("DummyBackend")

    expect_error(backendInitialize(dm), "not implemented for")
    expect_error(backendMerge(dm), "implemented for")
    expect_error(export(dm), "does not support")
    expect_error(acquisitionNum(dm), "implemented for")
    expect_error(peaksData(dm), "implemented for")
    expect_error(centroided(dm), "implemented for")
    expect_error(centroided(dm) <- TRUE, "implemented for")
    expect_error(collisionEnergy(dm), "implemented for")
    expect_error(collisionEnergy(dm) <- 3, "implemented for")
    expect_error(dataOrigin(dm), "implemented for")
    expect_error(dataOrigin(dm) <- "a", "implemented for")
    expect_error(dataStorage(dm), "implemented for")
    expect_error(dataStorage(dm) <- "a", "implemented for")
    expect_error(filterAcquisitionNum(dm), "implemented for")
    expect_error(intensity(dm), "implemented for")
    expect_error(intensity(dm) <- list(), "implemented for")
    expect_error(ionCount(dm), "implemented for")
    expect_error(isCentroided(dm), "implemented for")
    expect_error(isEmpty(dm), "implemented for")
    expect_error(isolationWindowLowerMz(dm), "implemented for")
    expect_error(isolationWindowLowerMz(dm) <- 3, "implemented for")
    expect_error(isolationWindowTargetMz(dm), "implemented for")
    expect_error(isolationWindowTargetMz(dm) <- 3, "implemented for")
    expect_error(isolationWindowUpperMz(dm), "implemented for")
    expect_error(isolationWindowUpperMz(dm) <- 3, "implemented for")
    expect_error(length(dm), "implemented for")
    expect_error(msLevel(dm), "implemented for")
    expect_error(mz(dm), "implemented for")
    expect_error(mz(dm) <- list(), "implemented for")
    expect_error(lengths(dm), "implemented for")
    expect_error(polarity(dm), "implemented for")
    expect_error(polarity(dm) <- 0, "implemented for")
    expect_error(precScanNum(dm), "implemented for")
    expect_error(precursorCharge(dm), "implemented for")
    expect_error(precursorIntensity(dm), "implemented for")
    expect_error(precursorMz(dm), "implemented for")
    expect_error(peaksData(dm) <- list(), "implemented for")
    expect_error(rtime(dm), "implemented for")
    expect_error(rtime(dm) <- 1.2, "implemented for")
    expect_error(scanIndex(dm), "implemented for")
    expect_error(selectSpectraVariables(dm), "implemented for")
    expect_error(smoothed(dm), "implemented for")
    expect_error(smoothed(dm) <- TRUE, "implemented for")
    expect_error(spectraData(dm), "implemented for")
    expect_error(spectraData(dm) <- DataFrame(), "implemented for")
    expect_error(spectraNames(dm), "implemented for")
    expect_error(spectraNames(dm) <- "a", "implemented for")
    expect_error(spectraVariables(dm), "implemented for")
    expect_error(split(dm, f = 2), "implemented for")
    expect_error(tic(dm), "implemented for")
    expect_error(dm[1], "implemented for")
    expect_error(dm$a, "implemented for")
    expect_error(dm$a <- "a", "implemented for")
    expect_error(extractByIndex(dm, 1), "implemented for")
})

test_that("extractByIndex not implemented fallback", {
    ## Backends that don't implement a dedicated `extractByIndex` method should
    ## fall back to the [ method.
    setClass("DummyBackend",
             contains = "MsBackend",
             slots = c(d = "integer"))
    dm <- new("DummyBackend")
    expect_error(extractByIndex(dm, 1L), "'extractByIndex' not implemented")

    dm@d <- 1:4

    ## Have an implementation for [ but not extractByIndex:
    setMethod("[", "DummyBackend", function(x, i, j, ..., drop = FALSE) {
        x@d <- x@d[i]
        x
    })

    res <- dm[c(3, 1)]
    expect_equal(res@d, c(3L, 1L))

    res <- extractByIndex(dm, c(3, 1))
    expect_equal(res@d, c(3L, 1L))
})

test_that("reset,MsBackend works", {
    res <- reset(sciex_mzr)
    expect_equal(res, sciex_mzr)
})

test_that("backendBpparam,MsBackend works", {
    expect_equal(backendBpparam(sciex_mzr), bpparam())
    mcp <- MulticoreParam(2)
    res <- backendBpparam(sciex_mzr, mcp)
    expect_equal(res, mcp)
    res <- backendBpparam(sciex_mzr, SerialParam())
    expect_equal(res, SerialParam())
})

test_that("backendParallelFactor,MsBackend works", {
    expect_equal(backendParallelFactor(MsBackendMemory()), factor())
})

test_that("dataStorageBasePath,MsExperiment works", {
    expect_identical(dataStorageBasePath(MsBackendMemory()), NA_character_)
    tmp <- MsBackendMemory()
    expect_warning(dataStorageBasePath(tmp) <- "/", "not support")
})
