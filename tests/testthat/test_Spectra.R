test_that("Spectra,DataFrame works", {
    df <- DataFrame()
    res <- Spectra(df)
    expect_true(validObject(res))
    expect_true(length(res) == 0)

    df <- DataFrame(msLevel = c(1L, 2L))
    res <- Spectra(df)
    expect_identical(msLevel(res), c(1L, 2L))
    expect_true(length(res) == 2)

    df$polarity <- "NEG"
    expect_error(Spectra(df), "wrong data type: polarity")
})

test_that("Spectra,missing works", {
    res <- Spectra()
    expect_true(length(res) == 0)

    be <- backendInitialize(MsBackendDataFrame(), DataFrame(msLevel = c(1L, 2L),
                                                            fromFile = 1L))
    res <- Spectra(backend = be)
    expect_true(length(res) == 2)
    expect_identical(msLevel(res), c(1L, 2L))
})

test_that("Spectra,MsBackend works", {
    res <- Spectra()
    expect_true(length(res) == 0)

    be <- backendInitialize(MsBackendDataFrame(), DataFrame(msLevel = c(1L, 2L),
                                                            fromFile = 1L))
    res <- Spectra(be)
    expect_true(length(res) == 2)
    expect_identical(msLevel(res), c(1L, 2L))
})

test_that("Spectra,character works", {
    res <- Spectra(sciex_file, backend = MsBackendMzR())
    expect_true(is(res@backend, "MsBackendMzR"))
    expect_equal(unique(res@backend$dataStorage), sciex_file)
    expect_identical(rtime(res), rtime(sciex_mzr))

    res_2 <- Spectra(sciex_file)
    expect_true(is(res@backend, "MsBackendDataFrame"))
    expect_identical(rtime(res), rtime(res_2))
})

test_that("setBackend,Spectra works", {
    df <- DataFrame(rtime = as.numeric(1:9),
                    fact = c(2L, 1L, 2L, 1L, 3L, 2L, 3L, 3L, 1L))
    sps <- Spectra(df)
    res <- setBackend(sps, MsBackendDataFrame())
    expect_true(ncol(sps@backend@spectraData) < ncol(res@backend@spectraData))
    expect_identical(sps@backend@spectraData$fact,
                     res@backend@spectraData$fact)
    expect_identical(rtime(res), rtime(sps))
    expect_identical(dataStorage(res), dataStorage(sps))
    expect_identical(dataOrigin(res), dataStorage(sps))

    ## Use a different factor.
    res <- setBackend(sps, MsBackendDataFrame(), f = df$fact)
    expect_true(ncol(sps@backend@spectraData) < ncol(res@backend@spectraData))
    expect_identical(as.vector(sps@backend@spectraData$fact),
                     as.vector(res@backend@spectraData$fact))
    expect_identical(dataStorage(res), dataStorage(sps))
    expect_identical(rtime(res), rtime(sps))

    ## switch from mzR to DataFrame
    sps <- Spectra(sciex_mzr)
    res <- setBackend(sps, MsBackendDataFrame())
    expect_identical(rtime(sps), rtime(res))
    expect_identical(mz(sps), mz(res))
    expect_true(is(res@backend@spectraData$msLevel, "Rle"))
    expect_true(is(sps@backend@spectraData$msLevel, "Rle"))
    expect_true(is.integer(res$msLevel))
    expect_identical(dataOrigin(res), dataStorage(sps))

    ## switch from DataFrame to hdf5
    tdir <- normalizePath(paste0(tempdir(), "/a"))
    res <- setBackend(sps, MsBackendHdf5Peaks(), hdf5path = tdir)
    expect_identical(rtime(sps), rtime(res))
    expect_identical(peaks(sps), peaks(res))
    expect_identical(dataOrigin(res), dataStorage(sps))

    ## from DataFrame to hdf5 providing file names - need to disable
    ## parallelization
    res <- setBackend(sps, MsBackendHdf5Peaks(),
                      files = c(tempfile(), tempfile()),
                      f = rep(1, length(sps)))
    expect_identical(rtime(sps), rtime(res))
    expect_identical(peaks(sps), peaks(res))

    ## errors:
    expect_error(setBackend(sps, MsBackendMzR()), "is read-only")
})

test_that("c,Spectra works", {
    df1 <- DataFrame(msLevel = c(1L, 1L, 1L))
    df1$mz <- list(c(1.1, 1.2), c(1.5), c(1.4, 1.5, 1.6))
    df1$intensity <- list(c(4.5, 23), 452.1, c(4.1, 342, 123))
    sp1 <- Spectra(df1)

    df2 <- DataFrame(msLevel = c(2L, 2L), rtime = c(1.2, 1.5))
    df2$mz <- list(1.5, 1.5)
    df2$intensity <- list(1234.1, 34.23)
    sp2 <- Spectra(df2)

    df3 <- DataFrame(msLevel = c(3L, 3L), other_col = "a")
    df3$mz <- list(c(1.4, 1.5, 1.6), c(1.8, 1.9))
    df3$intensity <- list(c(123.4, 12, 5), c(43.1, 5))
    sp3 <- Spectra(df3)

    df4 <- df3
    df4$mz <- NULL
    df4$intensity <- NULL
    sp4 <- Spectra(df4)

    res <- c(sp1, sp2, sp3)
    expect_true(is(res, "Spectra"))
    expect_equal(length(res), sum(nrow(df1), nrow(df2), nrow(df3)))
    expect_identical(msLevel(res), c(1L, 1L, 1L, 2L, 2L, 3L, 3L))
    expect_identical(res$other_col, c(NA, NA, NA, NA, NA, "a", "a"))
    expect_true(length(res@processingQueue) == 0)
    expect_true(length(res@processing) == 1)

    ## One Spectra without m/z and intensity
    res <- c(sp3, sp4)
    expect_true(is(res, "Spectra"))
    expect_identical(mz(res), NumericList(c(1.4, 1.5, 1.6), c(1.8, 1.9),
                                          numeric(), numeric(),
                                          compress = FALSE))
    expect_identical(msLevel(res), rep(3L, 4))
    expect_identical(intensity(res), NumericList(c(123.4, 12, 5), c(43.1, 5),
                                                 numeric(), numeric(),
                                                 compress = FALSE))
    res <- c(sp4, sp3)
    expect_true(is(res, "Spectra"))
    expect_identical(mz(res), NumericList(numeric(), numeric(),
                                          c(1.4, 1.5, 1.6), c(1.8, 1.9),
                                          compress = FALSE))
    expect_identical(msLevel(res), rep(3L, 4))
    expect_identical(intensity(res), NumericList(numeric(), numeric(),
                                                 c(123.4, 12, 5), c(43.1, 5),
                                                 compress = FALSE))

    ## Two Spectra without m/z and intensity
    res <- c(sp4, sp4)
    expect_true(is(res, "Spectra"))
    expect_identical(mz(res), NumericList(numeric(), numeric(), numeric(),
                                          numeric(), compress = FALSE))

    sp1@metadata <- list(version = "1.0.0", date = date())
    res <- c(sp1, sp2)
    expect_equal(res@metadata, sp1@metadata)

    sp1@processingQueue <- list(ProcessingStep(sum))
    expect_error(c(sp1, sp2), "with non-empty processing")

    ## Different backends
    s1 <- Spectra(sciex_mzr)
    s2 <- Spectra(sciex_hd5)
    expect_error(c(s1, s2), "backends of the same type")

    ## BackendMzR
    res <- c(Spectra(tmt_mzr), Spectra(sciex_mzr))
    expect_identical(msLevel(res), c(msLevel(tmt_mzr), msLevel(sciex_mzr)))
    expect_identical(msLevel(sciex_mzr), msLevel(res[dataStorage(res) %in%
                                                     sciex_file]))
    expect_identical(msLevel(tmt_mzr), msLevel(res[dataStorage(res) ==
                                                   dataStorage(tmt_mzr)[1]]))
})

test_that("acquisitionNum,Spectra works", {
    sps <- Spectra()
    res <- acquisitionNum(sps)
    expect_identical(res, integer())
    df <- DataFrame(msLevel = c(1L, 1L))
    sps <- Spectra(df)
    res <- acquisitionNum(sps)
    expect_identical(res, c(NA_integer_, NA_integer_))
    df$acquisitionNum <- 1:2
    sps <- Spectra(df)
    res <- acquisitionNum(sps)
    expect_identical(res, 1:2)
})

test_that("centroided,centroided<-,Spectra works", {
    sps <- Spectra()
    res <- centroided(sps)
    expect_identical(res, logical())
    df <- DataFrame(msLevel = c(1L, 2L, 1L), centroided = c(TRUE, FALSE, TRUE))
    sps <- Spectra(df)
    res <- centroided(sps)
    expect_identical(res, c(TRUE, FALSE, TRUE))
    centroided(sps) <- c(FALSE, TRUE, FALSE)
    res <- centroided(sps)
    expect_identical(res, c(FALSE, TRUE, FALSE))
    centroided(sps) <- FALSE
    expect_true(all(centroided(sps) == FALSE))
    expect_error(centroided(sps) <- 3, "'logical' of length")
    expect_error(centroided(sps) <- c(TRUE, FALSE, FALSE, TRUE), "length 1 or 3")
})

test_that("collisionEnergy,collisionEnergy<-,Spectra works", {
    sps <- Spectra()
    res <- collisionEnergy(sps)
    expect_identical(res, numeric())

    df <- DataFrame(msLevel = c(1L, 2L, 2L))
    sps <- Spectra(df)
    res <- collisionEnergy(sps)
    expect_true(all(is.na(res)))
    expect_equal(length(res), length(sps))

    collisionEnergy(sps) <- c(1.2, 1.4, 1.7)
    res <- collisionEnergy(sps)
    expect_identical(res, c(1.2, 1.4, 1.7))

    df$collisionEnergy <- c(3.4, 4.3, 2.3)
    sps <- Spectra(df)
    res <- collisionEnergy(sps)
    expect_identical(res, c(3.4, 4.3, 2.3))

    expect_error(collisionEnergy(sps) <- 4, "of length 3")
    expect_error(collisionEnergy(sps) <- c("a", "b", "c"), "'numeric'")
})

test_that("dataOrigin,Spectra works", {
    sps <- Spectra()
    res <- dataOrigin(sps)
    expect_identical(res, character())

    df <- DataFrame(msLevel = c(1L, 2L))
    sps <- Spectra(df)
    res <- dataOrigin(sps)
    expect_identical(res, rep(NA_character_, length(sps)))

    df <- DataFrame(msLevel = c(1L, 2L), dataOrigin = c("a", "b"))
    sps <- Spectra(df)
    res <- dataOrigin(sps)
    expect_identical(res, c("a", "b"))

    sps <- Spectra(backend = sciex_mzr)
    res <- dataOrigin(sps)
    expect_identical(res, dataStorage(sps))
})

test_that("dataStorage,Spectra works", {
    sps <- Spectra()
    res <- dataStorage(sps)
    expect_identical(res, character())

    df <- DataFrame(msLevel = c(1L, 2L))
    sps <- Spectra(df)
    res <- dataStorage(sps)
    expect_identical(res, rep("<memory>", 2))

    sps <- Spectra(sciex_mzr)
    res <- dataStorage(sps)
    expect_identical(res, rep(sciex_file, each = 931))
})

test_that("length,Spectra works", {
    sps <- Spectra()
    expect_equal(length(sps), 0)

    df <- DataFrame(msLevel = c(1L, 2L))
    sps <- Spectra(df)
    expect_true(length(sps) == 2)
})

test_that("intensity,Spectra works", {
    sps <- Spectra()
    res <- intensity(sps)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == 0)

    df <- DataFrame(msLevel = c(1L, 1L),
                    rtime = c(1.2, 1.4),
                    centroided = TRUE)
    df$intensity <- list(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                         c(9, 6, 0, 0, 3, 9, 4, 0, 0, 0))
    df$mz <- list(c(1:10), c(1:10))
    sps <- Spectra(df)
    res <- intensity(sps)
    expect_identical(res, NumericList(
                          c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                          c(9, 6, 0, 0, 3, 9, 4, 0, 0, 0),
                          compress = FALSE))
    sps <- clean(sps, all = TRUE)
    res <- intensity(sps)
    expect_identical(res, NumericList(c(1, 6, 3, 9, 1), c(9, 6, 3, 9, 4),
                                      compress = FALSE))
})

test_that("ionCount,Spectra works", {
    sps <- Spectra()
    res <- ionCount(sps)
    expect_identical(res, numeric())

    df <- DataFrame(msLevel = c(1L, 2L), centroided = TRUE)
    df$intensity <- list(c(5, 9, 3), c(9, 8, 2))
    df$mz <- list(1:3, 1:3)
    sps <- Spectra(df)

    res <- ionCount(sps)
    expect_identical(res, c(17, 19))

    sps <- removePeaks(sps, threshold = 4)
    res <- ionCount(sps)
    expect_identical(res, c(14, 17))
})

test_that("isCentroided,Spectra works", {
    sps <- Spectra()
    res <- isCentroided(sps)
    expect_identical(res, logical())

    df <- DataFrame(msLevel = c(1L, 1L))
    df$intensity <- list(c(5, 6, 1), c(5, 3, 1))
    df$mz <- list(1:3, 1:3)
    sps <- Spectra(df)

    res <- isCentroided(sps)
    expect_identical(res, c(NA, NA))

    sps <- Spectra(sciex_mzr)
    res <- isCentroided(sps)
    expect_true(length(res) == length(sps))
    expect_true(all(!res))
})

test_that("isEmpty,Spectra works", {
    sps <- Spectra()
    res <- isEmpty(sps)
    expect_identical(res, logical())

    df <- DataFrame(msLevel = c(2L, 2L), centroided = TRUE)
    df$intensity <- list(c(4, 6, 1), c(45, 2))
    df$mz <- list(1:3, 1:2)

    sps <- Spectra(df)
    res <- isEmpty(sps)
    expect_identical(res, c(FALSE, FALSE))

    sps <- removePeaks(sps, threshold = 100)
    res <- isEmpty(sps)
    expect_identical(res, c(FALSE, FALSE))

    sps <- clean(sps, all = TRUE)
    res <- isEmpty(sps)
    expect_identical(res, c(TRUE, TRUE))
})

test_that("isolationWindowLowerMz,Spectra works", {
    sps <- Spectra()
    expect_identical(isolationWindowLowerMz(sps), numeric())

    sps <- Spectra(sciex_mzr)
    expect_true(all(is.na(isolationWindowLowerMz(sps))))
    isolationWindowLowerMz(sps) <- as.numeric(1:length(sps))
    expect_identical(isolationWindowLowerMz(sps), as.numeric(1:length(sps)))

    sps <- Spectra(tmt_mzr)
    expect_true(all(is.na(isolationWindowLowerMz(sps)[msLevel(sps) == 1L])))
    expect_true(all(!is.na(isolationWindowLowerMz(sps)[msLevel(sps) == 2L])))
})

test_that("isolationWindowTargetMz,Spectra works", {
    sps <- Spectra()
    expect_identical(isolationWindowTargetMz(sps), numeric())

    sps <- Spectra(sciex_mzr)
    expect_true(all(is.na(isolationWindowTargetMz(sps))))
    isolationWindowTargetMz(sps) <- as.numeric(1:length(sps))
    expect_identical(isolationWindowTargetMz(sps), as.numeric(1:length(sps)))

    sps <- Spectra(tmt_mzr)
    expect_true(all(is.na(isolationWindowTargetMz(sps)[msLevel(sps) == 1L])))
    expect_true(all(!is.na(isolationWindowTargetMz(sps)[msLevel(sps) == 2L])))
    expect_true(all(isolationWindowTargetMz(sps)[msLevel(sps) == 2L] >
                    isolationWindowLowerMz(sps)[msLevel(sps) == 2L]))
})

test_that("isolationWindowUpperMz,Spectra works", {
    sps <- Spectra()
    expect_identical(isolationWindowUpperMz(sps), numeric())

    sps <- Spectra(sciex_mzr)
    expect_true(all(is.na(isolationWindowUpperMz(sps))))
    isolationWindowUpperMz(sps) <- as.numeric(1:length(sps))
    expect_identical(isolationWindowUpperMz(sps), as.numeric(1:length(sps)))

    sps <- Spectra(tmt_mzr)
    expect_true(all(is.na(isolationWindowUpperMz(sps)[msLevel(sps) == 1L])))
    expect_true(all(!is.na(isolationWindowUpperMz(sps)[msLevel(sps) == 2L])))
    expect_true(all(isolationWindowUpperMz(sps)[msLevel(sps) == 2L] >
                    isolationWindowTargetMz(sps)[msLevel(sps) == 2L]))
})

test_that("msLevel,Spectra works", {
    sps <- Spectra()
    expect_identical(msLevel(sps), integer())

    df <- DataFrame(msLevel = c(1L, 2L))
    sps <- Spectra(df)
    expect_identical(msLevel(sps), c(1L, 2L))
})

test_that("mz,Spectra works", {
    sps <- Spectra()
    res <- mz(sps)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == 0)

    df <- DataFrame(msLevel = c(1L, 1L),
                    rtime = c(1.2, 1.4),
                    centroided = TRUE)
    df$intensity <- list(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                         c(9, 6, 0, 0, 3, 9, 4, 0, 0, 0))
    df$mz <- list(c(1:10), c(1:10))
    sps <- Spectra(df)
    res <- mz(sps)
    expect_identical(res, NumericList(1:10, 1:10, compress = FALSE))
    sps <- clean(sps, all = TRUE)
    res <- mz(sps)
    expect_equal(res, NumericList(c(3, 4, 5, 8, 9), c(1, 2, 5, 6, 7),
                                  compress = FALSE))
})

test_that("peaks,Spectra works", {
    df <- DataFrame(msLevel = c(1L, 2L), fromFile = 1L)
    df$mz <- list(1:4, 1:5)
    df$intensity <- list(1:4, 1:5)
    be <- backendInitialize(MsBackendDataFrame(), file = NA_character_, df)
    sps <- Spectra(backend = be)
    res <- peaks(sps)
    expect_true(is(res, "SimpleList"))
    expect_equal(res[[1]][, 1], 1:4)
    expect_equal(res[[2]][, 1], 1:5)

    sps <- Spectra(backend = MsBackendDataFrame())
    res <- peaks(sps)
    expect_true(is(res, "SimpleList"))
    expect_true(length(res) == 0)
})

test_that("peaksCount,Spectra works", {
    sps <- Spectra()
    res <- peaksCount(sps)
    expect_identical(res, integer())

    df <- DataFrame(msLevel = c(2L, 2L), centroided = TRUE)
    df$intensity <- list(c(4, 6, 1), c(45, 2))
    df$mz <- list(1:3, 1:2)

    sps <- Spectra(df)
    res <- peaksCount(sps)
    expect_identical(res, c(3L, 2L))

    sps <- removePeaks(sps, threshold = 100)
    res <- peaksCount(sps)
    expect_identical(res, c(3L, 2L))

    sps <- clean(sps, all = TRUE)
    res <- peaksCount(sps)
    expect_identical(res, c(0L, 0L))
})

test_that("polarity,polarity<-,Spectra works", {
    sps <- Spectra()
    expect_identical(polarity(sps), integer())

    df <- DataFrame(msLevel = c(1L, 1L))
    sps <- Spectra(df)
    expect_identical(polarity(sps), c(NA_integer_, NA_integer_))
    polarity(sps) <- c(1L, -1L)
    expect_identical(polarity(sps), c(1L, -1L))
    polarity(sps) <- 2L
    expect_identical(polarity(sps), c(2L, 2L))

    expect_error(polarity(sps) <- c(1L, 2L, 1L), "of length 1 or 2")
    expect_error(polarity(sps) <- c("a", "b"), "an 'integer'")
})

test_that("precScanNum,Spectra works", {
    sps <- Spectra()
    expect_identical(precScanNum(sps), integer())

    sps <- Spectra(DataFrame(msLevel = c(2L, 2L)))
    expect_identical(precScanNum(sps), c(NA_integer_, NA_integer_))
    sps <- Spectra(DataFrame(msLevel = c(2L, 2L), precScanNum = c(1L, 2L)))
    expect_identical(precScanNum(sps), c(1L, 2L))
})

test_that("precursorCharge,Spectra works", {
    sps <- Spectra()
    expect_identical(precursorCharge(sps), integer())

    sps <- Spectra(DataFrame(msLevel = c(2L, 2L)))
    expect_identical(precursorCharge(sps), c(NA_integer_, NA_integer_))
    sps <- Spectra(DataFrame(msLevel = c(2L, 2L), precursorCharge = c(1L, 3L)))
    expect_identical(precursorCharge(sps), c(1L, 3L))
})

test_that("precursorIntensity,Spectra works", {
    sps <- Spectra()
    expect_identical(precursorIntensity(sps), numeric())

    sps <- Spectra(DataFrame(msLevel = c(2L, 2L)))
    expect_identical(precursorIntensity(sps), c(NA_real_, NA_real_))
    sps <- Spectra(DataFrame(msLevel = c(2L, 2L),
                             precursorIntensity = c(1.3, 3.2)))
    expect_identical(precursorIntensity(sps), c(1.3, 3.2))
})

test_that("precursorMz,Spectra works", {
    sps <- Spectra()
    expect_identical(precursorMz(sps), numeric())

    sps <- Spectra(DataFrame(msLevel = c(2L, 2L)))
    expect_identical(precursorMz(sps), c(NA_real_, NA_real_))
    sps <- Spectra(DataFrame(msLevel = c(2L, 2L),
                             precursorMz = c(234.2, 668.2)))
    expect_identical(precursorMz(sps), c(234.2, 668.2))
})

test_that("rtime,rtime<-,Spectra works", {
    sps <- Spectra()
    expect_identical(rtime(sps), numeric())

    sps <- Spectra(DataFrame(msLevel = c(1L, 2L)))
    expect_identical(rtime(sps), c(NA_real_, NA_real_))
    sps <- Spectra(DataFrame(msLevel = c(1L, 2L), rtime = c(2.1, 4.2)))
    expect_identical(rtime(sps), c(2.1, 4.2))
    rtime(sps) <- c(1.2, 1.3)
    expect_identical(rtime(sps), c(1.2, 1.3))

    expect_error(rtime(sps) <- c(TRUE, FALSE), "'numeric' of length 2")
    expect_error(rtime(sps) <- 1.3, "'numeric' of length 2")
})

test_that("scanIndex,Spectra works", {
    sps <- Spectra()
    expect_identical(scanIndex(sps), integer())
    sps <- Spectra(DataFrame(msLevel = c(1L, 2L)))
    expect_identical(scanIndex(sps), c(NA_integer_, NA_integer_))

    sps <- Spectra(DataFrame(msLevel = c(1L, 2L), scanIndex = c(1L, 2L)))
    expect_identical(scanIndex(sps), c(1L, 2L))
})

test_that("selectSpectraVariables,Spectra works", {
    sps <- Spectra()
    res <- selectSpectraVariables(sps, c("msLevel", "rtime"))
    expect_equal(sps, res)

    sps <- Spectra(DataFrame(msLevel = c(1L, 2L), rtime = c(1.2, 1.4),
                             other_col = "a"))
    res <- selectSpectraVariables(sps, c("rtime", "other_col"))
    expect_equal(msLevel(res), c(NA_integer_, NA_integer_))
    expect_equal(rtime(res), c(1.2, 1.4))
    expect_equal(spectraData(res, "other_col")[, 1], c("a", "a"))

    res <- selectSpectraVariables(sps, c("msLevel"))
    expect_identical(rtime(res), c(NA_real_, NA_real_))
    expect_identical(msLevel(res), 1:2)
    expect_true(!any(spectraVariables(res) %in% "other_col"))

    expect_error(selectSpectraVariables(sps, c("rtime", "something")), "not")
})

test_that("smoothed,smoothed<-,Spectra works", {
    sps <- Spectra()
    expect_identical(smoothed(sps), logical())

    sps <- Spectra(DataFrame(msLevel = 1:2))
    expect_identical(smoothed(sps), c(NA, NA))
    sps <- Spectra(DataFrame(msLevel = 1:2, smoothed = TRUE))
    expect_identical(smoothed(sps), c(TRUE, TRUE))

    smoothed(sps) <- c(FALSE, TRUE)
    expect_identical(smoothed(sps), c(FALSE, TRUE))
    smoothed(sps) <- FALSE
    expect_identical(smoothed(sps), c(FALSE, FALSE))

    expect_error(smoothed(sps) <- c("a", "b"), "of length 1 or 2")
    expect_error(smoothed(sps) <- c(TRUE, FALSE, TRUE), "of length 1 or 2")
})

test_that("spectraData,Spectra works", {
    sps <- Spectra()
    res <- spectraData(sps)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)
    expect_true(all(c("msLevel", "rtime", "mz", "intensity") %in% colnames(res)))

    df <- DataFrame(msLevel = c(1L, 1L), rtime = c(1.2, 1.3), centroided = TRUE)
    df$mz <- list(1:10, 1:10)
    df$intensity <- list(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                         c(9, 6, 0, 0, 3, 0, 0, 0, 3, 2))
    sps <- Spectra(df)
    res <- spectraData(sps)
    expect_true(is(res, "DataFrame"))
    expect_equal(res$msLevel, c(1L, 1L))
    expect_equal(res$rtime, c(1.2, 1.3))
    expect_equal(res$precScanNum, c(NA_integer_, NA_integer_))
    expect_equal(res$mz, NumericList(1:10, 1:10, compress = FALSE))
    expect_equal(res$intensity, NumericList(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                                            c(9, 6, 0, 0, 3, 0, 0, 0, 3, 2),
                                            compress = FALSE))
    sps <- clean(sps, all = TRUE)
    res <- spectraData(sps)
    expect_equal(res$mz, NumericList(c(3, 4, 5, 8, 9), c(1, 2, 5, 9, 10),
                                     compress = FALSE))
    expect_equal(res$intensity, NumericList(c(1, 6, 3, 9, 1), c(9, 6, 3, 3, 2),
                                            compress = FALSE))
})

test_that("spectraData<-,Spectra works", {
    df <- DataFrame(msLevel = c(1L, 1L), rtime = c(1.2, 1.3), centroided = TRUE)
    df$mz <- list(1:10, 1:10)
    df$intensity <- list(c(0, 0, 1, 6, 3, 0, 0, 9, 1, 0),
                         c(9, 6, 0, 0, 3, 0, 0, 0, 3, 2))
    sps <- Spectra(df)

    spectraData(sps)$add_col <- 3
    expect_true(any(spectraVariables(sps) == "add_col"))
    expect_equal(spectraData(sps)$add_col, c(3, 3))

    df$mz <- list(11:20, 11:20)
    df$precursorMz <- c(0, 0)

    spectraData(sps) <- df
    expect_true(!any(spectraVariables(sps) == "add_col"))
    expect_equal(mz(sps), NumericList(11:20, 11:20, compress = FALSE))

    expect_error(spectraData(sps) <- c(1, 2, 4, 5), "has to be 2")

    expect_error(spectraData(sps) <- df[c(1, 1, 2), ], "with 2 rows")

    tmp <- sciex_mzr
    sps <- Spectra(tmp)
    expect_warning(spectraData(sps)$some_col <- "yes")
    expect_true(any(spectraVariables(sps) == "some_col"))
    expect_true(all(spectraData(sps, "some_col")[, 1] == "yes"))
    expect_true(is(sps@backend@spectraData$some_col, "Rle"))
    sps$other_col <- "other_value"
    expect_true(any(spectraVariables(sps) == "other_col"))
    expect_identical(sps$other_col, rep("other_value", length(sps)))

    spd <- spectraData(sps)
    spd$msLevel <- 2L
    spectraData(sps) <- spd
    expect_true(all(msLevel(sps) == 2L))
})

test_that("spectraNames,spectraNames<-,Spectra works", {
    sps <- Spectra()
    expect_identical(spectraNames(sps), NULL)

    sps <- Spectra(DataFrame(msLevel = c(1L, 1L)))
    expect_identical(spectraNames(sps), NULL)

    spectraNames(sps) <- c("a", "b")
    expect_identical(spectraNames(sps), c("a", "b"))

    expect_error(spectraNames(sps) <- 1, "invalid rownames length")
})

test_that("spectraVariables,Spectra works", {
    sps <- Spectra()
    res <- spectraVariables(sps)
    exp_col <- c("msLevel", "rtime", "acquisitionNum", "scanIndex", "mz",
                 "intensity", "dataStorage", "centroided", "smoothed",
                 "polarity", "precScanNum", "precursorMz", "precursorIntensity",
                 "precursorCharge", "collisionEnergy")
    expect_true(all(exp_col %in% res))
    df <- DataFrame(other_col = "a", msLevel = c(1L, 1L))
    sps <- Spectra(df)
    res <- spectraVariables(sps)
    expect_true(all(c(exp_col, "other_col") %in% res))
})

test_that("tic,Spectra works", {
    sps <- Spectra()
    expect_identical(tic(sps), numeric())

    sps <- Spectra(DataFrame(msLevel = c(1L, 1L)))
    expect_identical(tic(sps), c(NA_real_, NA_real_))

    df <- DataFrame(msLevel = c(1L, 1L), totIonCurrent = c(5, 3))
    sps <- Spectra(df)
    expect_identical(tic(sps), c(5, 3))
    expect_identical(tic(sps, initial = FALSE), c(0, 0))

    df$intensity <- list(c(3, 3, 1), c(5, 3, 1))
    df$mz <- list(1:3, 1:3)
    sps <- Spectra(df)
    expect_identical(tic(sps, initial = FALSE), c(7, 9))
})

test_that("$, $<-,Spectra works", {
    sps <- Spectra()
    expect_identical(sps$msLevel, integer())

    sps <- Spectra(DataFrame(msLevel = c(1L, 2L), other_column = "a"))
    expect_identical(sps$msLevel, c(1L, 2L))
    expect_identical(sps$other_column, c("a", "a"))

    sps$second_col <- c(1, 4)
    expect_identical(sps$second_col, c(1, 4))
    expect_true(any(spectraVariables(sps) == "second_col"))

    df <- DataFrame(msLevel = c(1L, 2L), other_column = "a")
    df$mz <- list(1:3, 1:4)
    df$intensity <- list(c(3, 6, 3), c(56, 6, 3, 2))
    sps <- Spectra(df)

    sps$intensity <- list(c(4, 4, 4), c(2, 4, 6, 3))
    expect_equal(intensity(sps), NumericList(c(4, 4, 4), c(2, 4, 6, 3),
                                             compress = FALSE))

    sps <- Spectra(sciex_mzr)
    sps$add_col <- "something"
    expect_true(all(sps$add_col == "something"))

    expect_error(sps$mz <- mz(sps))

    expect_error(sps$add_col <- c(1, 2), "has to be either 1 or")
})

#### ---------------------------------------------------------------------------
##
##                      FILTERING AND SUBSETTING
##
#### ---------------------------------------------------------------------------

test_that("[,Spectra works", {
    df <- DataFrame(msLevel = c(1L, 2L, 1L, 2L),
                    rtime = c(1, 2, 1, 2))
    sps <- Spectra(df)
    res <- sps[c(2, 4), ]
    expect_identical(msLevel(res), c(2L, 2L))
    expect_identical(sps, sps[])

    expect_error(sps[, 4], "columns is not")

    res <- sps[integer()]
    expect_true(is(res, "Spectra"))
    expect_true(length(res) == 0)

    sps <- Spectra(sciex_mzr)
    tmp <- sps[dataStorage(sps) == sciex_file[2], ]
    expect_true(all(dataStorage(tmp) == sciex_file[2]))
    expect_equal(unique(tmp$dataStorage), sciex_file[2])
    expect_equal(rtime(tmp), rtime(sps)[dataStorage(sps) == sciex_file[2]])
})

test_that("filterAcquisitionNum,Spectra works", {
    sps <- Spectra()
    res <- filterAcquisitionNum(sps, n = 3L)
    expect_equal(length(res), 0)
    expect_equal(length(res@processing), 1)

    sps <- Spectra(sciex_mzr)
    res <- filterAcquisitionNum(sps, n = 1:10, dataStorage = sciex_file[2])
    expect_equal(acquisitionNum(res),
                 c(1:sum(dataStorage(sps) == sciex_file[1]), 1:10))
})

test_that("filterDataOrigin,Spectra works", {
    sps <- Spectra()
    res <- filterDataOrigin(sps)
    expect_true(length(res) == 0)
    expect_true(length(res@processing) == 1)

    sps <- Spectra(sciex_mzr)
    res <- filterDataOrigin(sps, dataOrigin = "2")
    expect_true(length(res) == 0)
    res <- filterDataOrigin(sps, sciex_file[1])
    expect_identical(rtime(res), rtime(sps)[1:931])
    expect_true(length(res@processing) > length(sps@processing))

    dorig <- rep(letters[1:7], each = length(sps)/7)
    dataOrigin(sps) <- dorig
    res <- filterDataOrigin(sps, dataOrigin = c("d", "a"))
    expect_equal(unique(dataOrigin(res)), c("d", "a"))
    expect_equal(rtime(res)[1:266], rtime(sps)[sps$dataOrigin == "d"])
    expect_equal(peaks(res)[1:266], SimpleList(sciex_pks[sps$dataOrigin == "d"]))
})

test_that("filterDataStorage,Spectra works", {
    sps <- Spectra()
    res <- filterDataStorage(sps)
    expect_true(length(res) == 0)
    expect_true(length(res@processing) == 1)

    sps <- Spectra(sciex_mzr)
    res <- filterDataStorage(sps, "2")
    expect_true(length(res) == 0)
    res <- filterDataStorage(sps, sciex_file[2])
    expect_identical(rtime(res), rtime(sps)[dataStorage(sps) == sciex_file[2]])
    expect_true(length(res@processing) > length(sps@processing))
    expect_identical(peaks(res),
                     SimpleList(sciex_pks[dataStorage(sps) == sciex_file[2]]))
})

test_that("filterEmptySpectra,Spectra works", {
    sps <- Spectra()
    res <- filterEmptySpectra(sps)
    expect_true(length(res) == 0)
    expect_true(length(res@processing) == 1)

    df <- DataFrame(msLevel = c(1L, 2L, 1L, 2L),
                    rtime = c(1, 2, 1, 2), centroided = TRUE)
    df$mz <- SimpleList(1:3, integer(), 1:5, integer())
    df$intensity <- SimpleList(c(4, 6, 3), numeric(), c(34, 2, 5, 9, 9), numeric())
    sps <- Spectra(df)

    res <- filterEmptySpectra(sps)
    expect_equal(length(res), 2)
    expect_equal(length(res@processing), 1)
    expect_equal(rtime(res), c(1, 1))

    sps <- clean(removePeaks(sps, threshold = 20), all = TRUE)
    res <- filterEmptySpectra(sps)
    expect_equal(length(res), 1)
    expect_equal(rtime(res), 1)
    expect_equal(length(res@processing), 3)

    sps <- clean(removePeaks(sps, threshold = 50), all = TRUE)
    res <- filterEmptySpectra(sps)
    expect_equal(length(res), 0)

    sps <- Spectra(sciex_mzr)
    res <- filterEmptySpectra(sps)
    expect_equal(rtime(res), rtime(sps))
})

test_that("filterIsolationWindow,Spectra works", {
    sps <- Spectra()
    res <- filterIsolationWindow(sps)
    expect_true(is(res, "Spectra"))
    expect_true(length(res@processing) == 1)

    sps <- Spectra(sciex_mzr)
    res <- filterIsolationWindow(sps, 123.323)
    expect_true(length(res) == 0)

    sps <- Spectra(tmt_mzr)
    res <- filterIsolationWindow(sps, 544)
    expect_true(length(res) == 3)
    expect_true(all(isolationWindowLowerMz(res) < 544))
    expect_true(all(isolationWindowUpperMz(res) > 544))
    expect_true(all(precursorMz(res) < 545 & precursorMz(res) > 543))
})

test_that("filterMsLevel,Spectra works", {
    sps <- Spectra()
    res <- filterMsLevel(sps)
    expect_true(length(res) == 0)
    expect_true(length(res@processing) == 1)

    sps <- Spectra(sciex_mzr)
    res <- filterMsLevel(sps, 2L)
    expect_true(length(res) == 0)

    sps <- Spectra(tmt_mzr)
    expect_true(all(1:2 %in% msLevel(sps)))
    res <- filterMsLevel(sps, 1L)
    expect_false(all(1:2 %in% msLevel(res)))
    expect_true(all(msLevel(res) == 1))

    res <- filterMsLevel(sps, c(2L, 4L))
    expect_false(all(1:2 %in% msLevel(res)))
    expect_true(all(msLevel(res) == 2))
})

test_that("filterPolarity,Spectra works", {
    sps <- Spectra()
    res <- filterPolarity(sps)
    expect_true(length(res) == 0)
    expect_true(length(res@processing) == 1)

    sps <- Spectra(sciex_mzr)
    res <- filterPolarity(sps, polarity = c(2, 0))
    expect_true(length(res) == 0)

    res <- filterPolarity(sps, 1)
    expect_true(all(polarity(res) == 1))
    expect_equal(rtime(res), rtime(sps))
    expect_true(length(res@processing) == 1)
})

test_that("filterPrecursorMz,Spectra works", {
    sps <- Spectra()
    res <- filterPrecursorMz(sps)
    expect_true(is(res, "Spectra"))
    expect_true(length(res@processing) == 1)

    sps <- Spectra(tmt_mzr)
    res <- filterPrecursorMz(sps, mz = 544.75)
    expect_true(length(res) == 0)
    res <- filterPrecursorMz(sps, mz = 544.75, ppm = 40)
    expect_true(length(res) == 2)
})

test_that("filterPrecursorScan,Spectra works", {
    sps <- Spectra()
    res <- filterPrecursorScan(sps, 3)
    expect_true(length(res) == 0)
    expect_true(is(res, "Spectra"))
    expect_true(length(res@processing) == 1)

    sps <- Spectra(tmt_mzr)
    res <- filterPrecursorScan(sps, c(1087L, 1214L))
    expect_true(sum(msLevel(res) == 1) == 2)
    expect_true(all(c(1087L, 1214L) %in% acquisitionNum(res)))
})

test_that("filterRt,Spectra works", {
    sps <- Spectra()
    res <- filterRt(sps, c(1, 2))
    expect_true(is(res, "Spectra"))
    expect_true(length(res) == 0)
    expect_true(length(res@processing) == 1)

    sps <- Spectra(sciex_mzr)
    res <- filterRt(sps, rt = c(100, 120))
    expect_true(all(rtime(res) >= 100 & rtime(res) <= 120))
    expect_error(filterRt(sps, rt = c(100)))
    expect_error(filterRt(sps, rt = c(120, 100)))
    expect_error(filterRt(sps, rt = c("100", "120")))
    expect_error(filterRt(sps, rt = c(100, 120), msLevel. = "1"))

    res <- filterRt(sps, rt = c(100, 120), msLevel = 2L)
    expect_equal(rtime(res), rtime(sps))
})

#### ---------------------------------------------------------------------------
##
##                      DATA MANIPULATION METHODS
##
#### ---------------------------------------------------------------------------

test_that("bin,Spectra works", {
    sps <- Spectra(tmt_mzr)
    pks <- peaks(sps)
    res <- bin(sps, binSize = 2)
    expect_true(length(res@processingQueue) == 1)
    res1 <- bin(sps, msLevel = 1, binSize = 2)

    expect_identical(peaks(res1)[res1$msLevel == 2],
                     pks[sps$msLevel == 2])

    mzr <- range(unlist(mz(sps)))
    brks <- seq(floor(mzr[1]), ceiling(mzr[2]), by = 2)
    res1 <- bin(sps, msLevel = 1, binSize = 2, breaks = brks)
    res1_pks <- peaks(res1)
    res_pks <- peaks(res)
    expect_identical(res1_pks[res1$msLevel == 1],
                     res_pks[res$msLevel == 1])
    expect_true(all(lengths(res_pks) != lengths(pks)))

    expect_warning(res <- bin(sps, msLevel = 3))
    expect_identical(res, sps)
})

test_that("clean,Spectra works", {
    sps <- Spectra()
    res <- clean(sps)
    expect_true(length(res@processingQueue) == 1)
    expect_equal(res@processingQueue[[1]],
                 ProcessingStep(.peaks_clean,
                                list(all = FALSE, msLevel = integer())))

    res <- clean(sps, all = TRUE, msLevel = 2L)
    expect_true(length(res@processingQueue) == 1)
    expect_equal(res@processingQueue[[1]],
                 ProcessingStep(.peaks_clean,
                                list(all = TRUE, msLevel = 2L)))

    res <- clean(Spectra(sciex_mzr), all = TRUE)
    pks_res <- lapply(sciex_pks, .peaks_clean, all = TRUE,
                      spectrumMsLevel = 1L)
    expect_identical(peaks(res), SimpleList(pks_res))
})

test_that("pickPeaks,Spectra works", {
    sps <- Spectra()
    expect_error(pickPeaks(sps, halfWindowSize = 1), "integer")
    expect_error(pickPeaks(sps, halfWindowSize = 1L:2L), "length 1")
    expect_error(pickPeaks(sps, halfWindowSize = -1L), "> 0")
    expect_error(pickPeaks(sps, method = "foo"), "MAD")
    expect_error(pickPeaks(sps, snr = "foo"), "numeric")
    expect_error(pickPeaks(sps, snr = 1L:2L), "length 1")
    expect_error(pickPeaks(sps, snr = -1L), ">= 0")
    expect_error(pickPeaks(sps, k = 1), "integer")
    expect_error(pickPeaks(sps, k = 1L:2L), "length 1")
    expect_error(pickPeaks(sps, k = -1L), ">= 0")
    expect_error(pickPeaks(sps, descending = NA), "TRUE or FALSE")
    expect_error(pickPeaks(sps, descending = c(TRUE, TRUE)), "TRUE or FALSE")
    expect_error(pickPeaks(sps, threshold = "foo"), "numeric")
    expect_error(pickPeaks(sps, threshold = 1L:2L), "length 1")
    expect_error(pickPeaks(sps, threshold = -1L), ">= 0")
    expect_error(pickPeaks(sps, threshold = 2L), "<= 1")

    res <- pickPeaks(sps)
    expect_true(length(res@processingQueue) == 1)
    expect_equal(res@processingQueue[[1]],
                 ProcessingStep(.peaks_pick,
                                list(halfWindowSize = 2L, method = "MAD",
                                     snr = 0, k = 0L, descending = FALSE,
                                     threshold = 0L, msLevel = integer())))
    expect_match(res@processing,
                 "Peak picking with MAD noise estimation, hws = 2, snr = 0 \\[")
    res <- pickPeaks(sps, k = 2L)
    expect_match(res@processing,
                 paste0("Peak picking with MAD noise estimation, hws = 2, ",
                        "snr = 0 and centroid refinement \\["))

    sps <- Spectra(sciex_mzr)
    expect_equal(pickPeaks(sps, msLevel. = 3), sps)

    res <- pickPeaks(sps)
    pks_res <- lapply(sciex_pks, .peaks_pick, spectrumMsLevel = 1L,
                      centroided = FALSE)
    expect_identical(peaks(res), SimpleList(pks_res))
})

test_that("removePeaks,Spectra works", {
    sps <- Spectra()
    res <- removePeaks(sps, threshold = 10)
    expect_true(length(res@processingQueue) == 1)
    expect_equal(res@processingQueue[[1]],
                 ProcessingStep(.peaks_remove,
                                list(threshold = 10, msLevel = integer())))

    sps <- Spectra(sciex_mzr)
    centroided(sps) <- TRUE
    res <- removePeaks(sps, threshold = 5000)
    pks_res <- lapply(sciex_pks, .peaks_remove, threshold = 5000,
                      spectrumMsLevel = 1L, centroided = TRUE)
    expect_identical(peaks(res), SimpleList(pks_res))
})
