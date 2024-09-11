test_that(".valid_processing_queue works", {
    expect_true(is.character(.valid_processing_queue(list(3, 5))))
    lst <- list(ProcessingStep(mean), ProcessingStep("max"))
    expect_true(is.null(.valid_processing_queue(lst)))
})

test_that("addProcessing works", {
    tst <- Spectra()
    tst <- addProcessing(tst, mean)
    expect_true(length(tst@processingQueue) == 1)
    expect_error(addProcessing(tst, "4"))
    tst <- addProcessing(tst, function(z, t) z * t, t = 4)
    expect_true(length(tst@processingQueue) == 2)
    show(tst)
})

test_that("applyProcessing works", {
    ## Initialize required objects.
    sps_mzr <- filterRt(Spectra(sciex_mzr), rt = c(10, 20))
    ## Add processings.
    centroided(sps_mzr) <- TRUE
    sps_mzr <- replaceIntensitiesBelow(sps_mzr, threshold = 5000,
                                       value = NA_real_)
    sps_mzr <- filterIntensity(sps_mzr)
    expect_true(length(sps_mzr@processingQueue) == 2)
    expect_error(applyProcessing(sps_mzr), "is read-only")

    ## Create writeable backends.
    sps_mem <- setBackend(sps_mzr, backend = MsBackendDataFrame())
    sps_h5 <- setBackend(sps_mzr, backend = MsBackendHdf5Peaks(),
                         files = c(tempfile(), tempfile()),
                         f = rep(1, length(sps_mzr)))
    expect_true(length(sps_mem@processingQueue) == 2)
    expect_true(length(sps_h5@processingQueue) == 2)
    expect_identical(peaksData(sps_mzr), peaksData(sps_mem))
    expect_identical(peaksData(sps_h5), peaksData(sps_mem))

    ## MsBackendDataFrame
    res <- applyProcessing(sps_mem)
    expect_true(length(res@processingQueue) == 0)
    expect_true(length(res@processing) > length(sps_mem@processing))
    expect_identical(rtime(res), rtime(sps_mem))
    expect_identical(peaksData(res), peaksData(sps_mem))

    ## MsBackendHdf5Peaks
    res <- applyProcessing(sps_h5)
    expect_true(length(res@processingQueue) == 0)
    expect_true(length(res@processing) > length(sps_h5@processing))
    expect_identical(rtime(res), rtime(sps_mem))
    expect_identical(peaksData(res), peaksData(sps_mem))
    expect_true(all(res@backend@modCount > sps_h5@backend@modCount))

    ## Applying the processing queue invalidated the original object!
    expect_error(peaksData(sps_h5))
    sps_h5 <- setBackend(sps_mzr, backend = MsBackendHdf5Peaks(),
                         files = c(tempfile(), tempfile()),
                         f = rep(1, length(sps_mzr)))

    ## Use an arbitrary splitting factor ensuring that the results are still OK.
    f <- rep(letters[1:9], 8)
    f <- sample(f)

    ## MsBackendHdf5Peaks
    res <- applyProcessing(sps_mem, f = f)
    expect_true(length(res@processingQueue) == 0)
    expect_true(length(res@processing) > length(sps_mem@processing))
    expect_identical(rtime(res), rtime(sps_mem))
    expect_identical(peaksData(res), peaksData(sps_mem))

    ## MsBackendHdf5Peaks: throws an error, because the factor f does not
    ## match the dataStorage.
    expect_error(applyProcessing(sps_h5, f = f))

    sps_h5 <- setBackend(sps_mzr, backend = MsBackendHdf5Peaks(),
                         files = c(tempfile(), tempfile()),
                         f = rep(1, length(sps_mzr)))
    res <- applyProcessing(sps_h5, f = rep(1, length(sps_h5)))
    expect_true(length(res@processingQueue) == 0)
    expect_true(length(res@processing) > length(sps_h5@processing))
    expect_identical(rtime(res), rtime(sps_mem))
    expect_identical(peaksData(res), peaksData(sps_mem))
    expect_true(all(res@backend@modCount > sps_h5@backend@modCount))

    expect_error(applyProcessing(sps_mem, f = 1:2), "has to be equal to the")
})

test_that(".check_ms_level works", {
    expect_true(.check_ms_level(sciex_mzr, 1))
    expect_warning(.check_ms_level(sciex_mzr, 2))
    expect_false(.check_ms_level(sciex_mzr, 2))
    expect_error(.check_ms_level(sciex_mzr, "a"), "must be numeric")

    expect_true(.check_ms_level(tmt_mzr, 1))
    expect_true(.check_ms_level(tmt_mzr, 2))
    expect_true(.check_ms_level(tmt_mzr, c(1, 2)))
    expect_true(.check_ms_level(tmt_mzr, c(1, 4)))
})

test_that(".compare_spectra_self work", {
    sps <- Spectra(sciex_hd5)[120:126]
    sps <- setBackend(sps, MsBackendDataFrame())

    res <- .compare_spectra_chunk(sps, sps)
    expect_true(ncol(res) == length(sps))
    expect_true(nrow(res) == length(sps))
    expect_equal(diag(res), rep(1, length(sps)))

    res_2 <- .compare_spectra_self(sps)
    expect_equal(dim(res), dim(res_2))
    expect_identical(diag(res), diag(res_2))
    expect_identical(res[!lower.tri(res)], res_2[!lower.tri(res_2)])
})

test_that(".compare_spectra_chunk works", {
    sps <- Spectra(sciex_hd5)[120:126]
    sps <- setBackend(sps, MsBackendDataFrame())

    res <- .compare_spectra_chunk(sps, sps)
    expect_true(ncol(res) == length(sps))
    expect_true(nrow(res) == length(sps))
    expect_equal(diag(res), rep(1, length(sps)))

    res_2 <- .compare_spectra_chunk(sps, sps[3])
    expect_true(ncol(res_2) == 1)
    expect_true(nrow(res_2) == length(sps))
    expect_identical(res_2[, 1], res[, 3])

    res_2 <- .compare_spectra_chunk(sps[5], sps)
    expect_true(ncol(res_2) == length(sps))
    expect_true(nrow(res_2) == 1)
    expect_identical(res_2[1, ], res[5, ])

    res_3 <- .compare_spectra_chunk(sps[5], sps, chunkSize = 2)
    expect_equal(res_2, res_3)

    cor_fun <- function(x, y, ...) {
        cor(x[, 2], y[, 2], use = "pairwise.complete.obs")
    }
    res <- .compare_spectra_chunk(sps[1], sps[1], FUN = cor_fun)
    expect_true(res[1, 1] == 1)
    res <- .compare_spectra_chunk(sps[1], sps[2], FUN = cor_fun)
    res_2 <- .compare_spectra_chunk(sps[1], sps[2])
    expect_true(res[1, 1] > res_2[1, 1])
})

test_that(".lapply works", {
    sps <- Spectra(sciex_mzr)[120:126]
    expect_error(.lapply(sps), "missing")
    res <- .lapply(sps, FUN = rtime)
    expect_identical(unlist(res, use.names = FALSE), rtime(sps))

    ## Effect of unsplit: get everything in right order.
    res <- .lapply(sps, FUN = rtime, f = c(4, 1, 6, 7, 2, 3, 5))
    expect_identical(unsplit(res, f = c(4, 1, 6, 7, 2, 3, 5)), rtime(sps))

    ## arbitrary function
    my_fun <- function(x, add) {
        x$rtime + add
    }
    res <- .lapply(sps, FUN = my_fun, add = 3)
    expect_identical(rtime(sps) + 3, unlist(res, use.names = FALSE))

    ## After clean and stuff
    spsc <- filterIntensity(replaceIntensitiesBelow(
        sps, threshold = 4000, value = NA_real_))
    res <- .lapply(spsc, FUN = function(z) sum(intensity(z)))
    res_2 <- vapply(intensity(spsc), sum, numeric(1))
    expect_identical(unlist(res, use.names = FALSE), res_2)
})

test_that(".concatenate_spectra works", {
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

    res <- .concatenate_spectra(list(sp1, sp2, sp3))
    expect_true(is(res, "Spectra"))
    expect_equal(length(res), sum(nrow(df1), nrow(df2), nrow(df3)))
    expect_identical(msLevel(res), c(1L, 1L, 1L, 2L, 2L, 3L, 3L))
    expect_identical(res$other_col, c(NA, NA, NA, NA, NA, "a", "a"))
    expect_true(length(res@processingQueue) == 0)
    expect_true(length(res@processing) == 1)

    ## One Spectra without m/z and intensity
    res <- .concatenate_spectra(list(sp3, sp4))
    expect_true(is(res, "Spectra"))
    expect_identical(mz(res), NumericList(c(1.4, 1.5, 1.6), c(1.8, 1.9),
                                          numeric(), numeric(),
                                          compress = FALSE))
    expect_identical(msLevel(res), rep(3L, 4))
    expect_identical(intensity(res), NumericList(c(123.4, 12, 5), c(43.1, 5),
                                                 numeric(), numeric(),
                                                 compress = FALSE))
    res <- .concatenate_spectra(list(sp4, sp3))
    expect_true(is(res, "Spectra"))
    expect_identical(mz(res), NumericList(numeric(), numeric(),
                                          c(1.4, 1.5, 1.6), c(1.8, 1.9),
                                          compress = FALSE))
    expect_identical(msLevel(res), rep(3L, 4))
    expect_identical(intensity(res), NumericList(numeric(), numeric(),
                                                 c(123.4, 12, 5), c(43.1, 5),
                                                 compress = FALSE))

    ## Two Spectra without m/z and intensity
    res <- .concatenate_spectra(list(sp4, sp4))
    expect_true(is(res, "Spectra"))
    expect_identical(mz(res), NumericList(numeric(), numeric(), numeric(),
                                          numeric(), compress = FALSE))

    sp1@metadata <- list(version = "1.0.0", date = date())
    res <- c(sp1, sp2)
    expect_equal(res@metadata, sp1@metadata)

    sp1@processingQueue <- list(ProcessingStep(sum))
    expect_error(c(sp1, sp2), "with non-empty processing")

    ## Empty Spectra
    sp1 <- Spectra()
    sp2 <- Spectra()
    res <- c(sp1, sp2)

    ## Different backends
    s1 <- Spectra(sciex_mzr)
    s2 <- Spectra(sciex_hd5)
    expect_error(c(s1, s2), "backends of the same type")

    ## BackendMzR
    res <- .concatenate_spectra(list(Spectra(tmt_mzr), Spectra(sciex_mzr)))
    expect_identical(msLevel(res), c(msLevel(tmt_mzr), msLevel(sciex_mzr)))
    expect_identical(msLevel(sciex_mzr), msLevel(res[dataStorage(res) %in%
                                                     sciex_file]))
    expect_identical(msLevel(tmt_mzr), msLevel(res[dataStorage(res) ==
                                                   dataStorage(tmt_mzr)[1]]))
})

test_that(".combine_spectra works", {
    spd <- DataFrame(msLevel = c(2L, 2L, 2L), rtime = c(1, 2, 3))
    spd$mz <- list(c(12, 14, 45, 56), c(14.1, 34, 56.1), c(12.1, 14.15, 34.1))
    spd$intensity <- list(c(10, 20, 30, 40), c(11, 21, 31), c(12, 22, 32))
    sps <- Spectra(spd)
    res <- .combine_spectra(sps)
    expect_true(length(res) == 1)
    expect_equal(res$mz[[1]], sort(unlist(spd$mz)))

    res <- .combine_spectra(sps, FUN = combinePeaksData, tolerance = 0.1)
    expect_true(length(res) == 1)
    expect_equal(res$mz[[1]], c(mean(c(12, 12.1)), mean(c(14, 14.1, 14.15)),
                                mean(c(34, 34.1)), 45, mean(c(56, 56.1))))
    expect_equal(res$intensity[[1]], c(mean(c(10, 12)), mean(c(20, 11, 22)),
                                       mean(c(21, 32)), 30, mean(c(40, 31))))
    res <- .combine_spectra(sps, FUN = combinePeaksData, tolerance = 0.1,
                            mzFun = max, intensityFun = median)
    expect_true(length(res) == 1)
    expect_equal(res$mz[[1]], c(max(c(12, 12.1)), max(c(14, 14.1, 14.15)),
                                max(c(34, 34.1)), 45, max(c(56, 56.1))))
    expect_equal(res$intensity[[1]],
                 c(median(c(10, 12)), median(c(20, 11, 22)),
                   median(c(21, 32)), 30, median(c(40, 31))))

    ## See if it works with MsBackendMzR
    sps <- Spectra(sciex_mzr)
    res <- .combine_spectra(sps, tolerance = 0.1, FUN = combinePeaksData)
    expect_true(length(res) == 2)
    expect_true(is(res, "Spectra"))
    expect_true(class(res@backend) == "MsBackendMemory")
    expect_true(length(unlist(res$mz)) < length(unlist(sps$mz)))
})

test_that("combineSpectra works", {
    spd <- DataFrame(msLevel = c(2L, 2L, 2L), rtime = c(1, 2, 3))
    spd$mz <- list(c(12, 14, 45, 56), c(14.1, 34, 56.1), c(12.1, 14.15, 34.1))
    spd$intensity <- list(c(10, 20, 30, 40), c(11, 21, 31), c(12, 22, 32))
    sps <- Spectra(spd)
    res <- combineSpectra(sps, tolerance = 0.1)
    expect_true(is(res, "Spectra"))
    expect_true(length(res) == 1)

    sps <- rev(Spectra(sciex_mzr))
    res <- combineSpectra(sps, tolerance = 0.1)

    expect_true(is(res, "Spectra"))
    expect_true(length(res) == 2)
    expect_true(class(res@backend) == "MsBackendMemory")
    expect_equal(res$dataOrigin, unique(sps$dataStorage))

    ## Different f
    sps$crude_rtime <- as.integer(rtime(sps))
    res <- combineSpectra(sps, tolerance = 0.1, f = sps$crude_rtime)
    expect_equal(unique(res$dataOrigin), unique(sps$dataStorage))
    fls <- unique(res$dataOrigin)
    expect_equal(res$crude_rtime[res$dataOrigin == fls[1]],
                 unique(sps$crude_rtime[sps$dataOrigin == fls[1]]))
    expect_equal(res$crude_rtime[res$dataOrigin == fls[2]],
                 unique(sps$crude_rtime[sps$dataOrigin == fls[2]]))

    ## With a non-empty processing queue.
    spd <- DataFrame(msLevel = c(2L, 2L, 2L), rtime = c(1, 2, 3))
    spd$mz <- list(c(12, 14, 45, 56), c(14.1, 34, 56.1), c(12.1, 14.15, 34.1))
    spd$intensity <- list(c(10, 20, 30, 40), c(11, 21, 31), c(12, 22, 32))
    sps <- Spectra(spd)
    sps <- filterIntensity(sps, 30)
    res <- combineSpectra(sps)
    expect_true(length(res) == 1)
    expect_true(all(intensity(res)[[1]] >= 30))
    expect_true(length(sps@processingQueue) > 0)

    ## Define p do check if parallel processing does the same.
    sps2 <- c(Spectra(spd), Spectra(spd))
    sps2 <- filterIntensity(sps2, 30)
    res2 <- combineSpectra(sps2, f = rep(1:2, each = 3), p = rep(1:2, each = 3))
    expect_equal(mz(res)[[1]], mz(res2)[[1]])
    expect_true(all(intensity(res2)[[1]] >= 30))

    ## With a non-empty processing queue and an hdf5 backend.
    h5p <- tempdir()
    sps2$scanIndex <- seq_along(sps2)
    sps2 <- setBackend(sps2, MsBackendHdf5Peaks(), hdf5path = h5p)
    expect_true(length(sps2@processingQueue) == 1)
    res3 <- combineSpectra(sps2, f = rep(1:2, each = 3), p = rep(1, 6))
    expect_equal(unname(intensity(res2)), intensity(res3))
    expect_equal(unname(mz(res2)), mz(res3))
    expect_true(length(res3@processingQueue) == 0)
    expect_true(length(sps2@processingQueue) == 1)
    expect_true(validObject(sps2))
    expect_error(mz(sps2), "have changed")
})

test_that("dropNaSpectraVariables works", {
    ## with a MsBackend
    res <- dropNaSpectraVariables(sciex_mzr)
    expect_true(all(vapply1l(res@spectraData, function(z) !any(is.na(z)))))
    ## with a Spectra
    sps <- Spectra(sciex_mzr)
    res <- dropNaSpectraVariables(sps)
    expect_true(all(vapply1l(res@backend@spectraData,
                             function(z) !any(is.na(z)))))
})

test_that(".has_mz works", {
    sps <- Spectra(sciex_mzr)[1:10]
    sps <- setBackend(sps, MsBackendDataFrame())
    mzs <- mz(sps)
    x <- c(mzs[[2]][5], mzs[[3]][8])

    res <- .has_mz(sps, mz = x, ppm = 0)
    expect_true(length(res) == length(sps))
    expect_true(is.logical(res))

    spd <- DataFrame(msLevel = c(2L, 2L, 2L), rtime = c(1, 2, 3))
    spd$mz <- list(c(12, 14, 45, 56), c(14.1, 34, 56.1), c(12.1, 14.15, 34.1))
    spd$intensity <- list(c(10, 20, 30, 40), c(11, 21, 31), c(12, 22, 32))
    sps <- Spectra(spd)

    res <- .has_mz(sps, mz = c(14, 34))
    expect_equal(res, c(TRUE, TRUE, FALSE))
    res <- .has_mz(sps, mz = c(14, 34), tolerance = 0.15)
    expect_equal(res, c(TRUE, TRUE, TRUE))

    res <- .has_mz(sps, mz = c(14, 34), condFun = all)
    expect_true(all(!res))
    res <- .has_mz(sps, mz = c(14, 34), condFun = all, tolerance = 0.15)
    expect_equal(res, c(FALSE, TRUE, TRUE))
})

test_that(".has_mz_each works", {
    spd <- DataFrame(msLevel = c(2L, 2L, 2L), rtime = c(1, 2, 3))
    spd$mz <- list(c(12, 14, 45, 56), c(14.1, 34, 56.1), c(12.1, 14.15, 34.1))
    spd$intensity <- list(c(10, 20, 30, 40), c(11, 21, 31), c(12, 22, 32))
    sps <- Spectra(spd)

    res <- .has_mz_each(sps, mz = c(14, 34, 12.1), ppm = 0)
    expect_true(is.logical(res))
    expect_true(length(res) == length(sps))
    expect_equal(res, c(TRUE, TRUE, TRUE))

    res <- .has_mz_each(sps, mz = c(NA, 34, 34))
    expect_equal(res, c(NA, TRUE, FALSE))

    res <- .has_mz_each(sps, mz = c(14, 14, 14), tolerance = 0.1)
    expect_equal(res, c(TRUE, TRUE, FALSE))
})


test_that("joinSpectraData works", {
    ms <- Spectra(msdata::proteomics(pattern = "TMT10", full.names = TRUE))
    k <- sample(length(ms), 10)
    spd2 <- spd1 <- DataFrame(id = ms$spectrumId[k],
                              X = 1:10,
                              Y = letters[1:10])
    ## Input errors
    expect_error(joinSpectraData(ms, spd2, by.x = 1:2),
                 "'by.x' and 'by.y' must be of length 1.")
    expect_error(joinSpectraData(ms, spd2, by.y = 1:2),
                 "'by.x' and 'by.y' must be of length 1.")
    expect_error(joinSpectraData(ms, spd2, by.x = 1),
                 "'by.x' and 'by.y' must be characters.")
    expect_error(joinSpectraData(ms, spd2, by.y = 1),
                 "'by.x' and 'by.y' must be characters.")
    expect_error(joinSpectraData(ms, spd2, by.x = "id"),
                 "'by.x' not found in spectra variables.")
    ## Works with by.y provided
    ms1 <- joinSpectraData(ms, spd1, by.y = "id")
    ## Error with wrong by.y
    expect_error(joinSpectraData(ms, spd2),
                 "'by.y' not found.")
    names(spd2)[1] <- "spectrumId"
    ## Works with default by.y
    ms2 <- joinSpectraData(ms, spd2)
    ## Expected results
    expect_identical(spectraData(ms1[k]),
                     spectraData(ms2[k]))
    expect_identical(c(spectraVariables(ms), "X", "Y"),
                     spectraVariables(ms1))
    expect_identical(spectraData(ms1)$X[k], spd1$X)
    expect_identical(spectraData(ms1)$Y[k], spd1$Y)
    ## Test suffix.y
    spd2$msLevel <- 2
    ms3 <- joinSpectraData(ms, spd2)
    expect_identical(c(spectraVariables(ms), "X", "Y", "msLevel.y"),
                     spectraVariables(ms3))
    ## Check classes after merging
    spd2$NumList <- List(mapply(rep, 1.0, 1:10))
    spd2$IntList <- List(mapply(rep, 1L, 1:10))
    spd2$CharList <- List(mapply(rep, "a", 1:10))
    ms2 <- joinSpectraData(ms, spd2)
    expect_true(is(ms2$X, "integer"))
    expect_true(is(ms2$Y, "character"))
    expect_true(is(ms2$NumList, "NumericList"))
    expect_true(is(ms2$IntList, "IntegerList"))
    expect_true(is(ms2$CharList, "CharacterList"))

    ## With MsBackendMemory
    ms <- Spectra(msdata::proteomics(pattern = "TMT10", full.names = TRUE))
    ms <- setBackend(ms, MsBackendMemory(), BPPARAM = SerialParam())
    k <- sample(length(ms), 10)
    spd2 <- spd1 <- DataFrame(id = ms$spectrumId[k],
                              X = 1:10,
                              Y = letters[1:10])
    ## Input errors
    expect_error(joinSpectraData(ms, spd2, by.x = 1:2),
                 "'by.x' and 'by.y' must be of length 1.")
    expect_error(joinSpectraData(ms, spd2, by.y = 1:2),
                 "'by.x' and 'by.y' must be of length 1.")
    expect_error(joinSpectraData(ms, spd2, by.x = 1),
                 "'by.x' and 'by.y' must be characters.")
    expect_error(joinSpectraData(ms, spd2, by.y = 1),
                 "'by.x' and 'by.y' must be characters.")
    expect_error(joinSpectraData(ms, spd2, by.x = "id"),
                 "'by.x' not found in spectra variables.")
    ## Works with by.y provided
    ms1 <- joinSpectraData(ms, spd1, by.y = "id")
    ## Error with wrong by.y
    expect_error(joinSpectraData(ms, spd2),
                 "'by.y' not found.")
    names(spd2)[1] <- "spectrumId"
    ## Works with default by.y
    ms2 <- joinSpectraData(ms, spd2)
    ## Expected results
    expect_identical(spectraData(ms1[k]),
                     spectraData(ms2[k]))
    expect_identical(c(spectraVariables(ms), "X", "Y"),
                     spectraVariables(ms1))
    expect_identical(spectraData(ms1)$X[k], spd1$X)
    expect_identical(spectraData(ms1)$Y[k], spd1$Y)
    ## Test suffix.y
    spd2$msLevel <- 2
    ms3 <- joinSpectraData(ms, spd2)
    expect_identical(c(spectraVariables(ms), "X", "Y", "msLevel.y"),
                     spectraVariables(ms3))
    ## Use S3 instead of S4 objects here.
    spd2$num_list <- mapply(rep, 1.0, 1:10)
    spd2$int_list <- mapply(rep, 1L, 1:10)
    spd2$char_list <- mapply(rep, "a", 1:10)

    ms2 <- joinSpectraData(ms, spd2)
    expect_true(is(ms2$X, "integer"))
    expect_true(is(ms2$Y, "character"))
    expect_true(is.list(ms2$num_list))
    expect_true(is.list(ms2$int_list))
    expect_true(is.list(ms2$char_list))

    ## Same with S4 objects
    spd2$NumList <- List(mapply(rep, 1.0, 1:10))
    spd2$IntList <- List(mapply(rep, 1L, 1:10))
    spd2$CharList <- List(mapply(rep, "a", 1:10))
    ms2 <- joinSpectraData(ms, spd2)
    expect_true(is(ms2$X, "integer"))
    expect_true(is(ms2$Y, "character"))
    expect_true(is(ms2$NumList, "NumericList"))
    expect_true(is(ms2$IntList, "IntegerList"))
    expect_true(is(ms2$CharList, "CharacterList"))
})

test_that("joinSpectraData key checks work", {
    spd <- DataFrame(msLevel = rep(2L, 3),
                     rtime = c(1.1, 1.2, 1.3),
                     key = paste0("sp", c(1, 1, 3)))
    spd$mz <- list(c(100, 103.2, 104.3, 106.5),
                   c(45.6, 120.4, 190.2),
                   c(45.6, 120.4, 190.2))
    spd$intensity <- list(c(200, 400, 34.2, 17),
                          c(12.3, 15.2, 6.8),
                          c(12.3, 15.2, 6.8))
    sp <- Spectra(spd)
    df <- DataFrame(key = paste0("sp", c(1, 1, 3)),
                    var1 = c(10, 20, 30))

    ## Duplicates in `x` key aren't allowed
    expect_error(joinSpectraData(sp, df, by.x = "key"))

    ## Duplicates in `y` key throw a warning
    sp$key <- paste0("sp", 1:3)
    expect_warning(res <- joinSpectraData(sp, df, by.x = "key"))
    expect_identical(res$var1, c(20, NA, 30))

    ## No duplicates in any key
    df$key <- paste0("sp", 1:3)
    res <- joinSpectraData(sp, df, by.x = "key")
    expect_identical(res$var1, c(10, 20, 30))
})


test_that("processingLog works", {
    sps <- Spectra()
    expect_equal(processingLog(sps), character())
    sps@processing <- c("a", "b", "c")
    expect_equal(processingLog(sps), c("a", "b", "c"))
})

test_that(".processingQueueVariables works", {
    sps <- Spectra()
    expect_equal(.processingQueueVariables(sps), character())
})

test_that(".peaksapply works", {
    ## Use a processing with precursorMz and msLevel
    fl <- dir(system.file("TripleTOF-SWATH", package = "msdata"),
              full.names = TRUE)[1]
    sps <- Spectra(backendInitialize(MsBackendMzR(), files = fl))

    loss <- function(x, spectrumMsLevel, precursorMz, ...) {
        if (spectrumMsLevel == 2L)
            x[, "mz"] <- x[, "mz"] - precursorMz
        x
    }
    sps_2 <- addProcessing(sps, loss,
                           spectraVariables = c("msLevel", "precursorMz"))
    res_2 <- .peaksapply(sps_2)
    mzs_2 <- IRanges::NumericList(lapply(res_2, function(z) z[, "mz"]),
                                  compress = FALSE)
    mzs <- mz(sps)
    ms2 <- sps$msLevel == 2L
    mzs[ms2] <- mzs[ms2] - precursorMz(sps[ms2])
    expect_equal(mzs, mzs_2)

    sps <- Spectra(backend = sciex_mzr)
    sps$centroided <- TRUE
    res <- .peaksapply(sps, FUN = .peaks_replace_intensity, t = 50000,
                       spectraVariables = c("msLevel", "centroided"))
    expect_true(is.list(res))
    expect_equal(length(res), length(sps))
    expect_true(all(vapply(res, is.matrix, logical(1))))

    ## Ensure that this works with arbitrary ordering of the factor f
    res2 <- .peaksapply(sps, FUN = .peaks_replace_intensity, t = 50000,
                        f = rep(1:2, length(sps)/2),
                        spectraVariables = c("msLevel", "centroided"))
    expect_identical(res, res2)

    sps@processingQueue <- list(
        ProcessingStep(.peaks_replace_intensity, list(t = 50000)))
    res_2 <- .peaksapply(sps, spectraVariables = c("msLevel", "centroided"))
    expect_equal(res, res_2)

    res_3 <- .peaksapply(sps, FUN = .peaks_filter_intensity,
                         intensity = c(0.1, Inf),
                         spectraVariables = c("msLevel", "centroided"))
    expect_true(all(vapply(res_3, nrow, integer(1)) <
                    vapply(res_2, nrow, integer(1))))
    expect_true(!any(vapply(res_3, function(z) any(z[, 2] == 0), logical(1))))

    sps@processingQueue <- c(sps@processingQueue,
                             list(ProcessingStep(.peaks_filter_intensity,
                                                 list(intensity = c(0.1, Inf)))))
    res_4 <- .peaksapply(sps, spectraVariables = c("msLevel", "centroided"))
    expect_equal(res_3, res_4)

    ## processing queue and f.
    a <- sps[1:10]
    ref <- .peaksapply(a, spectraVariables = c("msLevel", "centroided"),
                       f = rep(1, 10))
    res <- .peaksapply(a, spectraVariables = c("msLevel", "centroided"),
                       f = c(3, 3, 3, 3, 2, 2, 2, 2, 1, 1))
    expect_equal(ref, res)
    res <- .peaksapply(a, spectraVariables = c("msLevel", "centroided"),
                       f = NULL)
    expect_equal(ref, res)

    ## no processing queue and f.
    a@processingQueue <- list()
    ref <- .peaksapply(a, spectraVariables = c("msLevel", "centroided"),
                       f = rep(1, 10))
    res <- .peaksapply(a, spectraVariables = c("msLevel", "centroided"),
                       f = c(3, 3, 3, 3, 2, 2, 2, 2, 1, 1))
    expect_equal(ref, res)
    res <- .peaksapply(a, spectraVariables = c("msLevel", "centroided"),
                       f = NULL)
    expect_equal(ref, res)
})

test_that(".apply_processing_queue works", {
    inp <- list(1:5, 1:3, 5)
    expect_equal(.apply_processing_queue(inp), inp)
    res <- .apply_processing_queue(inp, queue = list(ProcessingStep("sum")))
    expect_equal(res, list(sum(1:5), sum(1:3), 5))

    q <- list(ProcessingStep(function(x, y, ...) x + y, ARGS = list(y = 3)),
              ProcessingStep(function(x, y, ...) x - y, ARGS = list(y = 1)))
    res <- .apply_processing_queue(inp, queue = q)
    expect_equal(res, list((1:5 + 2), (1:3 + 2), 7))

    be <- sciex_mzr
    pks <- peaksData(be)
    pq <- list(ProcessingStep(.peaks_replace_intensity, list(t = 50000)))
    spd <- spectraData(be, columns = c("msLevel", "centroided"))
    spd$centroided <- TRUE
    res <- .apply_processing_queue(
        pks, spectraData = as.data.frame(spd), queue = pq)
    expect_true(all(vapply(res, function(z) all(z[z[, 2] > 0, 2] > 50000),
                           logical(1))))
    expect_equal(vapply(res, nrow, integer(1)), vapply(pks, nrow, integer(1)))

    ## Length 2
    pq <- c(pq, list(ProcessingStep(.peaks_filter_intensity,
                                    list(intensity = c(0.1, Inf)))))
    res <- .apply_processing_queue(pks, spd, queue = pq)
    expect_true(all(vapply(res, function(z) all(z[z[, 2] > 0, 2] > 50000),
                           logical(1))))
    expect_true(all(vapply(res, nrow, integer(1)) <
                    vapply(pks, nrow, integer(1))))
})

test_that(".estimate_precursor_intensity works", {
    fl <- msdata::proteomics("MS3TMT11.mzML", full.names = TRUE)
    tmp <- Spectra(fl, backend = MsBackendMzR())

    ## previous
    res <- .estimate_precursor_intensity(tmp)
    expect_true(all(is.na(res[msLevel(tmp) == 1L])))
    expect_true(all(is.na(res[msLevel(tmp) == 3L])))
    expect_warning(.estimate_precursor_intensity(tmp, msLevel = 3L),
                   "not yet validated")
    expect_true(
        cor(res, precursorIntensity(tmp), use = "pairwise.complete.obs") > 0.9)

    ## arbitrary order
    idx <- sample(seq_along(tmp))
    res_2 <- .estimate_precursor_intensity(tmp[idx])
    expect_equal(res_2, res[idx])

    ## interpolation
    res_2 <- .estimate_precursor_intensity(tmp, method = "interpolation")
    expect_true(is.character(all.equal(res, res_2)))
    expect_true(cor(res, res_2, use = "pairwise.complete.obs") > 0.99)

    ## no MS1
    tmp_2 <- filterMsLevel(tmp, msLevel = 2:3)
    res <- .estimate_precursor_intensity(tmp_2)
    expect_true(all(is.na(res)))

    ## no MS2
    tmp_2 <- filterMsLevel(tmp, msLevel = 1L)
    res <- .estimate_precursor_intensity(tmp_2)
    expect_true(all(is.na(res)))
})

test_that(".chunk_factor works", {
    res <- .chunk_factor(10, chunkSize = 3)
    expect_equal(res, as.factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4)))
    res <- .chunk_factor(10)
    expect_equal(res, as.factor(rep(1L, 10)))
})

test_that("chunkapply works", {
    smem <- setBackend(Spectra(sciex_mzr), MsBackendMemory())
    res <- chunkapply(smem, lengths, chunkSize = 10)
    expect_equal(res, lengths(smem))
    a <- smem[1:10]
    chnks <- as.factor(c(1, 1, 1, 2, 2, 2, 3, 3))
    expect_error(chunkapply(a, lengths, chunks = chnks), "does not match")
    chnks <- as.factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4))
    res <- chunkapply(a, lengths, chunks = chnks)
    expect_equal(res, lengths(a))

    chnks <- as.factor(c(2, 2, 2, 1, 1, 1, 4, 3, 3, 3))
    res2 <- chunkapply(a, lengths, chunks = chnks)
    expect_equal(res2, res)
})

test_that("deisotopeSpectra works", {
    res <- deisotopeSpectra(sps_dia[1:10], ppm = 20)
    expect_true(length(res@processingQueue) > length(sps_dia@processingQueue))
    expect_true(all(lengths(res) < lengths(sps_dia[1:10])))
})

test_that("reduceSpectra works", {
    res <- reduceSpectra(sps_dia[14], ppm = 20, tolerance = 0.1)
    expect_true(length(res@processingQueue) > length(sps_dia@processingQueue))
    expect_true(nrow(peaksData(res)[[1L]]) < nrow(peaksData(sps_dia[14])[[1L]]))
})

test_that("filterPrecursorMaxIntensity works", {
    ms1 <- filterMsLevel(sps_dda, 1L)
    res <- filterPrecursorMaxIntensity(ms1)
    expect_equal(rtime(res), rtime(ms1))
    expect_true(length(res@processing) > length(ms1@processing))

    expect_equal(rtime(filterMsLevel(res, 1L)),
                 rtime(filterMsLevel(sps_dda, 1L)))

    res <- filterPrecursorMaxIntensity(filterRt(sps_dda, c(200, 300)))

    ## calculate manually:
    tmp <- filterRt(filterMsLevel(sps_dda, 2L), c(200, 300))
    idx <- order(precursorMz(tmp))
    tmp <- tmp[idx]

    grps <- MsCoreUtils::group(precursorMz(tmp), tolerance = 0, ppm = 20)
    tmpl <- split(tmp, as.factor(grps))
    ref <- lapply(tmpl, function(z) {
        z[which.max(precursorIntensity(z))]
    })
    ref <- concatenateSpectra(ref)      # maybe run on subset only...
    expect_equal(sort(rtime(ref)), rtime(filterMsLevel(res, 2L)))

    ## Artificial data.
    tmp <- sps_dda[1:12]
    tmp$precursorMz <- c(NA, 13, NA, 5, 5, NA, 13, 5, 3, NA, 13, NA)
    tmp$precursorIntensity <- c(NA, 4, NA, 3, 200, NA, 8, 100, 23, NA, 400, NA)
    res <- filterPrecursorMaxIntensity(tmp)
    expect_true(length(res) == 8)
    expect_equal(rtime(res), rtime(tmp[c(1, 3, 5, 6, 9, 10, 11, 12)]))
})

test_that("filterPrecursorIsotopes works", {
    x <- sps_dda[1:10]
    x$precursorMz <- c(NA, 803.4, 804.41, 115, NA, 805.41, 116, 123.2, NA, NA)
    x$precursorIntensity <- c(NA, 100, 42.7, 100, NA, 15.58, 4.7, 100, NA, NA)

    res <- filterPrecursorIsotopes(x, ppm = 40)
    expect_equal(precursorIntensity(res), c(NA, 100, 100, NA, 100, NA, NA))
    expect_equal(rtime(res), rtime(x)[c(1, 2, 4, 5, 8, 9, 10)])

    ## single isotopologue group
    res <- filterPrecursorIsotopes(x, ppm = 20)
    expect_equal(precursorIntensity(res), c(NA, 100, 100, NA, 4.7, 100, NA, NA))
    expect_equal(rtime(res), rtime(x)[c(1, 2, 4, 5, 7, 8, 9, 10)])

    ## DDA data will not have isotope peaks
    x <- filterRt(sps_dda, c(200, 300))
    x$precursorIntensity <- estimatePrecursorIntensity(x)
    res <- filterPrecursorIsotopes(x)
    expect_true(length(res) < length(x))

    x <- filterMsLevel(sps_dda, 1L)
    res <- filterPrecursorIsotopes(x)
    expect_equal(rtime(res), rtime(x))
})

test_that("scalePeaks works", {
    tmp <- data.frame(msLevel = c(1L, 2L, 2L, 3L), rtime = 1:4)
    tmp$mz <- list(1:3, 1:2, 1:4, 1:3)
    tmp$intensity <- list(c(12, 32.2, 12.1), c(34, 35.2),
                          c(1, 2, 3, 4), c(112, 341, 532))
    sps <- Spectra(tmp)
    res <- scalePeaks(sps)
    expect_true(length(res@processingQueue) > length(sps@processingQueue))
    expect_true(length(res@processing) > length(sps@processing))

    expect_true(all(sum(intensity(res)) == 1))
    expect_true(all(unlist(intensity(res) < 1)))

    res <- scalePeaks(sps, by = max)
    expect_true(all(max(intensity(res)) == 1))

    res <- scalePeaks(sps, by = sum, msLevel. = 2)
    expect_equal(res$intensity[[1L]], sps$intensity[[1L]])
    expect_true(sum(intensity(res)[[2L]]) == 1)
    expect_true(sum(intensity(res)[[3L]]) == 1)
    expect_equal(res$intensity[[4L]], sps$intensity[[4L]])
})

test_that("filterPrecursorPeaks,Spectra works", {
    x <- Spectra(tmt_mzr[5:15])
    expect_error(filterPrecursorPeaks(4), "instance of class")
    expect_error(filterPrecursorPeaks(x, tolerance = c(3, 4)), "length 1")
    expect_error(filterPrecursorPeaks(x, ppm = c(12.3, 24.4)), "length 1")
    expect_warning(res <- filterPrecursorPeaks(x, msLevel. = 5), "available")
    expect_equal(res, x)

    res <- filterPrecursorPeaks(x, mz = "==", tolerance = 0.2)
    expect_true(length(res@processing) > 0)
    expect_true(lengths(res)[1L] < lengths(x)[1L])

    res <- filterPrecursorPeaks(x, mz = ">=")
    expect_true(all(lengths(res)[msLevel(x) > 1L] <
                    lengths(x)[msLevel(x) > 1L]))
    expect_equal(lengths(res)[msLevel(x) == 1L], lengths(x)[msLevel(x) == 1L])
})

test_that("processingChunkSize works", {
    s <- Spectra()
    expect_equal(processingChunkSize(s), Inf)
    processingChunkSize(s) <- 1000
    expect_equal(processingChunkSize(s), 1000)
    expect_error(processingChunkSize(s) <- c(1, 2), "length 1")
    expect_error(processingChunkSize(s) <- "A", "character")
})

test_that("processingChunkFactor works", {
    s <- Spectra()
    expect_equal(processingChunkFactor(s), factor())
    tmp <- Spectra(sciex_mzr)

    expect_equal(length(processingChunkFactor(tmp)), length(tmp))
    expect_true(is.factor(processingChunkFactor(tmp)))

    processingChunkSize(tmp) <- 1000
    res <- processingChunkFactor(tmp)
    expect_true(is.factor(res))
    expect_true(length(res) == length(tmp))
    expect_equal(levels(res), c("1", "2"))

    expect_equal(.parallel_processing_factor(tmp), processingChunkFactor(tmp))

    expect_error(processingChunkFactor("a"), "Spectra")
})

test_that("filterPeaksRanges,Spectra works", {
    df <- data.frame(rtime = 123.3, new_var = 4, msLevel = 2L)
    df$mz <- list(c(100.1, 100.2, 100.3, 100.4, 200.1, 200.2, 200.3,
                    300.1, 300.3, 300.4, 300.5))
    df$intensity <- list(1:11)
    s <- Spectra(df)
    ## Check errors
    expect_error(filterPeaksRanges(3), "'Spectra' object")
    expect_error(filterPeaksRanges(s, rtime = c(1, 2), not_exist = c(1, 2)),
                 "valid spectra variables")
    expect_error(filterPeaksRanges(s, rtime = 2, mz = c(1, 2)),
                 "'numeric' of length 2")
    expect_error(filterPeaksRanges(
        s, rtime = rbind(c(1, 2), c(2, 3)), mz = c(1, 2)),
        "Number of rows of the range matrices")

    ## Single range per variable
    res <- filterPeaksRanges(s, rtime = c(100, 200), mz = cbind(200, 300))
    expect_true(inherits(res, "Spectra"))
    expect_true(length(res@processingQueue) > 0L)
    expect_equal(res@processingQueueVariables, c("rtime", "msLevel"))
    expect_equal(length(res@processing), 1L)
    a <- peaksData(res)[[1L]]
    expect_equal(a[, 2L], c(5:7))
    res <- filterPeaksRanges(s, rtime = c(100, 200), mz = cbind(200, 300),
                             keep = FALSE)
    a <- peaksData(res)[[1L]]
    expect_equal(a[, 2L], c(1:4, 8:11))

    ## Multiple ranges per variable
    res <- filterPeaksRanges(
        s, new_var = rbind(c(1, 8), c(1, 4), c(1, 5)),
        rtime = rbind(c(100, 200), c(400, 500), c(100, 200)),
        mz = rbind(c(100, 100.3), c(0, 500), c(300.3, 310)))
    expect_true(inherits(res, "Spectra"))
    a <- peaksData(res)[[1L]]
    expect_equal(a[, 2L], c(1:3, 9:11))
    res <- filterPeaksRanges(
        s, new_var = rbind(c(1, 8), c(1, 4), c(1, 5)),
        rtime = rbind(c(100, 200), c(400, 500), c(100, 200)),
        mz = rbind(c(100, 100.3), c(0, 500), c(300.3, 310)), keep = FALSE)
    expect_true(inherits(res, "Spectra"))
    a <- peaksData(res)[[1L]]
    expect_equal(a[, 2L], c(4:8))

    ## Filter also with msLevel; to have the same behaviour as with other
    ## filters we would need to add a second filter for e.g. MS level 2
    s <- c(s, s)
    s$msLevel <- c(1L, 2L)
    res <- filterPeaksRanges(s, rtime = c(100, 200), msLevel = c(1, 1),
                             mz = c(100, 200))
    a <- peaksData(res)[[1L]]
    expect_equal(a[, 2L], 1:4)
    a <- peaksData(res)[[2L]]
    expect_true(nrow(a) == 0L)
    res <- filterPeaksRanges(s, rtime = rbind(c(100, 200), c(100, 200)),
                             msLevel = rbind(c(1, 1), c(2, 2)),
                             mz = rbind(c(100, 200), c(0, 400)))
    a <- peaksData(res)[[1L]]
    expect_equal(a[, 2L], 1:4)
    a <- peaksData(res)[[2L]]
    expect_equal(a[, 2L], 1:11)
    res <- filterPeaksRanges(s, rtime = rbind(c(100, 200), c(100, 200)),
                             msLevel = rbind(c(1, 1), c(2, 2)),
                             mz = rbind(c(100, 200), c(0, 400)),
                             keep = FALSE)
    a <- peaksData(res)[[1L]]
    expect_equal(a[, 2L], 5:11)
    a <- peaksData(res)[[2L]]
    expect_true(nrow(a) == 0)
})

test_that(".spectra_to_spectrum_list works", {
    ## With MSnbase in Enhances we would need to call it like below.
    if (requireNamespace("MSnbase", quietly = TRUE)) {
        ## library(MSnbase) # MSnbase needs to be in Suggests
        a <- .spectra_to_spectrum_list(sps_dda, chunkSize = 5000)
        expect_true(is.list(a))
        expect_equal(length(a), length(sps_dda))
    }
})
