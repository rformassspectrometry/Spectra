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
})

test_that(".apply_processing_queue works", {
    inp <- list(1:5, 1:3, 5)
    expect_equal(.apply_processing_queue(inp), inp)
    res <- .apply_processing_queue(inp, msLevel = rep(0, 3),
                                   centroided = rep(FALSE, 3),
                                   list(ProcessingStep("sum")))
    expect_equal(res, list(sum(1:5), sum(1:3), 5))

    q <- list(ProcessingStep(function(x, y, ...) x + y, ARGS = list(y = 3)),
              ProcessingStep(function(x, y, ...) x - y, ARGS = list(y = 1)))
    res <- .apply_processing_queue(inp, msLevel = rep(0, 3),
                                   centroided = rep(FALSE, 3), q)
    expect_equal(res, list((1:5 + 2), (1:3 + 2), 7))

    be <- sciex_mzr
    pks <- as.list(be)
    pq <- list(ProcessingStep(.peaks_replace_intensity, list(t = 50000)))
    res <- .apply_processing_queue(pks, msLevel(be),
                                   rep(TRUE, length(be)), pq)
    expect_true(all(vapply(res, function(z) all(z[z[, 2] > 0, 2] > 50000),
                           logical(1))))
    expect_equal(vapply(res, nrow, integer(1)), vapply(pks, nrow, integer(1)))

    ## Length 2
    pq <- c(pq, list(ProcessingStep(.peaks_filter_intensity,
                                    list(intensity = c(0.1, Inf)))))
    res <- .apply_processing_queue(pks, msLevel(be),
                                   rep(TRUE, length(be)), pq)
    expect_true(all(vapply(res, function(z) all(z[z[, 2] > 0, 2] > 50000),
                           logical(1))))
    expect_true(all(vapply(res, nrow, integer(1)) <
                    vapply(pks, nrow, integer(1))))
})

test_that(".peaksapply works", {
    sps <- Spectra(backend = sciex_mzr)
    res <- .peaksapply(sps, FUN = .peaks_replace_intensity, t = 50000)
    expect_true(is.list(res))
    expect_equal(length(res), length(sps))
    expect_true(all(vapply(res, is.matrix, logical(1))))

    ## Ensure that this works with arbitrary ordering of the factor f
    res2 <- .peaksapply(sps, FUN = .peaks_replace_intensity, t = 50000,
                        f = rep(1:2, length(sps)/2))
    expect_identical(res, res2)

    sps@processingQueue <- list(
        ProcessingStep(.peaks_replace_intensity, list(t = 50000)))
    res_2 <- .peaksapply(sps)
    expect_equal(res, res_2)

    res_3 <- .peaksapply(sps, FUN = .peaks_filter_intensity,
                         intensity = c(0.1, Inf))
    expect_true(all(vapply(res_3, nrow, integer(1)) <
                    vapply(res_2, nrow, integer(1))))
    expect_true(!any(vapply(res_3, function(z) any(z[, 2] == 0), logical(1))))

    sps@processingQueue <- c(sps@processingQueue,
                             list(ProcessingStep(.peaks_filter_intensity,
                                                 list(intensity = c(0.1, Inf)))))
    res_4 <- .peaksapply(sps)
    expect_equal(res_3, res_4)
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
    expect_identical(as.list(sps_mzr), as.list(sps_mem))
    expect_identical(as.list(sps_h5), as.list(sps_mem))

    ## MsBackendDataFrame
    res <- applyProcessing(sps_mem)
    expect_true(length(res@processingQueue) == 0)
    expect_true(length(res@processing) > length(sps_mem@processing))
    expect_identical(rtime(res), rtime(sps_mem))
    expect_identical(as.list(res), as.list(sps_mem))

    ## MsBackendHdf5Peaks
    res <- applyProcessing(sps_h5)
    expect_true(length(res@processingQueue) == 0)
    expect_true(length(res@processing) > length(sps_h5@processing))
    expect_identical(rtime(res), rtime(sps_mem))
    expect_identical(as.list(res), as.list(sps_mem))
    expect_true(all(res@backend@modCount > sps_h5@backend@modCount))

    ## Applying the processing queue invalidated the original object!
    expect_error(as.list(sps_h5))
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
    expect_identical(as.list(res), as.list(sps_mem))

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
    expect_identical(as.list(res), as.list(sps_mem))
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

test_that(".compare_spectra, .compare_spectra_self work", {
    sps <- Spectra(sciex_hd5)[120:126]
    sps <- setBackend(sps, MsBackendDataFrame())

    res <- .compare_spectra(sps, sps)
    expect_true(ncol(res) == length(sps))
    expect_true(nrow(res) == length(sps))
    expect_equal(diag(res), rep(1, length(sps)))

    res_2 <- .compare_spectra(sps, sps[3])
    expect_true(ncol(res_2) == 1)
    expect_true(nrow(res_2) == length(sps))
    expect_identical(res_2[, 1], res[, 3])

    res_2 <- .compare_spectra(sps[5], sps)
    expect_true(ncol(res_2) == length(sps))
    expect_true(nrow(res_2) == 1)
    expect_identical(res_2[1, ], res[5, ])

    res_2 <- .compare_spectra_self(sps)
    expect_equal(dim(res), dim(res_2))
    expect_identical(diag(res), diag(res_2))
    expect_identical(res[!lower.tri(res)], res_2[!lower.tri(res_2)])

    cor_fun <- function(x, y, ...) {
        cor(x[, 2], y[, 2], use = "pairwise.complete.obs")
    }
    res <- .compare_spectra(sps[1], sps[1], FUN = cor_fun)
    expect_true(res[1, 1] == 1)
    res <- .compare_spectra(sps[1], sps[2], FUN = cor_fun)
    res_2 <- .compare_spectra(sps[1], sps[2])
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

    res <- .combine_spectra(sps, FUN = combinePeaks, tolerance = 0.1)
    expect_true(length(res) == 1)
    expect_equal(res$mz[[1]], c(mean(c(12, 12.1)), mean(c(14, 14.1, 14.15)),
                                mean(c(34, 34.1)), 45, mean(c(56, 56.1))))
    expect_equal(res$intensity[[1]], c(mean(c(10, 12)), mean(c(20, 11, 22)),
                                       mean(c(21, 32)), 30, mean(c(40, 31))))
    res <- .combine_spectra(sps, FUN = combinePeaks, tolerance = 0.1,
                            mzFun = max, intensityFun = median)
    expect_true(length(res) == 1)
    expect_equal(res$mz[[1]], c(max(c(12, 12.1)), max(c(14, 14.1, 14.15)),
                                max(c(34, 34.1)), 45, max(c(56, 56.1))))
    expect_equal(res$intensity[[1]], c(median(c(10, 12)), median(c(20, 11, 22)),
                                       median(c(21, 32)), 30, median(c(40, 31))))

    ## See if it works with MsBackendMzR
    sps <- Spectra(sciex_mzr)
    res <- .combine_spectra(sps, tolerance = 0.1, FUN = combinePeaks)
    expect_true(length(res) == 2)
    expect_true(is(res, "Spectra"))
    expect_true(class(res@backend) == "MsBackendDataFrame")
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
    expect_true(class(res@backend) == "MsBackendDataFrame")
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
