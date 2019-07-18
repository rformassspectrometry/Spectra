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
    pks <- peaks(be)
    pq <- list(ProcessingStep(.peaks_remove, list(t = 50000)))
    res <- .apply_processing_queue(pks, msLevel(be),
                                   rep(TRUE, length(be)), pq)
    expect_true(all(vapply(res, function(z) all(z[z[, 2] > 0, 2] > 50000),
                           logical(1))))
    expect_equal(vapply(res, nrow, integer(1)), vapply(pks, nrow, integer(1)))

    ## Length 2
    pq <- c(pq, list(ProcessingStep(.peaks_clean, list(all = TRUE))))
    res <- .apply_processing_queue(pks, msLevel(be),
                                   rep(TRUE, length(be)), pq)
    expect_true(all(vapply(res, function(z) all(z[z[, 2] > 0, 2] > 50000),
                           logical(1))))
    expect_true(all(vapply(res, nrow, integer(1)) <
                    vapply(pks, nrow, integer(1))))
})

test_that(".peaksapply works", {
    sps <- Spectra(backend = sciex_mzr)
    res <- .peaksapply(sps, FUN = .peaks_remove, t = 50000)
    expect_true(is.list(res))
    expect_equal(length(res), length(sps))
    expect_true(all(vapply(res, is.matrix, logical(1))))

    ## Ensure that this works with arbitrary ordering of the factor f
    res2 <- .peaksapply(sps, FUN = .peaks_remove, t = 50000,
                        f = rep(1:2, length(sps)/2))
    expect_identical(res, res2)

    sps@processingQueue <- list(ProcessingStep(.peaks_remove, list(t = 50000)))
    res_2 <- .peaksapply(sps)
    expect_equal(res, res_2)

    res_3 <- .peaksapply(sps, FUN = .peaks_clean, all = TRUE)
    expect_true(all(vapply(res_3, nrow, integer(1)) <
                    vapply(res_2, nrow, integer(1))))
    expect_true(!any(vapply(res_3, function(z) any(z[, 2] == 0), logical(1))))

    sps@processingQueue <- c(sps@processingQueue,
                             list(ProcessingStep(.peaks_clean, list(all = TRUE))))
    res_4 <- .peaksapply(sps)
    expect_equal(res_3, res_4)
})

test_that("applyProcessing works", {
    ## Initialize required objects.
    sps_mzr <- filterRt(Spectra(sciex_mzr), rt = c(10, 20))
    ## Add processings.
    centroided(sps_mzr) <- TRUE
    sps_mzr <- removePeaks(sps_mzr, t = 5000)
    sps_mzr <- clean(sps_mzr, all = TRUE)
    expect_true(length(sps_mzr@processingQueue) == 2)
    expect_error(applyProcessing(sps_mzr), "is read-only")

    ## Create writeable backends.
    sps_mem <- setBackend(sps_mzr, backend = MsBackendDataFrame())
    sps_h5 <- setBackend(sps_mzr, backend = MsBackendHdf5Peaks(),
                         files = c(tempfile(), tempfile()),
                         f = rep(1, length(sps_mzr)))
    expect_true(length(sps_mem@processingQueue) == 2)
    expect_true(length(sps_h5@processingQueue) == 2)
    expect_identical(peaks(sps_mzr), peaks(sps_mem))
    expect_identical(peaks(sps_h5), peaks(sps_mem))

    ## MsBackendDataFrame
    res <- applyProcessing(sps_mem)
    expect_true(length(res@processingQueue) == 0)
    expect_true(length(res@processing) > length(sps_mem@processing))
    expect_identical(rtime(res), rtime(sps_mem))
    expect_identical(peaks(res), peaks(sps_mem))

    ## MsBackendHdf5Peaks
    res <- applyProcessing(sps_h5)
    expect_true(length(res@processingQueue) == 0)
    expect_true(length(res@processing) > length(sps_h5@processing))
    expect_identical(rtime(res), rtime(sps_mem))
    expect_identical(peaks(res), peaks(sps_mem))
    expect_true(all(res@backend@modCount > sps_h5@backend@modCount))

    ## Applying the processing queue invalidated the original object!
    expect_error(peaks(sps_h5))
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
    expect_identical(peaks(res), peaks(sps_mem))

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
    expect_identical(peaks(res), peaks(sps_mem))
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
    sps <- Spectra(sciex_hd5)[120:130]
    sps <- setBackend(sps, MsBackendDataFrame())

    res <- Spectra:::.compare_spectra(sps, sps)
    expect_true(ncol(res) == length(sps))
    expect_true(nrow(res) == length(sps))
    expect_equal(diag(res), rep(1, length(sps)))

    res_2 <- Spectra:::.compare_spectra(sps, sps[3])
    expect_true(ncol(res_2) == 1)
    expect_true(nrow(res_2) == length(sps))
    expect_identical(res_2[, 1], res[, 3])

    res_2 <- Spectra:::.compare_spectra(sps[5], sps)
    expect_true(ncol(res_2) == length(sps))
    expect_true(nrow(res_2) == 1)
    expect_identical(res_2[1, ], res[5, ])

    res_2 <- Spectra:::.compare_spectra_self(sps)
    expect_equal(dim(res), dim(res_2))
    expect_identical(diag(res), diag(res_2))
    expect_identical(res[!lower.tri(res)], res_2[!lower.tri(res_2)])
})
