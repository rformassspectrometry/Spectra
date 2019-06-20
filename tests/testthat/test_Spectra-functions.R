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

test_that(".peaks_remove works", {
    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = 1:length(int), intensity = int)
    res <- .peaks_remove(x, 1L, centroided = FALSE)
    expect_equal(res, x)
    res <- .peaks_remove(x, 1L, centroided = FALSE, t = 3)
    int_exp <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 0, 0, 0, 0,
                 1, 5, 10, 5, 1)
    expect_equal(res[, "intensity"], int_exp)

    res <- .peaks_remove(x, 1L, centroided = TRUE)
    int_exp <- c(0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 3, 10, 6, 2, 0, 0, 0, 2, 0, 0,
                 0, 5, 10, 5, 0)
    expect_equal(res[, "intensity"], int_exp)

    res <- .peaks_remove(x, 1L, centroided = TRUE, msLevel = 2L)
    expect_equal(res[, "intensity"], int)
})

test_that(".peaks_clean works", {
    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = 1:length(int), intensity = int)

    res <- .peaks_clean(x, 1L, all = FALSE)
    expect_true(is.matrix(res))
    expect_equal(res[, "intensity"], c(0, 1, 2, 3, 1, 0, 0, 1, 3, 10, 6, 2,
                                       1, 0, 1, 2, 0, 0, 1, 5, 10, 5, 1))
    res <- .peaks_clean(x, 1L, all = TRUE)
    expect_equal(res[, "intensity"], c(1, 2, 3, 1, 1, 3, 10, 6, 2, 1, 1, 2,
                                       1, 5, 10, 5, 1))

    expect_equal(.peaks_clean(x, 2L, msLevel = 1L), x)
    expect_equal(.peaks_clean(x, 1L, msLevel = 2L), x)
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

test_that(".peaks_bin works", {
    s1 <- cbind(mz = 1:5, intensity = 1:5)
    s2 <- cbind(mz = 1:5 + 0.1, intensity = 1:5)
    r1 <- cbind(mz = 1:5 + 0.5, intensity = 1:5)
    r2 <- cbind(mz = c(2, 4, 6), intensity = c(3, 7, 5))
    r3 <- cbind(mz = c(2, 4, 6), intensity = c(1.5, 3.5, 5))
    r31 <- cbind(mz = c(2, 4, 6), intensity = c(1.5, 3.5, 5))
    r4 <- cbind(mz = c(1, 3, 5), intensity = c(1, 5, 9))
    expect_equal(.peaks_bin(s1, 1, binSize = 1), r1)
    expect_equal(.peaks_bin(s1, 1, binSize = 2), r2)
    expect_equal(.peaks_bin(s1, 1, binSize = 2, fun = mean), r3)
    expect_equal(.peaks_bin(s1, 1, breaks = seq(0, 7, by = 2)), r4)
    expect_equal(.peaks_bin(s2, 1, binSize = 1), r1)
    expect_equal(.peaks_bin(s2, 1, binSize = 2, fun = mean), r31)
    expect_equal(.peaks_bin(s2, 1, breaks = seq(0, 7, by = 2)), r4)
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
