test_that(".valid_processing_queue works", {
    expect_true(is.character(.valid_processing_queue(list(3, 5))))
    lst <- list(ProcessingStep(mean), ProcessingStep("max"))
    expect_true(is.null(.valid_processing_queue(lst)))
})

test_that(".combine_data_frame works", {
    a <- DataFrame(col_1 = 1:5, col_2 = letters[1:5])
    b <- DataFrame(col_3 = 5:1, col_2 = 11:15)
    res <- .combine_data_frame(a, b)
    expect_true(ncol(res) == 3)
    expect_equal(colnames(res), c("col_1", "col_2", "col_3"))
    expect_equal(res$col_2, letters[1:5])
    res <- .combine_data_frame(b, a)
    expect_equal(res$col_2, 11:15)
    expect_error(.combine_data_frame(a, b[1:4, ]), "'x' and 'y'")

    res <- .combine_data_frame(a, DataFrame())
    expect_equal(res, a)
})

test_that("addProcessingStep works", {
    tst <- Spectra()
    tst <- addProcessingStep(tst, mean)
    expect_true(length(tst@processingQueue) == 1)
    expect_error(addProcessingStep(tst, "4"))
    tst <- addProcessingStep(tst, function(z, t) z * t, t = 4)
    expect_true(length(tst@processingQueue) == 2)
})

test_that(".remove_peaks works", {
    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = 1:length(int), intensity = int)
    res <- .remove_peaks(x, 1L, centroided = FALSE)
    expect_equal(res, x)
    res <- .remove_peaks(x, 1L, centroided = FALSE, t = 3)
    int_exp <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 0, 0, 0, 0,
                 1, 5, 10, 5, 1)
    expect_equal(res[, "intensity"], int_exp)

    res <- .remove_peaks(x, 1L, centroided = TRUE)
    int_exp <- c(0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 3, 10, 6, 2, 0, 0, 0, 2, 0, 0,
                 0, 5, 10, 5, 0)
    expect_equal(res[, "intensity"], int_exp)

    res <- .remove_peaks(x, 1L, centroided = TRUE, msLevel. = 2L)
    expect_equal(res[, "intensity"], int)
})

test_that(".clean_peaks works", {
    int <- c(0, 1, 2, 3, 1, 0, 0, 0, 0, 1, 3, 10, 6, 2, 1, 0, 1, 2, 0,
             0, 1, 5, 10, 5, 1)
    x <- cbind(mz = 1:length(int), intensity = int)

    res <- .clean_peaks(x, 1L, all = FALSE)
    expect_true(is.matrix(res))
    expect_equal(res[, "intensity"], c(0, 1, 2, 3, 1, 0, 0, 1, 3, 10, 6, 2,
                                       1, 0, 1, 2, 0, 0, 1, 5, 10, 5, 1))
    res <- .clean_peaks(x, 1L, all = TRUE)
    expect_equal(res[, "intensity"], c(1, 2, 3, 1, 1, 3, 10, 6, 2, 1, 1, 2,
                                       1, 5, 10, 5, 1))

    expect_equal(.clean_peaks(x, 2L, msLevel. = 1L), x)
    expect_equal(.clean_peaks(x, 1L, msLevel. = 2L), x)
})

test_that(".apply_processing_queue works", {
    inp <- list(1:5, 1:3, 5)
    expect_equal(Spectra:::.apply_processing_queue(inp), inp)
    res <- Spectra:::.apply_processing_queue(inp, msLevel = rep(0, 3),
                                             centroided = rep(FALSE, 3),
                                             list(ProcessingStep("sum")))
    expect_equal(res, list(sum(1:5), sum(1:3), 5))

    q <- list(ProcessingStep(function(x, y, ...) x + y, ARGS = list(y = 3)),
              ProcessingStep(function(x, y, ...) x - y, ARGS = list(y = 1)))
    res <- Spectra:::.apply_processing_queue(inp, msLevel = rep(0, 3),
                                   centroided = rep(FALSE, 3), q)
    expect_equal(res, list((1:5 + 2), (1:3 + 2), 7))

    be <- sciex_mzr
    pq <- list(ProcessingStep(.remove_peaks, list(t = 50000)))
    res <- .apply_processing_queue(peaks(be), msLevel(be),
                                   rep(TRUE, length(be)), pq)
    expect_true(all(vapply(res, function(z) all(z[z[, 2] > 0, 2] > 50000),
                           logical(1))))

    ## TODO: test with processing step of length 2.
})

test_that(".peaksapply works", {
    sps <- Spectra(backend = sciex_mzr)
    res <- .peaksapply(sps, FUN = .remove_peaks, t = 50000)
    expect_true(is.list(res))
    expect_equal(length(res), length(sps))
    expect_true(all(vapply(res, is.matrix, logical(1))))

    sps@processingQueue <- list(ProcessingStep(.remove_peaks, list(t = 50000)))
    res_2 <- .peaksapply(sps)
    expect_equal(res, res_2)

    ## TODO: test with processing step of length 2
})
