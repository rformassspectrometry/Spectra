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
    ## expect_error(addProcessingStep(tst, "4"))
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
