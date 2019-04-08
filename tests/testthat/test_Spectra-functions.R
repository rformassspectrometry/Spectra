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

test_that("Spectra works", {
    mzs <- list(1:4, 2:6, 1:7)
    ints <- list(abs(rnorm(4)), abs(rnorm(5)), abs(rnorm(7)))
    sp <- Spectra()
    expect_true(validObject(sp))
    expect_equal(length(sp), 0)

    expect_error(Spectra(msLevel = "a"), "'msLevel' needs to be")
    expect_error(Spectra(msLevel = 1, acquisitionNum = c(1, 2)),
                 "Different lengths of parameters")
    expect_error(Spectra(msLevel = c(1L, 1L), scanIndex = 1:2, mz = mzs),
                 "Length of 'mz'")
    expect_error(Spectra(msLevel = c(1L, 1L), scanIndex = 1:2, mz = 1:2),
                 "wrong data type: mz")
    expect_error(Spectra(msLevel = c(1L, 1L), scanIndex = 1:2,
                         mz = list(c(2, 1), 1:3)),
                 "mz values have to be sorted")
    expect_error(Spectra(msLevel = 1L, intensity = 3),
                 "wrong data type: intensity")

    sp <- Spectra(msLevel = c(1L, 1L, 1L), mz = mzs, intensity = ints,
                  scanIndex = 1:3)
    expect_true(validObject(sp))
    expect_true(validObject(sp@backend))
    expect_true(nrow(sp@backend@spectraData) == 3)
})
