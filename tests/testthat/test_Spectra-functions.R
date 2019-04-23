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
