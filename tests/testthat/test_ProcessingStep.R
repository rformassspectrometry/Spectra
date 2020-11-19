test_that("ProcessingStep constructor", {
    ps <- new("ProcessingStep", FUN = "mean", ARGS = list(na.rm = TRUE))
    expect_match(capture.output(show(ps))[1],
                 "Object of class \"ProcessingStep\"")
    expect_error(new("ProcessingStep", FUN = "aaaa"))
    ps <- ProcessingStep(mean, list(na.rm = TRUE))
    expect_true(validObject(ps))
})

test_that("ProcessingStep executeProcessingStep", {
    ps <- ProcessingStep("sum", list(c(1, 4, 5)))
    expect_identical(executeProcessingStep(ps), 10)

    ps <- ProcessingStep("sum", list(c(1, 4, 5, NA)))
    ## Pass optional arguments to ...
    expect_identical(executeProcessingStep(ps, na.rm = TRUE), 10)
})

test_that(".cat_fun works", {
    res <- .cat_fun("sum")
    expect_equal(res, "sum")
    res <- .cat_fun(sum)
    expect_equal(res, "user-provided function")
})

test_that("show,ProcessingStep works", {
    ps <- ProcessingStep(sum, list(1:4))
    show(ps)

    ps <- ProcessingStep(function(z) z + 4)
    show(ps)
    ps <- ProcessingStep(sum, list(a = sqrt))
    show(ps)
})
