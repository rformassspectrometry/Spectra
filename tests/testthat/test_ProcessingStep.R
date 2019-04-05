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

test_that(".apply_processing_queue works", {
    inp <- list(1:5, 1:3, 5)
    expect_equal(.apply_processing_queue(inp), inp)
    res <- .apply_processing_queue(inp, list(ProcessingStep("sum")))
    expect_equal(res, list(sum(1:5), sum(1:3), 5))

    q <- list(ProcessingStep(function(x, y) x + y, ARGS = list(y = 3)),
              ProcessingStep(function(x, y) x - y, ARGS = list(y = 1)))
    res <- .apply_processing_queue(inp, q)
    expect_equal(res, list((1:5 + 2), (1:3 + 2), 7))
})
