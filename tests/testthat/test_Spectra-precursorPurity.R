sps_dda_mem <- setBackend(sps_dda, MsBackendMemory())

test_that("precursorPurity works", {
    expect_error(precursorPurity(4), "needs to be a 'Spectra' object")

    res <- precursorPurity(sps_dda, tolerance = 0.1)
    res_2 <- precursorPurity(sps_dda, useReportedIsolationWindow = TRUE)
    ## purity is expected to be smaller for larger isolation window
    expect_true(length(which(res_2 < res)) > length(which(res_2 > res)))
})

test_that(".precursorPurity works", {
    ## error
    expect_error(.precursorPurity(sps_dda[c(34, 1, 34)]), "not increasingly")
    
    ## Small tolerance
    res <- .precursorPurity(sps_dda, tolerance = 0.05, ppm = 10)
    expect_true(is.numeric(res))
    expect_equal(length(res), length(sps_dda))
    expect_true(all(is.na(res[msLevel(sps_dda) == 1L])))
    expect_true(all(!is.na(res[msLevel(sps_dda) == 2L])))
    
    ## Large tolerance
    res_2 <- .precursorPurity(sps_dda, tolerance = 1, ppm = 10)
    expect_true(is.numeric(res_2))
    expect_equal(length(res_2), length(sps_dda))
    expect_true(all(is.na(res_2[msLevel(sps_dda) == 1L])))
    expect_true(all(!is.na(res_2[msLevel(sps_dda) == 2L])))
    ## Purity is generally lower for larger isolation window
    expect_true(length(which(res_2 < res)) > length(which(res_2 > res)))

    ## reported isolation window
    res_3 <- .precursorPurity(sps_dda, useIsolationWindow = TRUE)
    expect_false(all(res_3 == res_2))
    expect_false(all(res_3 == res))
    expect_true(length(which(res_3 < res)) > length(which(res_3 > res)))

    ## No MS2 data
    tmp <- filterMsLevel(sps_dda, 1L)
    expect_true(all(is.na(.precursorPurity(tmp, tolerance = 0.3))))
})
