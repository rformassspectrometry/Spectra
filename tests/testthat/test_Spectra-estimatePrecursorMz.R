test_that("estimatePrecursorMz works", {
    expect_error(estimatePrecursorMz(3), "'Spectra' object")
})

test_that(".adjust_dda_precursor_mz works", {
    tmp <- setBackend(sps_dda, MsBackendMemory())
    res <- .adjust_dda_precursor_mz(tmp, tolerance = 0.5)
    expect_true(length(res) == length(tmp))
    expect_false(identical(res, precursorMz(tmp)))
    expect_equal(is.na(res), is.na(precursorMz(tmp)))

    res_2 <- estimatePrecursorMz(tmp, tolerance = 0.5)
    expect_identical(res, res_2)

    tmp <- tmp[sample(seq_along(tmp))]
    expect_error(estimatePrecursorMz(tmp, "not increasingly"))
})
