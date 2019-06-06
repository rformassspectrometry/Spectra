test_that(".mzR_header works", {
    sc <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
    expect_error(.mzR_header(sc), "should have length 1")
    expect_error(.mzR_header(), "should have length 1")
    hdr <- .mzR_header(sc[1])
    expect_true(inherits(hdr, "DataFrame"))
    expect_equal(nrow(hdr), 931)
    expect_equal(hdr$scanIndex, 1:931)

    fl <- dir(system.file("cdf", package = "msdata"), full.names = TRUE)
    hdr <- .mzR_header(fl)
    expect_true(inherits(hdr, "DataFrame"))
    expect_equal(nrow(hdr), 1278)
    expect_equal(hdr$scanIndex, 1:1278)
})

test_that(".mzR_peaks work", {
    fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
    hdr <- .mzR_header(fls[1])

    res <- .mzR_peaks(fls[1], 1)
    expect_true(is(res, "list"))
    expect_true(is.matrix(res[[1]]))
    expect_equal(colnames(res[[1]]), c("mz", "intensity"))

    res_all <- .mzR_peaks(fls[1], hdr$scanIndex)
    expect_true(is(res_all, "list"))
    expect_true(is.matrix(res_all[[1]]))
    expect_equal(colnames(res_all[[1]]), c("mz", "intensity"))
    expect_equal(res[[1]], res_all[[1]])

    res_13 <- .mzR_peaks(fls[1], 13)
    expect_true(is(res_13, "list"))
    expect_true(is.matrix(res_13[[1]]))
    expect_equal(colnames(res_13[[1]]), c("mz", "intensity"))
    expect_equal(res_13[[1]], res_all[[13]])

    expect_error(.mzR_peaks(fls, 13), "length 1")
})
