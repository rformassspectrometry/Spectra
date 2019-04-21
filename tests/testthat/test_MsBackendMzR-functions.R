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

test_that(".compress_spectra_data works", {
    df <- DataFrame()
    res <- .compress_spectra_data(df)
    expect_equal(df, res)
    df <- DataFrame(msLevel = c(NA_integer_, NA_integer_), a = "a", b = 1:2)
    res <- .compress_spectra_data(df)
    expect_equal(res$msLevel, Rle(NA_integer_, 2))
    expect_equal(res$a, Rle("a", 2))
    res$fromFile <- c(1L, 1L)
    res <- .compress_spectra_data(res)
    expect_equal(res$msLevel, Rle(NA_integer_, 2))
    expect_equal(res$a, Rle("a", 2))
    expect_equal(res$b, 1:2)
    expect_equal(res$fromFile, Rle(1L, 2))
    res$fromFile <- c(1L, 2L)
    res <- .compress_spectra_data(res)
    expect_equal(res$fromFile, Rle(1:2))
})

test_that(".uncompress_spectra_data works", {
    df <- DataFrame(msLevel = c(NA_integer_, NA_integer_),
                    fromFile = Rle(1L, 2), other_col = Rle("a", 2))
    res <- .uncompress_spectra_data(df)
    expect_true(is.integer(res$msLevel))
    expect_true(is.integer(res$fromFile))
    expect_true(is(res$other_col, "character"))
})

test_that(".get_rle_column works", {
    df <- DataFrame(msLevel = c(NA_integer_, NA_integer_),
                    fromFile = Rle(1L, 2), other_col = Rle("a", 2))
    res <- .get_rle_column(df, column = "fromFile")
    expect_equal(res, c(1L, 1L))
    res <- .get_rle_column(df, column = "other_col")
    expect_equal(res, c("a", "a"))
    res <- .get_rle_column(df, column = "msLevel")
    expect_equal(res, c(NA_integer_, NA_integer_))
    res <- .get_rle_column(df, column = "precursorCharge")
    expect_equal(res, c(NA_integer_, NA_integer_))
    expect_error(.get_rle_column(df, column = "a"), "not available")
    df <- DataFrame()
    res <- .get_rle_column(df, column = "msLevel")
    expect_equal(res, integer())
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
