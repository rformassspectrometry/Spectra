
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

test_that(".combine_backend_dframe works", {
    df <- DataFrame(msLevel = c(1L, 2L, 2L), fromFile = 1L,
                    rtime = as.numeric(1:3))
    df2 <- DataFrame(msLevel = c(2L, 1L), fromFile = 1L,
                     rtime = c(4.1, 5.2), scanIndex = 1:2)
    df3 <- DataFrame(msLevel = c(1L, 2L), fromFile = 1L,
                     precScanNum = 1L, other_col = "z")
    be <- backendInitialize(MsBackendDFrame(), df)
    be2 <- backendInitialize(MsBackendDFrame(), df2)
    be3 <- backendInitialize(MsBackendDFrame(), df3)

    expect_equal(.combine_backend_data_frame(list(be)), be)
    expect_error(backendMerge(list(be, 4)), "backends of the same type")

    res <- .combine_backend_data_frame(list(be, be2, be3))
    expect_true(is(res, "MsBackendDFrame"))
    expect_identical(res@spectraData$dataStorage, Rle(rep("<memory>", 7)))
    expect_identical(dataStorage(res), rep("<memory>", 7))
    expect_identical(msLevel(res), c(1L, 2L, 2L, 2L, 1L, 1L, 2L))
    expect_identical(rtime(res), c(1:3, 4.1, 5.2, NA, NA))
    expect_identical(res@spectraData$other_col,
                     Rle(c(rep(NA_character_, 5), "z", "z")))
    expect_true(is(be3@spectraData$precScanNum, "integer"))
    expect_true(is(res@spectraData$precScanNum, "Rle"))

    ## One backend with and one without m/z
    df2$mz <- list(c(1.1, 1.2), c(1.1, 1.2))
    df2$intensity <- list(c(12.4, 3), c(123.4, 1))
    be2 <- backendInitialize(MsBackendDFrame(), df2)
    res <- Spectra:::.combine_backend_data_frame(list(be, be2, be3))
    expect_equal(lengths(mz(res)), c(0, 0, 0, 2, 2, 0, 0))

    ## With different dataStorage
    be$dataStorage <- c("a", "a", "a")
    be3$dataStorage <- c("z", "b")

    res <- .combine_backend_data_frame(list(be, be2, be3))
    expect_identical(res$dataStorage,
                     c("a", "a", "a", "<memory>", "<memory>", "z", "b"))
    expect_identical(res@spectraData$dataStorage,
                     Rle(c("a", "a", "a", "<memory>", "<memory>", "z", "b")))
    expect_identical(rtime(res), c(1:3, 4.1, 5.2, NA, NA))
})
