
test_that(".valid_ms_backend_files_from_file works", {
    expect_null(.valid_ms_backend_files_from_file(
        c("a", "b", "c"), c(1, 1, 2, 2, 3, 3, 3)))
    expect_match(.valid_ms_backend_files_from_file(
        c("a", "b", "c"), c(1, 1, 2, 2, 3, 3, 4)),
        "Index in 'fromFile'")
    expect_match(.valid_ms_backend_files_from_file(character(), 1:3),
                 "'files' can not be empty")
    expect_match(.valid_ms_backend_files_from_file("a", integer()),
                 "'fromFile' can not be empty")
})

test_that(".valid_spectra_data_required_columns works", {
    df <- DataFrame()
    expect_null(.valid_spectra_data_required_columns(df))
    df <- DataFrame(msLevel = 1L)
    expect_match(.valid_spectra_data_required_columns(df),
                 "Required column")
    df$dataStorage <- "some"
    expect_null(.valid_spectra_data_required_columns(df))
})

test_that(".valid_column_datatype works", {
    expect_null(.valid_column_datatype(DataFrame(a = 4)))
    expect_match(.valid_column_datatype(
        DataFrame(msLevel = "a", acquisitionNum = 2.3, mz = list(1:4))),
        "type: msLevel, acquisitionNum.")
    expect_null(.valid_column_datatype(
        DataFrame(msLevel = 1L, acquisitionNum = 2L, centroided = TRUE)))
})

test_that(".valid_mz_column works", {
    df <- DataFrame(msLevel = c(1L, 1L))
    df$mz <- IRanges::NumericList(1:3, 1:12)
    expect_null(.valid_mz_column(df))
    expect_null(.valid_mz_column(DataFrame(msLevel = 4L)))
    df$mz <- list(1:3, letters[1:4])
    expect_match(.valid_mz_column(df),
                 "mz column should be of type NumericList")
    df$mz <- IRanges::NumericList(1:4, c(3, 2, 5))
    expect_match(.valid_mz_column(df),
                 "sorted increasingly")
})

test_that(".valid_intensity_column works", {
    df <- DataFrame(msLevel = c(1L, 1L))
    df$intensity <- IRanges::NumericList(4, 2:5)
    expect_null(.valid_intensity_column(df))
    df$intensity <- list("g", TRUE)
    expect_match(.valid_intensity_column(df),
                 "should be of type NumericList")
})

test_that(".valid_intensity_mz_columns works", {
    be <- MsBackendDataFrame()
    expect_null(.valid_intensity_mz_columns(be@spectraData))
    be <- backendInitialize(be, DataFrame(fromFile = c(1L, 1L)))
    expect_null(.valid_intensity_mz_columns(be@spectraData))
    be@spectraData$mz <- list(1:3, 1:2)
    be@spectraData$intensity <- list(1:3, 2)
    expect_match(.valid_intensity_mz_columns(be@spectraData),
                 "Length of mz and intensity")
    be@spectraData$intensity <- list(1:3, 1:2)
    expect_null(.valid_intensity_mz_columns(be@spectraData))
})

test_that(".get_column works", {
    df <- DataFrame(msLevel = c(NA_integer_, NA_integer_),
                    fromFile = rep(1L, 2), other_col = rep("a", 2))
    res <- .get_column(df, column = "fromFile")
    expect_equal(res, c(1L, 1L))
    res <- .get_column(df, column = "other_col")
    expect_equal(res, c("a", "a"))
    res <- .get_column(df, column = "msLevel")
    expect_equal(res, c(NA_integer_, NA_integer_))
    res <- .get_column(df, column = "precursorCharge")
    expect_equal(res, c(NA_integer_, NA_integer_))
    expect_error(.get_column(df, column = "a"), "not available")
    df <- DataFrame()
    res <- .get_column(df, column = "msLevel")
    expect_equal(res, integer())
})

test_that(".sel_file works", {
    df <- DataFrame(msLevel = 1L,
                    dataOrigin = c("a", "a", "a", "a", "b", "c"))
    be <- backendInitialize(MsBackendDataFrame(), df)
    expect_error(.sel_file(be, dataStorage = 5), "should be an integer")
    expect_error(.sel_file(be, dataOrigin = 5), "should be an integer")
    res <- .sel_file(be)
    expect_identical(res, rep(TRUE, length(be)))
    dataStorage(be) <- c("a", "a", "b", "b", "c", "c")
    res <- .sel_file(be, dataStorage = c("c", "a"))
    expect_identical(res, c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
    res <- .sel_file(be, dataStorage = "z")
    expect_identical(res, rep(FALSE, length(be)))
    res <- .sel_file(be, dataStorage = 3)
    expect_identical(res, c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))
    res <- .sel_file(be, dataStorage = NA_character_)
    expect_identical(res, rep(FALSE, length(be)))
    expect_error(.sel_file(be, dataStorage = TRUE), "integer with the index")

    res <- .sel_file(be, dataOrigin = c("c", "a"))
    expect_identical(res, c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE))
    res <- .sel_file(be, dataOrigin = "z")
    expect_identical(res, rep(FALSE, length(be)))
    res <- .sel_file(be, dataOrigin = 3)
    expect_identical(res, c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))
    res <- .sel_file(be, dataOrigin = NA_character_)
    expect_identical(res, rep(FALSE, length(be)))
    expect_error(.sel_file(be, dataOrigin = TRUE), "integer with the index")
    res <- .sel_file(be, dataOrigin = c("b", "z"))
    expect_identical(res, c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
})

test_that(".combine_backend_data_frame works", {
    df <- DataFrame(msLevel = c(1L, 2L, 2L), fromFile = 1L,
                    rtime = as.numeric(1:3))
    df2 <- DataFrame(msLevel = c(2L, 1L), fromFile = 1L,
                     rtime = c(4.1, 5.2), scanIndex = 1:2)
    df3 <- DataFrame(msLevel = c(1L, 2L), fromFile = 1L,
                     precScanNum = 1L, other_col = "z")
    be <- backendInitialize(MsBackendDataFrame(), df)
    be2 <- backendInitialize(MsBackendDataFrame(), df2)
    be3 <- backendInitialize(MsBackendDataFrame(), df3)

    expect_equal(.combine_backend_data_frame(list(be)), be)
    expect_error(backendMerge(list(be, 4)), "backends of the same type")

    res <- .combine_backend_data_frame(list(be, be2, be3))
    expect_true(is(res, "MsBackendDataFrame"))
    expect_identical(res@spectraData$dataStorage, rep("<memory>", 7))
    expect_identical(dataStorage(res), rep("<memory>", 7))
    expect_identical(msLevel(res), c(1L, 2L, 2L, 2L, 1L, 1L, 2L))
    expect_identical(rtime(res), c(1:3, 4.1, 5.2, NA, NA))
    expect_identical(res@spectraData$other_col,
                     c(rep(NA_character_, 5), "z", "z"))
    expect_true(is(be3@spectraData$precScanNum, "integer"))

    ## One backend with and one without m/z
    df2$mz <- list(c(1.1, 1.2), c(1.1, 1.2))
    df2$intensity <- list(c(12.4, 3), c(123.4, 1))
    be2 <- backendInitialize(MsBackendDataFrame(), df2)
    res <- .combine_backend_data_frame(list(be, be2, be3))
    expect_equal(lengths(mz(res)), c(0, 0, 0, 2, 2, 0, 0))

    df2$pk_ann <- list(c("a", "b"), c("c", "d"))
    be4 <- backendInitialize(MsBackendDataFrame(), df2,
                             peaksVariables = c("mz", "intensity", "pk_ann"))
    expect_error(.combine_backend_data_frame(list(be, be2, be4)),
                 "different peaks variables")

    ## With different dataStorage
    be$dataStorage <- c("a", "a", "a")
    be3$dataStorage <- c("z", "b")

    res <- .combine_backend_data_frame(list(be, be2, be3))
    expect_identical(res$dataStorage,
                     c("a", "a", "a", "<memory>", "<memory>", "z", "b"))
    expect_identical(res@spectraData$dataStorage,
                     c("a", "a", "a", "<memory>", "<memory>", "z", "b"))
    expect_identical(rtime(res), c(1:3, 4.1, 5.2, NA, NA))
})

test_that(".subset_backend_data_frame works", {
    be <- MsBackendDataFrame()
    expect_identical(.subset_backend_data_frame(be), be)

    be <- as(sciex_hd5, "MsBackendDataFrame")
    expect_identical(.subset_backend_data_frame(be), be)
    res <- .subset_backend_data_frame(be, 13)
    expect_true(is(res, "MsBackendDataFrame"))
    expect_equal(res$rtime, be$rtime[13])
})

test_that(".peaks_variables works", {
    b <- MsBackendDataFrame()
    expect_equal(.peaks_variables(b), c("mz", "intensity"))
    b@peaksVariables <- c("a", "b")
    expect_equal(.peaks_variables(b), c("a", "b"))
})

test_that(".valid_peaks_variable_columns works", {
    tmp <- data.frame(rtime = 1:3)
    tmp$mz <- list(1:3, 1:2, 1:4)
    tmp$intensity <- list(1:3, 1:2, 1:2)
    expect_match(.valid_peaks_variable_columns(tmp, c("mz", "intensity")),
                 "differ")
    tmp$intensity <- NULL
    expect_equal(.valid_peaks_variable_columns(tmp, c("mz", "intensity")), NULL)
    tmp$intensity <- list(1:3, 1:2, 1:4)
    tmp$other <- list(1:3, 1:2, 1:4)
    expect_equal(.valid_peaks_variable_columns(
        tmp, c("mz", "intensity", "other")), NULL)
    tmp$other <- list("a", "b", "c")
    expect_match(.valid_peaks_variable_columns(
        tmp, c("mz", "intensity", "other")), "differ")

})
