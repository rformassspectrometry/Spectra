
## test_that(".valid_ms_backend_files_from_file works", {
##     expect_null(.valid_ms_backend_files_from_file(
##         c("a", "b", "c"), c(1, 1, 2, 2, 3, 3, 3)))
##     expect_match(.valid_ms_backend_files_from_file(
##         c("a", "b", "c"), c(1, 1, 2, 2, 3, 3, 4)),
##         "Index in 'fromFile'")
##     expect_match(.valid_ms_backend_files_from_file(character(), 1:3),
##                  "'files' can not be empty")
##     expect_match(.valid_ms_backend_files_from_file("a", integer()),
##                  "'fromFile' can not be empty")
## })

## test_that(".valid_spectra_data_required_columns works", {
##     df <- DataFrame()
##     expect_null(.valid_spectra_data_required_columns(df))
##     df <- DataFrame(msLevel = 1L)
##     expect_match(.valid_spectra_data_required_columns(df),
##                  "Required column")
##     df$fromFile <- 1L
##     expect_null(.valid_spectra_data_required_columns(df))
## })

## test_that(".valid_column_datatype works", {
##     expect_null(.valid_column_datatype(DataFrame(a = 4)))
##     expect_match(.valid_column_datatype(
##         DataFrame(msLevel = "a", acquisitionNum = 2.3, mz = list(1:4))),
##         "type: msLevel, acquisitionNum.")
##     expect_null(.valid_column_datatype(
##         DataFrame(msLevel = 1L, acquisitionNum = 2L, centroided = TRUE)))
## })

## test_that(".valid_mz_column works", {
##     df <- DataFrame(msLevel = c(1L, 1L))
##     df$mz <- list(1:3, 1:12)
##     expect_null(.valid_mz_column(df))
##     expect_null(.valid_mz_column(DataFrame(msLevel = 4L)))
##     df$mz <- list(1:3, letters[1:4])
##     expect_match(.valid_mz_column(df),
##                  "mz column should contain a list of numeric")
##     df$mz <- list(1:4, c(3, 2, 5))
##     expect_match(.valid_mz_column(df),
##                  "sorted increasingly")
## })

## test_that(".valid_intensity_column works", {
##     df <- DataFrame(msLevel = c(1L, 1L))
##     df$intensity <- list(4, 2:5)
##     expect_null(.valid_intensity_column(df))
##     df$intensity <- list("g", TRUE)
##     expect_match(.valid_intensity_column(df),
##                  "contain a list of numeric")
## })

## test_that(".valid_intensity_mz_columns works", {
##     be <- MsBackendDataFrame()
##     expect_null(.valid_intensity_mz_columns(be))
##     be <- backendInitialize(be, files = NA_character_,
##                             DataFrame(fromFile = c(1L, 1L)))
##     expect_null(.valid_intensity_mz_columns(be))
##     be@spectraData$mz <- list(1:3, 1:2)
##     be@spectraData$intensity <- list(1:3, 2)
##     expect_match(.valid_intensity_mz_columns(be),
##                  "Length of mz and intensity")
##     be@spectraData$intensity <- list(1:3, 1:2)
##     expect_null(.valid_intensity_mz_columns(be))
## })

## test_that(".get_spectra_data_column works", {
##     be <- MsBackendDataFrame()
##     expect_equal(.get_spectra_data_column(be, "rt"), numeric())
##     df <- DataFrame(fromFile = c(1L, 1L), scanIndex = c(1L, 2L), other_col = "a")
##     be <- backendInitialize(be, files = NA_character_, df)
##     expect_equal(.get_spectra_data_column(be, "fromFile"), c(1L, 1L))
##     expect_equal(.get_spectra_data_column(be, "scanIndex"), 1:2)
##     expect_equal(.get_spectra_data_column(be, "other_col"), c("a", "a"))
##     expect_equal(.get_spectra_data_column(be, "precScanNum"), c(NA_integer_,
##                                                                 NA_integer_))
##     expect_error(.get_spectra_data_column(be, c("a", "b")), "of length 1")
##     expect_error(.get_spectra_data_column(be, "something"), "No column")
## })

test_that(".compress_spectra_data works", {
    df <- DataFrame()
    res <- .compress_spectra_data(df)
    expect_equal(df, res)
    df <- DataFrame(msLevel = c(NA_integer_, NA_integer_), a = "a")
    res <- .compress_spectra_data(df)
    expect_equal(res$msLevel, Rle(NA_integer_, 2))
    res$fromFile <- c(1L, 1L)
    res <- .compress_spectra_data(res)
    expect_equal(res$msLevel, Rle(NA_integer_, 2))
    expect_equal(res$fromFile, Rle(1L, 2))
})

test_that(".uncompress_spectra_data works", {
    df <- DataFrame(msLevel = c(NA_integer_, NA_integer_),
                    fromFile = Rle(1L, 2), other_col = Rle("a", 2))
    res <- .uncompress_spectra_data(df)
    expect_true(is.integer(res$msLevel))
    expect_true(is.integer(res$fromFile))
    expect_true(is(res$other_col, "Rle"))
})

test_that(".initialize_spectra_data works", {
    res <- .initialize_spectra_data(DataFrame())
    expect_true(all(colnames(res) %in% names(.SPECTRA_DATA_COLUMNS)))
    df <- DataFrame(fromFile = c(1L, 1L, 2L), other_col = "a")
    res <- .initialize_spectra_data(df)
    expect_true(all(names(.SPECTRA_DATA_COLUMNS) %in% colnames(res)))
    expect_equal(res$other_col, c("a", "a", "a"))
    expect_equal(res$msLevel, Rle(NA_integer_, 3))
})

test_that(".is_class works", {
    expect_true(.is_class(5L, "integer"))
    expect_true(.is_class("a", "character"))
    expect_true(.is_class(Rle("a", 5), "character"))
    expect_true(.is_class(Rle(3.4, 2), "numeric"))
})

test_that(".valid_column_rle_datatype works", {
    df <- DataFrame(fromFile = c(1L, 2L), msLevel = Rle(NA_integer_, 2))
    expect_null(.valid_column_rle_datatype(df))
    df$msLevel <- c("a", "b")
    expect_match(.valid_column_rle_datatype(df), "expected: integer")
    df$msLevel <- Rle(1L, 2)
    df$rtime <- TRUE
    expect_match(.valid_column_rle_datatype(df), "expected: numeric")
})
