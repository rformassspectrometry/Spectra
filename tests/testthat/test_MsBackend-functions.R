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
    df$mz <- list(1:3, 1:12)
    expect_null(.valid_mz_column(df))
    expect_null(.valid_mz_column(DataFrame(msLevel = 4L)))
    df$mz <- list(1:3, letters[1:4])
    expect_match(.valid_mz_column(df),
                 "mz column should contain a list of numeric")
    df$mz <- list(1:4, c(3, 2, 5))
    expect_match(.valid_mz_column(df),
                 "sorted increasingly")
})

test_that(".valid_intensity_column works", {
    df <- DataFrame(msLevel = c(1L, 1L))
    df$intensity <- list(4, 2:5)
    expect_null(.valid_intensity_column(df))
    df$intensity <- list("g", TRUE)
    expect_match(.valid_intensity_column(df),
                 "contain a list of numeric")
})

test_that(".valid_ms_backend_files works", {
    expect_null(.valid_ms_backend_files(c("a", "b")))
    expect_match(.valid_ms_backend_files(c("a", "a")))
})

test_that(".valid_ms_backend_mod_count works", {
    expect_null(.valid_ms_backend_mod_count(1, "a"))
    expect_match(.valid_ms_backend_mod_count(integer(), "a"),
                 "Different number of")
})
