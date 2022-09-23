test_that(".df_pdata_column works", {
    tmp <- list(cbind(a = 1:3, b = 2:4), cbind(a = 1:10, b = 1:10),
                cbind(a = 2, b = 3))
    res <- .df_pdata_column(tmp, "a")
    expect_equal(res, list(1:3, 1:10, c(a = 2)))

    expect_error(.df_pdata_column(tmp, "c"), "peaks variable")
    res <- .df_pdata_column(tmp, "b")
    expect_equal(res, list(2:4, 1:10, c(b = 3)))
})

test_that(".get_peaks_columns_data_frame works", {
    df <- data.frame(a = 1:4, b = "b")
    expect_equal(.get_peaks_columns_data_frame(df), character())
    df$lst <- list(1:3, 1:10, 3:54, 1:4)
    expect_equal(.get_peaks_columns_data_frame(df), character())
    df$mz <- list(1:3, 1:10, 3:54, 1:4)
    expect_equal(.get_peaks_columns_data_frame(df), c("lst", "mz"))
    df$intensity <- list(1:3, 1:10, 3:54, 1:4)
    df$mz <- NULL
    expect_equal(.get_peaks_columns_data_frame(df), c("lst", "intensity"))
    df$mz <- list(1:3, 1:10, 3:54, 1:4)

    df <- DataFrame(df)
    df$Lst <- IRanges::NumericList(1:3, 1:10, 2:40, 5:12, compress = FALSE)
    expect_equal(.get_peaks_columns_data_frame(df),
                 c("lst", "intensity", "mz"))
    df$Lst <- IRanges::NumericList(1:3, 1:10, 3:54, 1:4, compress = FALSE)
    expect_equal(.get_peaks_columns_data_frame(df),
                 c("lst", "intensity", "mz", "Lst"))
})
