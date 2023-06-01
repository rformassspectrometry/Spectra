test_df <- DataFrame(msLevel = c(1L, 2L, 2L), scanIndex = 4:6)
test_df$mz <- list(c(1.1, 1.3, 1.5), c(4.1, 5.1), c(1.6, 1.7, 1.8, 1.9))
test_df$intensity <- list(c(45.1, 34, 12), c(234.4, 1333), c(42.1, 34.2, 65, 6))

test_that("MsBackendMemory works", {
    res <- MsBackendMemory()
    expect_s4_class(res, "MsBackendMemory")
})

test_that(".df_pdata_column works", {
    tmp <- list(cbind(a = 1:3, b = 2:4), cbind(a = 1:10, b = 1:10),
                cbind(a = 2, b = 3))
    res <- .df_pdata_column(tmp, "a")
    expect_equal(res, list(1:3, 1:10, c(a = 2)))

    expect_error(.df_pdata_column(tmp, "c"), "peaks variable")
    res <- .df_pdata_column(tmp, "b")
    expect_equal(res, list(2:4, 1:10, c(b = 3)))
})

test_that(".df_peaks_columns_data_frame works", {
    df <- data.frame(a = 1:4, b = "b")
    expect_equal(.df_peaks_columns_data_frame(df), character())
    df$lst <- list(1:3, 1:10, 3:54, 1:4)
    expect_equal(.df_peaks_columns_data_frame(df), character())
    df$mz <- list(1:3, 1:10, 3:54, 1:4)
    expect_equal(.df_peaks_columns_data_frame(df), c("lst", "mz"))
    df$intensity <- list(1:3, 1:10, 3:54, 1:4)
    df$mz <- NULL
    expect_equal(.df_peaks_columns_data_frame(df), c("lst", "intensity"))
    df$mz <- list(1:3, 1:10, 3:54, 1:4)

    df <- DataFrame(df)
    df$Lst <- IRanges::NumericList(1:3, 1:10, 2:40, 5:12, compress = FALSE)
    expect_equal(.df_peaks_columns_data_frame(df),
                 c("lst", "intensity", "mz"))
    df$Lst <- IRanges::NumericList(1:3, 1:10, 3:54, 1:4, compress = FALSE)
    expect_equal(.df_peaks_columns_data_frame(df),
                 c("lst", "intensity", "mz", "Lst"))
})

test_that(".df_spectra_data works", {
    be <- new("MsBackendMemory")

    res <- .df_spectra_data(be)
    expect_true(is.data.frame(res))
    expect_true(all(colnames(res) %in% names(coreSpectraVariables())))

    be <- backendInitialize(be, test_df)
    res <- .df_spectra_data(be)
    expect_true(all(names(coreSpectraVariables()) %in% colnames(res)))

    expect_equal(res$msLevel, test_df$msLevel)
    expect_equal(res$scanIndex, test_df$scanIndex)
    expect_true(all(res$dataStorage == "<memory>"))
    expect_equal(res$mz, IRanges::NumericList(test_df$mz, compress = FALSE))
    expect_equal(res$intensity,
                 IRanges::NumericList(test_df$intensity, compress = FALSE))

    tmp <- test_df
    tmp$pk_anno <- list(c("a", "b", "c"), c("", "d"), letters[12:15])
    be <- backendInitialize(
        be, tmp, peaksVariables = c("mz", "intensity", "pk_anno"))
    expect_true(length(be@peaksDataFrame) == 3)
    res <- .df_spectra_data(be)
    expect_equal(res$msLevel, test_df$msLevel)
    expect_equal(res$scanIndex, test_df$scanIndex)
    expect_true(all(res$dataStorage == "<memory>"))
    expect_equal(res$mz, IRanges::NumericList(test_df$mz, compress = FALSE))
    expect_equal(res$intensity,
                 IRanges::NumericList(test_df$intensity, compress = FALSE))
    expect_equal(res$pk_anno, tmp$pk_anno)

    tmp$add_anno <- list(c(1:3), 1:2, 1:4)
    be <- backendInitialize(be, tmp)
    res <- .df_spectra_data(be)
    expect_equal(res$pk_anno, tmp$pk_anno)
    expect_equal(res$add_anno, tmp$add_anno)

    res <- .df_spectra_data(be, "mz")
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), "mz")
    expect_equal(res$mz, IRanges::NumericList(tmp$mz, compress = FALSE))
    res <- .df_spectra_data(be, "msLevel")
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), "msLevel")
    expect_equal(res$msLevel, tmp$msLevel)
    res <- .df_spectra_data(be, "rtime")
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), "rtime")
    expect_equal(res$rtime, rep(NA_real_, 3))
    res <- .df_spectra_data(be, "pk_anno")
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), "pk_anno")
    expect_equal(res$pk_anno, tmp$pk_anno)
})

test_that(".df_subset works", {
    be <- new("MsBackendMemory")
    expect_error(.df_subset(be, 1:3), "out of bounds")

    be <- backendInitialize(be, test_df)
    res <- .df_subset(be, c(2, 1, 1))
    expect_equal(res$msLevel, test_df$msLevel[c(2, 1, 1)])
    expect_equal(res$mz, IRanges::NumericList(test_df$mz[c(2, 1, 1)],
                                              compress = FALSE))
    expect_equal(res$scanIndex, test_df$scanIndex[c(2, 1, 1)])

    vals <- list(letters[1:3], letters[1:2], letters[1:4])
    be$peak_anno <- vals
    res <- .df_subset(be, c(3, 1))
    expect_equal(res$peak_anno, vals[c(3, 1)])
})

test_that(".df_combine works", {
    be <- new("MsBackendMemory")
    be <- backendInitialize(be, test_df)
    tmp <- list(be, be, be)
    res <- .df_combine(tmp)
    expect_s4_class(res, "MsBackendMemory")
    expect_equal(length(res), length(be) * 3)
    expect_equal(res$msLevel, rep(test_df$msLevel, 3))
    expect_equal(res@peaksData[1:3], be@peaksData[1:3])
    expect_equal(res@peaksData[4:6], be@peaksData[1:3])
    expect_equal(res@peaksData[7:9], be@peaksData[1:3])

    ## With empty peak matrix
    df2 <- data.frame(msLevel = c(2L, 1L), rtime = c(1.2, 1.4))
    be2 <- backendInitialize(be, df2)
    res <- .df_combine(list(be, be2))
    expect_s4_class(res, "MsBackendMemory")
    expect_equal(length(res), 5)
    expect_equal(res$msLevel, c(test_df$msLevel, df2$msLevel))
    expect_equal(res$rtime, c(NA, NA, NA, 1.2, 1.4))
    expect_equal(res@peaksData[1:3], be@peaksData[1:3])
    expect_equal(res@peaksData[4:5], be2@peaksData[1:2])

    ## With peaksDataFrame
    tmp <- test_df
    tmp$peak_ann <- list(letters[1:3], letters[1:2], letters[1:4])
    be3 <- backendInitialize(
        be, tmp, peaksVariables = c("mz", "intensity", "peak_ann"))
    expect_error(.df_combine(list(be, be3)), "peak variables")
    res <- .df_combine(list(be3, be3))
    expect_s4_class(res, "MsBackendMemory")
    expect_equal(length(res), length(be) * 2)
    expect_equal(res$msLevel, rep(test_df$msLevel, 2))
    expect_equal(res@peaksData[1:3], be@peaksData[1:3])
    expect_equal(res@peaksData[4:6], be@peaksData[1:3])
    expect_equal(res@peaksDataFrame[1:3], be3@peaksDataFrame[1:3])
    expect_equal(res@peaksDataFrame[4:6], be3@peaksDataFrame[1:3])

    be3 <- new("MsBackendDataFrame")
    be3 <- backendInitialize(be3, test_df)
    expect_error(.df_combine(list(be, be3)), "same type")
})

test_that(".check_peaks_data_value works", {
    expect_error(.check_peaks_data_value("a", 3), "list-like")
    lst <- list("a")
    expect_error(.check_peaks_data_value(lst, 3), "match length")
    expect_error(.check_peaks_data_value(lst, 1), "matrix")
    lst <- list(data.frame(a = 3, b = 4), data.frame(a = 4, d = 4))
    expect_error(.check_peaks_data_value(lst, 2), "same column")
    lst <- list(data.frame(a = 3, b = 4), data.frame(a = 4))
    expect_error(.check_peaks_data_value(lst, 2), "same column")
})
