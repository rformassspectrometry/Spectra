test_that(".filterSpectraHierarchy works", {
    msLevel <- c(1, 1, 2, 2, 2, 3, 3, 3, 4, 4)
    acquisitionNum <- 1:10
    precursorScanNum <- c(0, 0, 200, 2, 2, 4, 4, 5, 6, 20)

    expect_error(.filterSpectraHierarchy(1:3, 2, 1), "have to be the same")

    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 1)), 1)
    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 2)), c(2, 4:9))
    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 1:2)), c(1:2, 4:9))
    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 8)), c(2, 5, 8))
    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 9)), c(2, 4, 6, 9))
    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 10)), 10)
    expect_equal(which(.filterSpectraHierarchy(
        acquisitionNum, precursorScanNum, 11)), integer())
})

test_that("sanitize_file_name works", {
    a <- c("<memory>", "/other/path")
    expect_warning(res <- sanitize_file_name(a), "file")
    expect_equal(basename(res)[1], "memory")
    expect_equal(basename(res)[2], "path")
})

test_that(".values_match_mz works", {
    pmz <- c(12.4, 15, 3, 12.4, 3, 1234, 23, 5, 12.4, NA, 3)
    mz <- c(200, 12.4, 3)

    res <- .values_match_mz(pmz, mz)
    expect_true(all(pmz[res] %in% mz))
    expect_false(any(pmz[-res] %in% mz))

    pmz <- rev(pmz)
    res <- .values_match_mz(pmz, mz)
    expect_true(all(pmz[res] %in% mz))
    expect_false(any(pmz[-res] %in% mz))

    res <- .values_match_mz(c(NA, NA), mz)
    expect_identical(res, integer())

    res <- .values_match_mz(pmz, c(NA, 3))
    expect_true(all(pmz[res] == 3))
})

test_that(".long_spectra_data2 works", {
    sd <- as.data.frame(spectraData(sciex_mzr))
    lf <- .long_spectra_data2(sd, peaksVariables(sciex_mzr))
    expect_true(is.data.frame(lf))
    expect_equal(colnames(sd), colnames(lf))
    expect_true(is.numeric(lf$mz))
    expect_true(is.numeric(lf$intensity))
    expect_equal(nrow(lf), sum(lengths(sd$mz)))

    lf <- .long_spectra_data2(sd[, c("msLevel", "rtime")],
                              peaksVariables = character())
    expect_equal(lf, sd[, c("msLevel", "rtime")])

    lf <- .long_spectra_data2(sd[, c("rtime", "mz")], "mz")
    expect_true(is.data.frame(lf))
    expect_equal(colnames(lf), c("rtime", "mz"))
})

test_that(".long_spectra_data3 works", {
    a <- as.data.frame(spectraData(sciex_mzr, c("msLevel", "rtime")))
    b <- peaksData(sciex_mzr)
    lf <- .long_spectra_data3(a, b)
    expect_true(is.data.frame(lf))
    expect_equal(colnames(lf), c("msLevel", "rtime", "mz", "intensity"))
    expect_equal(nrow(lf), sum(lengths(b) / 2))
})
