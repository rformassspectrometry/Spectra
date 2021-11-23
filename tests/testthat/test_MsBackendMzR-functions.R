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

    expect_equal(.mzR_peaks(fls[1L], hdr$scanIndex),
                 .mzR_peaks(fls[1L], hdr$scanIndex, readHeader = TRUE))
})

test_that(".pattern_to_cv works", {
    expect_equal(.pattern_to_cv("unknown"), NA_character_)
    expect_equal(.pattern_to_cv("peak picking"), "MS:1000035")
    expect_equal(.pattern_to_cv("centroid"), "MS:1000035")
    expect_equal(.pattern_to_cv("Alignment/retention time adjustment"),
                 "MS:1000745")
})

test_that(".guess_software_processing works", {
    df <- DataFrame(msLevel = c(1L, 2L, 2L), scanIndex = 4:6)
    df$mz <- list(c(1.1, 1.3, 1.5), c(4.1, 5.1), c(1.6, 1.7, 1.8, 1.9))
    df$intensity <- list(c(45.1, 34, 12), c(234.4, 1333), c(42.1, 34.2, 65, 6))

    sps <- Spectra(df)
    res <- .guess_software_processing(sps)
    expect_equal(unname(res[[1]][1]), "Spectra")
    expect_equal(unname(res[[1]][2]), as.character(packageVersion("Spectra")))
    expect_equal(unname(res[[1]][3]), "MS:-1")

    sps <- filterMsLevel(sps, 2L)
    res <- .guess_software_processing(sps)
    expect_equal(unname(res[[1]][4]), "MS:1001486")

    sps <- pickPeaks(sps)
    res <- .guess_software_processing(sps)
    expect_equal(unname(res[[1]][5]), "MS:1000035")
})

test_that(".write_ms_data_mzR works", {
    df <- DataFrame(msLevel = c(1L, 2L, 2L), scanIndex = 4:6)
    df$mz <- list(c(1.1, 1.3, 1.5), c(4.1, 5.1), c(1.6, 1.7, 1.8, 1.9))
    df$intensity <- list(c(45.1, 34, 12), c(234.4, 1333), c(42.1, 34.2, 65, 6))
    sps <- Spectra(df)

    fl <- tempfile()
    .write_ms_data_mzR(sps, fl)
    res <- Spectra(backendInitialize(MsBackendMzR(), file = fl))
    expect_equal(rtime(res), rtime(sps))
    expect_equal(mz(res), mz(sps))
    expect_equal(intensity(res), intensity(sps))
    expect_warning(.write_ms_data_mzR(sps, fl, copy = TRUE), "not found")

    ## With an mzML file AND copy option.
    sps <- Spectra(sciex_mzr)
    sps <- sps[20:50]
    sps <- filterRt(sps, c(1, 500))
    .write_ms_data_mzR(sps, fl)
    res <- Spectra(backendInitialize(MsBackendMzR(), file = fl))
    expect_equal(rtime(sps), rtime(res))
    expect_equal(mz(sps), mz(res))

    .write_ms_data_mzR(sps, fl, copy = TRUE)
    res <- Spectra(backendInitialize(MsBackendMzR(), file = fl))
    expect_equal(rtime(sps), rtime(res))
    expect_equal(mz(sps), mz(res))
})
