test_that(".spectra_from_file_mzR works", {
    mzf <- mzR::openMSfile(sciex_file[1])
    hdr <- mzR::header(mzf)
    hdr$spIdx <- 1:nrow(hdr)
    mzR::close(mzf)

    expect_error(.spectra_from_file_mzR())
    res_all <- Spectra:::.spectra_from_file_mzR(sciex_file[1],
                                      hdr)
    expect_true(length(res_all) == 931)
    res <- .spectra_from_file_mzR(sciex_file[1],
                                  hdr[13, ])
    expect_true(nrow(res[[1]]) == 1650)
    expect_equal(res[[1]], res_all[[13]])
})

test_that("backendReadSpectra,BackendMzR works", {

    spd <- DataFrame(spIdx = c(1:931, 1:931),
                     fileIdx = rep(c(1, 2), each = 931))
    be <- BackendMzR()
    be <- Spectra:::backendInitialize(be, sciex_file, spd)

    ## Get spectra from a single file
    res_1 <- Spectra:::backendReadSpectra(be, spd[spd$fileIdx == 1, ])
    expect_true(all(vapply(res_1, is.data.frame, logical(1))))
    expect_true(all(vapply(res_1, function(z)
        all(colnames(z) == c("mz", "intensity")), logical(1))))
    res_2 <- Spectra:::backendReadSpectra(be, spd[spd$fileIdx == 2, ])
    expect_true(all(vapply(res_2, is.data.frame, logical(1))))

    res <- Spectra:::backendReadSpectra(be, spd)
    expect_equal(res, c(res_1, res_2))
})
