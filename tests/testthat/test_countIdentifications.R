spdf <- data.frame(msLevel = as.integer(rep(rep(c(1, 2), c(1, 5)), 3)),
                   acquisitionNum = 1:18)
spdf$precScanNum <- ifelse(spdf$msLevel == 1, NA,
                           rep(spdf$acquisitionNum[spdf$msLevel == 1], each = 6))

test_that("countIdentifications() works (1)", {
    sp <- Spectra(spdf)
    expect_error(countIdentifications(sp))
})


test_that("countIdentifications() works (2)", {
    spdf$sequence <- ifelse(spdf$msLevel == 1, NA, "PEPTIDE")
    sp <- Spectra(spdf)
    ans <- countIdentifications(sp)
    ## All 5 MS2 scans are identified
    expect_identical(ans$countIdentifications,
                     as.integer(rep(rep(c(5, 1), c(1, 5)), 3)))
    ## 4 MS2 scans are identified
    sp$sequence[which(msLevel(sp) == 1) + 1] <- NA
    ans <- countIdentifications(sp)
    expect_identical(ans$countIdentifications,
                     as.integer(rep(c(4, 0, 1, 1, 1, 1), 3)))
    ## 3 MS2 scans are identified
    sp$sequence[which(msLevel(sp) == 1) + 2] <- NA
    ans <- countIdentifications(sp)
    expect_identical(ans$countIdentifications,
                     as.integer(rep(c(3, 0, 0, 1, 1, 1), 3)))
    ## 2 MS2 scans are identified
    sp$sequence[which(msLevel(sp) == 1) + 3] <- NA
    ans <- countIdentifications(sp)
    expect_identical(ans$countIdentifications,
                     as.integer(rep(c(2, 0, 0, 0, 1, 1), 3)))
    ## 1 MS2 scan is identified
    sp$sequence[which(msLevel(sp) == 1) + 4] <- NA
    ans <- countIdentifications(sp)
    expect_identical(ans$countIdentifications,
                     as.integer(rep(c(1, 0, 0, 0, 0, 1), 3)))
    ## No MS2 scans are identified
    sp$sequence[which(msLevel(sp) == 1) + 5] <- NA
    ans <- countIdentifications(sp)
    expect_identical(ans$countIdentifications,
                     as.integer(rep(c(0, 0, 0, 0, 0, 0), 3)))
})


test_that("countIdentifications() works (3)", {
    spdf$sequence <- ifelse(spdf$msLevel == 1, NA, "PEPTIDE")
    sp <- Spectra(spdf)
    ## All 5 MS2 scans are identified but they don't match the MS1
    ## acquistion numbers
    sp$acquisitionNum[which(msLevel(sp) == 1)] <- 101:103
    ans <- countIdentifications(sp)
    expect_identical(ans$countIdentifications,
                     as.integer(rep(c(0, 1, 1, 1, 1, 1), 3)))
    ## No sequence and no matching acquistion numbers
    sp$sequence <- NA
    ans <- countIdentifications(sp)
    expect_identical(ans$countIdentifications, rep(0L, 18))
})
