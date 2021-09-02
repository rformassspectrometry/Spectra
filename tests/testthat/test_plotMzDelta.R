context("test_plotMzDelta")

spd <- DataFrame(msLevel = rep(2L, 2))
spd$mz <- list(c(100, 200, 300), c(100, 200, 300))
spd$intensity <- list(c(100, 101, 102), c(100, 101, 102))
sp <- Spectra(spd)

test_that("computeMzDeltas works", {
    res <- computeMzDeltas(sp, percentage = 1, xlim = c(20, 1000))
    expect_identical(length(res), length(sp))
    expect_identical(res[[1]], res[[2]])
    expect_identical(res[[1]], c(100, 200, 100))

    res <- computeMzDeltas(sp, percentage = 1, xlim = c(20, 150))
    expect_identical(length(res), length(sp))
    expect_identical(res[[1]], res[[2]])
    expect_identical(res[[1]], c(100, 100))

    res <- computeMzDeltas(sp, percentage = 0.6, xlim = c(20, 1000))
    expect_identical(length(res), length(sp))
    expect_identical(res[[1]], res[[2]])
    expect_identical(res[[1]], 100)
})


test_that("computeMzDeltas works", {
    f <- msdata::proteomics(pattern = "TMT.+20141210.mzML.gz", full.names = TRUE)
    sp <- Spectra(f)
    d <- computeMzDeltas(sp[1:1000])
    vdiffr::expect_doppelganger(
                "plotMzDelta_1000",
                function() plotMzDelta(d))
})
