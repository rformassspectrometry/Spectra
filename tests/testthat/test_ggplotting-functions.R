#' # each base::graphics plot function must be wrapped by an anonymous function
#' # that could be called by `vdiffr::expect_doppelganger()`
#' Run devtools::test_active_file(file = "tests/testthat/test_ggplotting-functions.R")

## context("test_ggplotting-functions")

ints <- list(c(4.3412, 12, 8, 34, 23.4),
             c(8, 25, 16, 32))
mzs <- list(c(13.453421, 43.433122, 46.6653553, 129.111212, 322.24432),
            c(13.452, 43.5122, 129.112, 322.245))
df <- DataFrame(msLevel = c(1L, 1L), rtime = c(123.12, 124))
df$mz <- mzs
df$intensity <- ints
s <- Spectra(df)

test_that("ggplotSpectra works", {
    vdiffr::expect_doppelganger(
                "ggplotSpectra-color-each",
                ggplotSpectra(s, col = c("green", "blue"))
            )
    vdiffr::expect_doppelganger(
                "ggplotSpectra-color-peaks",
                ggplotSpectra(s, col = list(1:5, 1:4))
            )
    vdiffr::expect_doppelganger(
                "ggplotSpectra-color-peaks-label",
                ggplotSpectra(
                               s, labels = function(z) (mz(z)),
                               labelAngle = -30, labelHjust = 0,
                               col = list(1:5, 1:4))
            )
    vdiffr::expect_doppelganger(
                "ggplotSpectra-color-peaks-label-labelCol",
                ggplotSpectra(
                               s, labels = function(z) (mz(z)),
                               labelAngle = -30, labelHjust = 0,
                               col = list(1:5, 1:4), labelCol = "red")
            )
    vdiffr::expect_doppelganger(
                "ggplotSpectra-asp05",
                ggplotSpectra(s, asp = 1/2)
            )
    vdiffr::expect_doppelganger(
                "ggplotSpectra-asp2",
                ggplotSpectra(s, asp = 2)
            )
})

test_that("ggplotSpectraOverlay works", {
    vdiffr::expect_doppelganger(
                "ggplotSpectraOverlay-basic",
                ggplotSpectraOverlay(s, col = c("red", "green"))
            )

    vdiffr::expect_doppelganger(
                "ggplotSpectraOverlay-xlim",
                ggplotSpectraOverlay(s, xlim = c(0, 500))
            )

    vdiffr::expect_doppelganger(
                "ggplotSpectraOverlay-no-axes",
                ggplotSpectraOverlay(
                               s, axes = FALSE,
                               labels = function(z) mz(z))
            )
})

test_that("ggplotSpectraMirror works", {
    vdiffr::expect_doppelganger(
                "ggplotSpectraMirror-plain",
                ggplotSpectraMirror(s[1], s[2], main = "Comparison"))
    vdiffr::expect_doppelganger(
                "ggplotSpectraMirror-same",
                ggplotSpectraMirror(s[1], s[1], ppm = 0,
                                             tolerance = 0,
                                             frame.plot = FALSE))
    vdiffr::expect_doppelganger(
                "ggplotSpectraMirror-match-color",
                ggplotSpectraMirror(s[2], s[1], ppm = 0,
                                             tolerance = 0.1,
                                             labels = function(z) mz(z),
                                             matchCol = "red",
                                             matchLwd = 2, axes = FALSE,
                                             matchPch = 17))
    vdiffr::expect_doppelganger(
                "ggplotSpectraMirror-match-color-labelCol",
                ggplotSpectraMirror(s[2], s[1], ppm = 0,
                                             tolerance = 0.1,
                                             labels = function(z) mz(z),
                                             matchCol = "red",
                                             matchLwd = 2, axes = FALSE,
                                             matchPch = 17,
                                             labelCol = "blue"))
    expect_error(ggplotSpectraMirror(s), "have to be of length")
    expect_error(ggplotSpectraMirror(s[1], s[1], labels = list(c("a"))),
                 "occurs either because")
})

test_that(".plot_single_spectrum_ggplot works", {
    vdiffr::expect_doppelganger(
                "plot_single_spectrum_ggplot-basic",
                .plot_single_spectrum_ggplot(s[1])
            )
    vdiffr::expect_doppelganger(
                "plot_single_spectrum_ggplot-xlim",
                .plot_single_spectrum_ggplot(s[1], xlim = c(0, 700),
                                                ylim = c(0, 400))
            )

    vdiffr::expect_doppelganger(
                "plot_single_spectrum_ggplot-labels",
                .plot_single_spectrum_ggplot(
                              s[1], labels = format(mz(s)[[1]], digits = 4),
                              labelAngle = -30, labelHjust = 0)
            )
})

test_that("ggplotSpectra works with single peak spectrum", {
    df <- DataFrame(rtime = 132.2, msLevel = 1L)
    df$mz <- list(123)
    df$intensity <- list(4000)
    s <- Spectra(df)
    vdiffr::expect_doppelganger(
                "spectrum_single_peak_with_label",
                ggplotSpectra(s, labels = "long label",
                            labelAngle = -30, labelHjust = 0)
            )
})
