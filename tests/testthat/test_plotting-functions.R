#' # each base::graphics plot function must be wrapped by an anonymous function
#' # that could be called by `vdiffr::expect_doppelganger()`
#' # use `manage_cases()` to add new/verify changed plots
#' vdiffr::manage_cases(filter = "plotting")

ints <- list(c(4.3412, 12, 8, 34, 23.4),
             c(8, 25, 16, 32))
mzs <- list(c(13.453421, 43.433122, 46.6653553, 129.111212, 322.24432),
            c(13.452, 43.5122, 129.112, 322.245))
df <- DataFrame(msLevel = c(1L, 1L), rtime = c(123.12, 124))
df$mz <- mzs
df$intensity <- ints
s <- Spectra(df)

test_that("plotSpectra works", {
    vdiffr::expect_doppelganger(
                "plotSpectra-color-each",
                function() plotSpectra(s, col = c("green", "blue"))
            )
    vdiffr::expect_doppelganger(
                "plotSpectra-color-peaks",
                function() plotSpectra(s, col = list(1:5, 1:4))
            )
    vdiffr::expect_doppelganger(
                "plotSpectra-color-peaks-label",
                function() plotSpectra(
                               s, labels = function(z) unlist(mz(z)),
                               labelPos = 2, labelOffset = 0.1, labelSrt = -30,
                               col = list(1:5, 1:4))
            )

})

test_that("plotSpectraOverlay works", {
    vdiffr::expect_doppelganger(
                "plotSpectraOverlay-basic",
                function() plotSpectraOverlay(s, col = c("red", "green"))
            )

    vdiffr::expect_doppelganger(
                "plotSpectraOverlay-xlim",
                function() plotSpectraOverlay(s, xlim = c(0, 500))
            )

    vdiffr::expect_doppelganger(
                "plotSpectraOverlay-no-axes",
                function() plotSpectraOverlay(
                               s, axes = FALSE,
                               labels = function(z) mz(z)[[1L]])
            )
})

test_that(".plot_single_spectrum works", {
    vdiffr::expect_doppelganger(
                "plot_single_spectrum-basic",
                function().plot_single_spectrum(s[1])
            )
    vdiffr::expect_doppelganger(
                "plot_single_spectrum-xlim",
                function().plot_single_spectrum(s[1], xlim = c(0, 700),
                                                ylim = c(0, 400))
            )

    vdiffr::expect_doppelganger(
                "plot_single_spectrum-labels",
                function().plot_single_spectrum(
                              s[1], labels = format(mz(s)[[1]], digits = 4),
                              labelPos = 2, labelOffset = 0.1, labelSrt = -30)
            )

    vdiffr::expect_doppelganger(
                "plot_single_spectrum-labels-ass",
                function() {
                    .plot_single_spectrum(
                        s[1], labels = format(mz(s)[[1]], digits = 4),
                        labelPos = 2, labelOffset = 0.1, labelSrt = -30)
                    .plot_single_spectrum(s[1], add = TRUE, col = 2, type = "p")
                }
            )

    .plot_single_spectrum(
        s[1], labels = function(z) mz(z)[[1L]],
        main = "Spectrum with labels",
        labelSrt = -30, labelPos = 2)

    plot(3, 3, xlim = c(0, 1000), ylim = c(0, 1000))
    .plot_single_spectrum(
        s[1], labels = function(z) mz(z)[[1L]],
        main = "Spectrum with labels",
        labelSrt = -30, labelPos = 2, add = TRUE)

})
