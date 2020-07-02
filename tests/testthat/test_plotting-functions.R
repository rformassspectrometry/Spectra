#' @examples
#' # each base::graphics plot function must be wrapped by an anonymous function
#' # that could be called by `vdiffr::expect_doppelganger()`
#' # use `manage_cases()` to add new/verify changed plots
#' vdiffr::manage_cases(filter = "plot")

ints <- c(4.3412, 12, 8, 34, 23.4)
mzs <- c(13.453421, 43.433122, 46.6653553, 129.111212, 322.24432)
df <- DataFrame(msLevel = 1L, rtime = 123.12)
df$mz <- list(mzs)
df$intensity <- list(ints)
s <- Spectra(df)

test_that(".plot_single_spectrum works", {
    ## Just checking here that we don't get failures. vdiffr doesn't work
    ## because we're producing multiple pages.
    vdiffr::expect_doppelganger(
                "plot_single_spectrum-basic",
                function().plot_single_spectrum(s[1])
            )
    vdiffr::expect_doppelganger(
                "plot_single_spectrum-xlim",
                function().plot_single_spectrum(s[1], xlim = c(0, 700),
                                                ylim = c(0, 400))
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
