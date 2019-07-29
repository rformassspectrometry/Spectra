#' @title Peak data functions
#'
#' @description
#'
#' *Peaks data* functions are functions to work with peak matrices, two column
#' matrices with m/z and intensity values for peaks of a spectrum (one matrix
#' containing the peak data of one spectrum). The functions are supposed to
#' provide data manipulation methods to peak data, taking a (numeric) peak
#' matrix as input and returning a (numeric) peak matrix.
#'
#' @section Implementation notes:
#'
#' - The name of a peak data function should start with `.peaks_`.
#' - The function needs to have the `...` parameters to allow passing additional
#'   arguments and don't throw an error if e.g. `msLevel`, which will be passed
#'   to any `.peaks_` function, is not used by the function.
#'
#' For an example implementation see `.peaks_remove` below.
#'
#' @noRd
NULL

#' @description
#'
#' Remove peaks from spectrum data if intensity below `t`.
#'
#' @note
#'
#' This function takes base R types as input (`matrix`, `logical(1)` etc). It
#' might eventually be better if we used a single parameter being a `DataFrame`?
#'
#' @param x `matrix` with spectrum data (columns `"mz"` and `"intensity"`).
#'
#' @param spectrumMsLevel `integer(1)` defining the MS level of the spectrum.
#'
#' @param centroided `logical(1)` defining whether the spectrum data is
#'     centroided
#'
#' @param msLevel optional `integer` defining the MS level(s) to which the
#'     function should be applied.
#'
#' @param ... required. Support for passing additional parameters without
#'     throwing an error.
#'
#' @return `matrix` with columns `"mz"` and `"intensity"`.
#'
#' @importMethodsFrom IRanges extractList replaceROWS
#'
#' @importClassesFrom IRanges IRanges
#'
#' @noRd
.peaks_remove <- function(x, spectrumMsLevel, centroided = NA, t = "min",
                          msLevel = spectrumMsLevel, ...) {
    if (!spectrumMsLevel %in% msLevel || !length(x))
        return(x)
    if (is.na(centroided)) {
        warning("Centroided undefined (NA): keeping spectrum as is.")
        return(x)
    }
    if (t == "min")
        t <- min(x[x[, "intensity"] > 0, "intensity"])
    if (!is.numeric(t))
        stop("'t' must either be 'min' or numeric.")
    if (centroided) {
        x[x[, "intensity"] <= t, "intensity"] <- 0
    } else {
        ints <- x[, "intensity"]
        peakRanges <- as(ints > 0L, "IRanges")
        toLow <- max(extractList(ints, peakRanges)) <= t
        x[, "intensity"] <- replaceROWS(ints, peakRanges[toLow], 0)
    }
    x
}

#' @description
#'
#' Clean spectrum by removing 0-intensity peaks.
#'
#' @inheritParams .peaks_remove
#'
#' @return `matrix` with columns `"mz"` and `"intensity"`.
#'
#' @noRd
.peaks_clean <- function(x, spectrumMsLevel, all = FALSE,
                         msLevel = spectrumMsLevel, ...) {
    if (!spectrumMsLevel %in% msLevel || !length(x))
        return(x)
    x[utils.clean(x[, "intensity"], all), , drop = FALSE]
}

#' @description
#'
#' Bin spectrum. Code taken from MSnbase:::bin_Spectrum
#'
#' @inheritParams .peaks_remove
#'
#' @return `matrix` with columns `"mz"` and `"intensity"`
#'
#' @noRd
.peaks_bin <- function(x, spectrumMsLevel, binSize = 1L,
                       breaks = seq(floor(min(x[, 1])),
                                    ceiling(max(x[, 1])),
                                    by = binSize), fun = sum,
                       msLevel = spectrumMsLevel, ...) {
    if (!(spectrumMsLevel %in% msLevel))
        return(x)
    bins <- .bin_values(x[, 2], x[, 1], binSize = binSize, breaks = breaks,
                        fun = fun)
    cbind(mz = bins$mids, intensity = bins$x)
}

#' @importFrom stats quantile
#'
#' @description
#'
#' Simple function to estimate whether the spectrum is centroided. Was formerly
#' called `isCentroided`.
#'
#' @return `logical(1)`
#'
#' @noRd
.peaks_is_centroided <- function(x, spectrumMsLevel, centroided = NA,
                                 k = 0.025, qtl = 0.9) {
    .qtl <- quantile(x[, 2], qtl)
    x <- x[x[, 2] > .qtl, 1]
    quantile(diff(x), 0.25) > k
}

#' @importFrom MsCoreUtils localMaxima noise refineCentroids
#'
#' @description
#'
#' Simple peak detection based on local maxima above SNR * noise.
#'
#' @inheritParams .peaks_remove
#' @param `integer(1)`, half window size, the resulting window reaches from
#' `(i - halfWindowSize):(i + halfWindowSize)`.
#' @param k `integer(1)`, similar to `halfWindowSize`, number of values left
#'  and right of the peak that should be considered in the weighted mean
#'  calculation. If zero no refinement is done.
#' @param threshold `double(1)`, proportion of the maximal peak intensity.
#'  Just values above are used for the weighted mean calclulation.
#' @param descending `logical`, if `TRUE` just values between the nearest
#'  valleys around the peak centroids are used.
#'
#' @return `matrix` with columns `"mz"` and `"intensity"`.
#'
#' @noRd
.peaks_pick <- function(x, spectrumMsLevel, centroided = NA,
                        halfWindowSize = 2L, method = c("MAD", "SuperSmoother"),
                        SNR = 0L, k = 0L, descending = FALSE, threshold = 0,
                        msLevel = spectrumMsLevel, ...) {
    if (!(spectrumMsLevel %in% msLevel) || isTRUE(centroided))
        return(x)
    if (!nrow(x)) {
        warning("Spectrum is empty. Nothing to pick.")
        return(x)
    }

    n <- noise(x[, 1L], x[, 2L], method = method, ...)

    l <- localMaxima(x[, 2L], hws = halfWindowSize)

    p <- which(l & x[, 2L] > (SNR * n))

    if (k > 0L) {
        cbind(mz = refineCentroids(x = x[, 1L], y = x[, 2L], p = p,
                                   k = k, threshold = threshold,
                                   descending = descending),
              intensity = x[p, 2L])
    } else {
        x[p, , drop = FALSE]
    }
}
