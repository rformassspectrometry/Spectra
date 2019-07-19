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
#' Bin peaks of a spectrum.
#'
#' @inheritParams .peaks_remove
#'
#' @importFrom MsCoreUtils bin
#'
#' @return `matrix` with columns `"mz"` and `"intensity"`
#'
#' @noRd
.peaks_bin <- function(x, spectrumMsLevel, binSize = 1L,
                       breaks = seq(floor(min(x[, 1])),
                                    ceiling(max(x[, 1])),
                                    by = binSize), FUN = sum,
                       msLevel = spectrumMsLevel, ...) {
    if (!(spectrumMsLevel %in% msLevel))
        return(x)
    bins <- bin(x[, 2], x[, 1], size = binSize, breaks = breaks,
                FUN = FUN)
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

#' @description Compare peaks between 2 spectra
#'
#' @description
#'
#' Compare peaks from one spectrum to peaks from another spectrum with the
#' provided function (`FUN`). `.peaks_compare_intensities` ensures that peaks
#' are first matched between the two spectra with the [closest()] function.
#'
#' @param x peaks `matrix` with columns `"mz"` and `"intensity"`.
#'
#' @param y peaks `matrix` with columns `"mz"` and `"intensity"`.
#'
#' @param FUN `function` to be applied to the intensity values of matched
#'     peaks between peak lists `x` and `y`. The function should take two
#'     `numeric` vectors as their first argument and should return a single
#'     `numeric`. Parameter `...` will be passed to the function.
#'
#' @param tolerance `numeric(1)` with the acceptable (constant) difference in
#'     m/z values for matching peaks.
#'
#' @param ppm `numeric(1)` specifying a peak-specific, relative,
#'     tolerance for the peak matching. Peaks between the two spectra with a
#'     difference in their m/z smaller than *ppm* of the m/z of the peak (of
#'     the second spectrum) are matched. See also [ppm()].
#'
#' @param ... additional parameters passed to the `FUN` function.
#'
#' @return The result from `FUN`, which should be a `numeric(1)`.
#'
#' @author Johannes Rainer
#'
#' @noRd
#'
#' @importFrom stats cor
#'
#' @importFrom MsCoreUtils closest ppm
#'
#' @examples
#'
#' x <- cbind(c(31.34, 34.43, 54.65, 300.12, 514.5, 856.12),
#'      c(13, 15, 56, 7, 45, 321))
#'
#' y <- cbind(c(34.431, 35.34, 50.24, 300.121, 514.51, 856.122),
#'      c(123, 34, 323, 432, 12, 433))
#'
#' Spectra:::.peaks_compare_intensities(x, y)
#' Spectra:::.peaks_compare_intensities(x, y, ppm = 40)
#' Spectra:::.peaks_compare_intensities(x, y, ppm = 0)
#'
#' ## To get the number of common peaks
#' Spectra:::.peaks_compare_intensities(x, y, FUN = function(x, y) length(x))
#' Spectra:::.peaks_compare_intensities(x, y, FUN = function(x, y) length(x), ppm = 40)
#' ## NA if there are none common
#' Spectra:::.peaks_compare_intensities(x, y, FUN = function(x, y) length(x), ppm = 0)
.peaks_compare_intensities <- function(x, y, FUN = cor, tolerance = 0,
                                       ppm = 20, ...) {
    matches <- closest(x[, 1], y[, 1], tolerance = tolerance + ppm(y[, 1], ppm),
                       duplicates = "closest")
    not_na <- !is.na(matches)
    if (any(not_na)) {
        FUN(x[not_na, 2], y[matches[not_na], 2], ...)
    } else NA_real_
}
