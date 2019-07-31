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

#' @title Join (map) peaks of two spectra
#'
#' @name joinPeaks
#'
#' @description
#'
#' These functions map peaks from two spectra with each other if the difference
#' between their m/z values is smaller than defined with parameters `tolerance`
#' and `ppm`. All functions take two matrices
#'
#' - `joinPeaks`: maps peaks from two spectra allowing to specify the type of
#'   *join* that should be performed: `type = "outer"` each peak in `x` will be
#'   matched with each peak in `y`, for peaks that do not match any peak in the
#'   other spectra an `NA` intensity is returned. With `type = "left"` all peaks
#'   from the left spectrum (`x`) will be matched with peaks in `y`. Peaks in
#'   `y` that do not match any peak in `x` are omitted. `type = "right"` is the
#'   same as `type = "left"` only for `y`. Only peaks that can be matched
#'   between `x` and `y` are returned by `type = "inner"`, i.e. only
#'   peaks present in both spectra are reported.
#'
#' @section Implementation notes:
#'
#' A mapping function must take two numeric matrices `x` and `y` as input and
#' must return `list` with two elements named `"x"` and `"y"` that represent
#' the aligned input matrices. The function should also have `...` in its
#' definition. Parameters `ppm` and `tolerance` are suggested but not required.
#'
#' @param ppm `numeric(1)` defining a relative, m/z-dependent, maximal accepted
#'     difference between m/z values of peaks from the two spectra to be
#'     matched/mapped.
#'
#' @param tolerance `numeric(1)` defining a constant maximal accepted difference
#'     between m/z values of peaks from the two spectra to be matched/mapped.
#'
#' @param type For `joinPeaks`: `character(1)` specifying the type of join that
#'     should be performed. See function description for details.
#'
#' @param x `matrix` with two columns `"mz"` and `"intensity"` containing the
#'     m/z and intensity values of the mass peaks of a spectrum.
#'
#' @param y `matrix` with two columns `"mz"` and `"intensity"` containing the
#'     m/z and intensity values of the mass peaks of a spectrum.
#'
#' @param ... option parameters.
#' @return
#'
#' All functions return a `list` of elements `"x"` and `"y"` each being a two
#' column matrix with m/z (first column) and intensity values (second column).
#' The two matrices contain the matched peaks between input matrices `x` and `y`
#' and hence have the same number of rows. Peaks present in `x` but not in the
#' `y` input matrix have m/z and intensity values of `NA` in the result matrix
#' for `y` (and *vice versa*).
#'
#' @author Johannes Rainer
#'
#' @importFrom MsCoreUtils join ppm
#'
#' @export
#'
#' @examples
#'
#' x <- cbind(c(31.34, 50.14, 60.3, 120.9, 230, 514.13, 874.1),
#'     1:7)
#' y <- cbind(c(12, 31.35, 70.3, 120.9 + ppm(120.9, 5),
#'     230 + ppm(230, 10), 315, 514.14, 901, 1202),
#'     1:9)
#'
#' ## No peaks with identical m/z
#' joinPeaks(x, y, ppm = 0, type = "inner")
#'
#' ## With ppm 10 two peaks are overlapping
#' joinPeaks(x, y, ppm = 10, type = "inner")
#'
#' ## Outer join: contain all peaks from x and y
#' joinPeaks(x, y, ppm = 10, type = "outer")
#'
#' ## Left join: keep all peaks from x and those from y that match
#' joinPeaks(x, y, ppm = 10, type = "left")
#'
#' ## Right join: keep all peaks from y and those from x that match. Using
#' ## a constant tolerance of 0.01
#' joinPeaks(x, y, tolerance = 0.01, type = "right")
joinPeaks <- function(x, y, type = "outer", tolerance = 0, ppm = 10, ...) {
    map <- join(x[, 1], y[, 1], type = type, tolerance = tolerance, ppm = ppm)
    list(x = x[map$x, , drop = FALSE], y = y[map$y, , drop = FALSE])
}
