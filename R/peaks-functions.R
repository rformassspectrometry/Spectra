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
#' Replace peak intensities if intensity below `threshold` with `value`.
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
.peaks_replace_intensity <- function(x, spectrumMsLevel, centroided = NA,
                                     threshold = min, value = 0,
                                     msLevel = spectrumMsLevel, ...) {
    if (!spectrumMsLevel %in% msLevel || !length(x))
        return(x)
    if (is.na(centroided)) {
        warning("Centroided undefined (NA): keeping spectrum as is.")
        return(x)
    }
    if (is.function(threshold))
        threshold <- threshold(x[, "intensity"], na.rm = TRUE)
    if (centroided) {
        x[x[, "intensity"] <= threshold, "intensity"] <- value
    } else {
        ints <- x[, "intensity"]
        peakRanges <- as(ints > 0L, "IRanges")
        toLow <- max(extractList(ints, peakRanges)) <= threshold
        x[, "intensity"] <- replaceROWS(ints, peakRanges[toLow], value)
    }
    x
}

#' @description
#'
#' Filtering the spectrum keeping only peaks which are within the provided
#' intensity range. Note that also peaks with `NA` intensities are removed.
#'
#' @inheritParams .peaks_remove
#'
#' @param intensity `numeric(2)` with the lower and upper range.
#'
#' @importFrom MsCoreUtils between
#'
#' @noRd
.peaks_filter_intensity <- function(x, spectrumMsLevel, intensity = c(0, Inf),
                                    msLevel = spectrumMsLevel, ...) {
    if (!spectrumMsLevel %in% msLevel || !length(x))
        return(x)
    x[which(between(x[, "intensity"], intensity)), , drop = FALSE]
}

#' @description
#'
#' Filter peaks with an intensity-based function.
#'
#' @inheritParams .peaks_remove
#'
#' @param intfun function which takes intensities as first parameter and
#'     returns a `logical` of length equal to the number of peaks in the
#'     spectrum.
#'
#' @author Johannes Rainer
#'
#' @noRd
.peaks_filter_intensity_function <- function(x, spectrumMsLevel, intfun,
                                             msLevel = spectrumMsLevel,
                                             args = list(), ...) {
    if (!spectrumMsLevel %in% msLevel || !length(x))
        return(x)
    keep <- do.call(intfun, args = c(list(x[, "intensity"]), args))
    if (!is.logical(keep) || length(keep) != nrow(x))
        stop("Error in filterIntensity: the provided function does not return ",
             "a logical vector of length equal to the number of peaks.",
             call. = FALSE)
    x[which(keep), , drop = FALSE]
}

#' @description
#'
#' Filter the spectrum keeping only peaks that match the provided m/z value(s).
#'
#' @inheritParams .peaks_remove
#'
#' @param mz `numeric` with the m/z values to match against. This needs to be
#'     a sorted `numeric` (has to be checked upstream).
#'
#' @param tolerance `numeric` with the tolerance. Can be of length 1 or equal
#'     length `mz`.
#'
#' @param ppm `numeric` with the ppm. Can be of length 1 or equal length `mz`.
#'
#' @author Johannes Rainer
#'
#' @importFrom MsCoreUtils closest
#'
#' @noRd
.peaks_filter_mz_value <- function(x, spectrumMsLevel, mz = numeric(),
                                   tolerance = 0, ppm = 20,
                                   msLevel = spectrumMsLevel, ...) {
    if (!spectrumMsLevel %in% msLevel || !length(x))
        return(x)
    idx <- closest(mz, x[, "mz"], tolerance = tolerance, ppm = ppm,
                   duplicates = "closest", .check = FALSE)
    x[idx[!is.na(idx)], , drop = FALSE]
}

#' @description
#'
#' Filter the spectrum keeping only peaks that are within the provided m/z
#' range.
#'
#' @inheritParams .peaks_remove
#'
#' @param mz `numeric(2)` with the lower and upper m/z value.
#'
#' @author Johannes Rainer
#'
#' @noRd
.peaks_filter_mz_range <- function(x, spectrumMsLevel, mz = numeric(),
                                   msLevel = spectrumMsLevel, ...) {
    if (!spectrumMsLevel %in% msLevel || !length(x))
        return(x)
    x[between(x[, "mz"], mz), , drop = FALSE]
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
#' - `joinPeaksGnps` matches/maps peaks between spectra with an approach
#'   similar to the one used in GNPS: peaks are considered matching if a) the
#'   difference in their m/z values is smaller than defined by `tolerance`
#'   and `ppm` (this is the same as `joinPeaks`) **and** b) the difference of
#'   their m/z *adjusted* for the difference of the spectras' precursor is
#'   smaller than defined by `tolerance` and `ppm`. Based on this definition,
#'   peaks in `x` can match up to two peaks in `y` hence peaks in the returned
#'   matrices might be reported multiple times. Note that if one of
#'   `xPrecursorMz` or `yPrecursorMz` are `NA` or if both are the same, the
#'   results are the same as with [joinPeaks()]. To calculate GNPS similarity
#'   scores, [gnps()] should be called on the aligned peak matrices (i.e.
#'   `compareSpectra` should be called with `MAPFUN = joinPeaksGnps` and
#'   `FUN = gnps`).
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
#' @param type For `joinPeaks` and `joinPeaksGnps`: `character(1)` specifying
#'     the type of join that should be performed. See function description for
#'     details.
#'
#' @param x `matrix` with two columns `"mz"` and `"intensity"` containing the
#'     m/z and intensity values of the mass peaks of a spectrum.
#'
#' @param xPrecursorMz for `joinPeaksGnps`: `numeric(1)` with the precursor m/z
#'     of the spectrum `x`.
#'
#' @param y `matrix` with two columns `"mz"` and `"intensity"` containing the
#'     m/z and intensity values of the mass peaks of a spectrum.
#'
#' @param yPrecursorMz for `joinPeaksGnps`: `numeric(1)` with the precursor m/z
#'     of the spectrum `y`.
#'
#' @param ... optional parameters passed to the [MsCoreUtils::join()] function.
#'
#' @return
#'
#' All functions return a `list` of elements `"x"` and `"y"` each being a two
#' column matrix with m/z (first column) and intensity values (second column).
#' The two matrices contain the matched peaks between input matrices `x` and `y`
#' and hence have the same number of rows. Peaks present in `x` but not in the
#' `y` input matrix have m/z and intensity values of `NA` in the result matrix
#' for `y` (and *vice versa*).
#'
#' @author Johannes Rainer, Michael Witting
#'
#' @seealso [gnps()]
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
#'
#' ## GNPS-like peak matching
#'
#' ## Define spectra
#' x <- cbind(mz = c(10, 36, 63, 91, 93), intensity = c(14, 15, 999, 650, 1))
#' y <- cbind(mz = c(10, 12, 50, 63, 105), intensity = c(35, 5, 16, 999, 450))
#' ## The precursor m/z
#' pmz_x <- 91
#' pmz_y <- 105
#'
#' ## Plain joinPeaks identifies only 2 matching peaks: 1 and 5
#' joinPeaks(x, y)
#'
#' ## joinPeaksGnps finds 4 matches
#' joinPeaksGnps(x, y, pmz_x, pmz_y)
#'
#' ## with one of the two precursor m/z being NA, the result are the same as
#' ## with joinPeaks (with type = "left").
#' joinPeaksGnps(x, y, pmz_x, yPrecursorMz = NA)
joinPeaks <- function(x, y, type = "outer", tolerance = 0, ppm = 10, ...) {
    map <- join(x[, 1], y[, 1], type = type, tolerance = tolerance,
                ppm = ppm, ...)
    list(x = x[map$x, , drop = FALSE], y = y[map$y, , drop = FALSE])
}

#' @export
#'
#' @rdname joinPeaks
joinPeaksGnps <- function(x, y, xPrecursorMz = NA_real_,
                          yPrecursorMz = NA_real_, tolerance = 0,
                          ppm = 0, type = "outer", ...) {
    pdiff <- yPrecursorMz - xPrecursorMz
    map <- join(x[, 1L], y[, 1L], tolerance = tolerance, ppm = ppm,
                type = type, ...)
    if (is.finite(pdiff) && pdiff != 0) {
        pmap <- join(x[, 1L] + pdiff, y[, 1L], tolerance = tolerance,
                     ppm = ppm, type = type, ...)
        ## Keep only matches here
        nona <- !(is.na(pmap[[1L]]) | is.na(pmap[[2L]]))
        if (any(nona)) {
            map[[1L]] <- c(map[[1L]], pmap[[1L]][nona])
            map[[2L]] <- c(map[[2L]], pmap[[2L]][nona])
            idx <- order(map[[1L]])
            map[[1L]] <- map[[1L]][idx]
            map[[2L]] <- map[[2L]][idx]
        }
    }
    list(x = x[map[[1L]], , drop = FALSE], y = y[map[[2L]], , drop = FALSE])
}

#' @title GNPS spectra similarity score
#'
#' @description
#'
#' The GNPS approach considers in its similarity score, in addition to directly
#' matching peaks between two spectra, also matches of peaks which difference in
#' m/z values matches the difference of the spectra's precursor m/z. For peaks
#' that match multiple peaks in the other spectrum the matching peak pair with
#' the higher value/similarity is considered in the final similarity score.
#'
#' To calculate GNPS scores, the [joinPeaksGnps()] function should be used first
#' to match the peaks between the compared spectra and `gnps` should then be
#' called on the resulting aligned peak matrices. To use the GNPS score on a
#' [Spectra()] object `compareSpectra` should be used with
#' `MAPFUN = joinPeaksGnps` and `FUN = gnps`.
#'
#' @details
#'
#' The implementation bases on the R code from the publication listed in the
#' references.
#'
#' @param x `numeric` peak `matrix` with m/z (first column) and intensity
#'     values (second column). Rows in `x` and `y` are considered to contain
#'     matching peaks, i.e. the first row in `x` contains a peak matching the
#'     first peak (row) of `y` or `NA_real_` if the peak does not match.
#'
#' @param y `numeric` peak `matrix` with m/z (first column) and intensity
#'     values (second column). Rows in `x` and `y` are considered to contain
#'     matching peaks, i.e. the first row in `x` contains a peak matching the
#'     first peak (row) of `y` or `NA_real_` if the peak does not match.
#'
#' @param ... ignored
#'
#' @author Johannes Rainer, Michael Witting, based on the code from the ref
#'
#' @importFrom clue solve_LSAP
#'
#' @importFrom stats complete.cases
#'
#' @references
#'
#' Xing S, Hu Y, Yin Z, Liu M, Tang X, Fang M, Huan T. Retrieving and Utilizing
#' Hypothetical Neutral Losses from Tandem Mass Spectra for Spectral Similarity
#' Analysis and Unknown Metabolite Annotation. *Anal Chem.*
#' 2020 Nov 3;92(21):14476-14483. doi: 10.1021/acs.analchem.0c02521.
#'
#' @export
#'
#' @examples
#'
#' ## Define spectra
#' x <- cbind(mz = c(10, 36, 63, 91, 93), intensity = c(14, 15, 999, 650, 1))
#' y <- cbind(mz = c(10, 12, 50, 63, 105), intensity = c(35, 5, 16, 999, 450))
#' ## The precursor m/z
#' pmz_x <- 91
#' pmz_y <- 105
#'
#' ## joinPeaksGnps finds 4 matches
#' map <- joinPeaksGnps(x, y, pmz_x, pmz_y)
#' map
#'
#' ## Calculate similarity
#' gnps(map$x, map$y)
gnps <- function(x, y, ...) {
    if (nrow(x) != nrow(y))
        stop("'x' and 'y' are expected to be aligned peak matrices (i.e. ",
             "having the same number of rows).")
    ## Scale intensities; !duplicated because we can have duplicated matches.
    x_sum <- sum(x[!duplicated(x[, 1]), 2], na.rm = TRUE)
    y_sum <- sum(y[!duplicated(y[, 1]), 2], na.rm = TRUE)
    ## is 0 if only NAs in input - avoids division through 0
    if (x_sum == 0 || y_sum == 0)
        return(0)
    ## Keep only matches.
    keep <- which(complete.cases(cbind(x[, 1], y[, 1])))
    l <- length(keep)
    if (!l)
        return(0)
    x <- x[keep, , drop = FALSE]
    y <- y[keep, , drop = FALSE]
    scores <- sqrt(x[, 2]) / sqrt(x_sum) * sqrt(y[, 2]) / sqrt(y_sum)

    x_idx <- as.integer(factor(x[, 1]))
    y_idx <- as.integer(factor(y[, 1]))
    score_mat <- matrix(0, nrow = l, ncol = l)
    seq_l <- seq_len(l)
    for (i in seq_l) {
        score_mat[x_idx[i], y_idx[i]] <- scores[i]
    }
    best <- solve_LSAP(score_mat, maximum = TRUE)
    sum(score_mat[cbind(seq_l, as.integer(best))], na.rm = TRUE)
}


#' @importFrom MsCoreUtils localMaxima noise refineCentroids
#'
#' @description
#'
#' Simple peak detection based on local maxima above snr * noise.
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
                        snr = 0L, k = 0L, descending = FALSE, threshold = 0,
                        msLevel = spectrumMsLevel, ...) {
    if (!(spectrumMsLevel %in% msLevel))
        return(x)
    if (!nrow(x)) {
        warning("Spectrum is empty. Nothing to pick.")
        return(x)
    }

    n <- noise(x[, 1L], x[, 2L], method = method, ...)

    l <- localMaxima(x[, 2L], hws = halfWindowSize)

    p <- which(l & x[, 2L] > (snr * n))

    if (k > 0L) {
        cbind(mz = refineCentroids(x = x[, 1L], y = x[, 2L], p = p,
                                   k = k, threshold = threshold,
                                   descending = descending),
              intensity = x[p, 2L])
    } else {
        x[p, , drop = FALSE]
    }
}

#' @description
#'
#' Intensity smoothing.
#'
#' @inheritParams .peaks_remove
#' @param coef `matrix` smoothing coefficients (generated by
#' [MsCoreUtils::coefMA()], [MsCoreUtils::coefWMA()], or [MsCoreUtils::coefSG].
#'
#' @return `matrix` with columns `"mz"` and `"intensity"`.
#'
#' @author Sebastian Gibb
#'
#' @noRd
.peaks_smooth <- function(x, spectrumMsLevel, msLevel = spectrumMsLevel,
                          coef, ...) {
    if (!(spectrumMsLevel %in% msLevel))
        return(x)
    if (!nrow(x)) {
        warning("Spectrum is empty. Nothing to smooth.")
        return(x)
    }
    x[, 2L] <- MsCoreUtils::smooth(x[, 2L], cf = coef)
    x
}
