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
    x[which(MsCoreUtils::between(x[, "intensity"], intensity)), , drop = FALSE]
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
#' Note that **all** peaks in `x` matching any of the m/z values in `mz` will
#' be considered.
#'
#' @inheritParams .peaks_remove
#'
#' @param mz `numeric` with the m/z values to match against. This needs to be
#'     a sorted `numeric` (**has to be checked upstream**).
#'
#' @param tolerance `numeric` with the tolerance. Can be of length 1 or equal
#'     length `mz`.
#'
#' @param ppm `numeric` with the ppm. Can be of length 1 or equal length `mz`.
#'
#' @param keep `logical(1)` whether the matching peaks should be retained
#'     (`keep = TRUE`, the default`) or dropped (`keep = FALSE`).
#'
#' @author Johannes Rainer
#'
#' @importFrom MsCoreUtils closest
#'
#' @noRd
.peaks_filter_mz_value <- function(x, spectrumMsLevel, mz = numeric(),
                                   tolerance = 0, ppm = 20,
                                   msLevel = spectrumMsLevel,
                                   keep = TRUE, ...) {
    if (!spectrumMsLevel %in% msLevel || !length(x))
        return(x)
    no_match <- is.na(MsCoreUtils::closest(x[, "mz"], mz, tolerance = tolerance,
                                           ppm = ppm, duplicates = "keep",
                                           .check = FALSE))
    if (keep) x[!no_match, , drop = FALSE]
    else x[no_match, , drop = FALSE]
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
                                   msLevel = spectrumMsLevel, keep = TRUE,
                                   ...) {
    if (!spectrumMsLevel %in% msLevel || !length(x))
        return(x)
    if (keep)
        x[MsCoreUtils::between(x[, "mz"], mz), , drop = FALSE]
    else x[!MsCoreUtils::between(x[, "mz"], mz), , drop = FALSE]
}

#' @description
#'
#' Bin peaks of a spectrum.
#'
#' @param mids mid points. This parameter is mandatory.
#'
#' @param zero.rm `logical`  indicating whether to remove bins with zero
#'     intensity. Defaults to `TRUE`, meaning the function will discard bins
#'     created with an intensity of 0 to enhance memory efficiency.
#'
#'
#' @inheritParams .peaks_remove
#'
#' @return `matrix` with columns `"mz"` and `"intensity"`
#'
#' @noRd
.peaks_bin <- function(x, spectrumMsLevel, binSize = 1L,
                       breaks = seq(floor(min(x[, 1])),
                                    ceiling(max(x[, 1])),
                                    by = binSize),
                       agg_fun = sum,
                       mids,
                       msLevel = spectrumMsLevel, zero.rm = TRUE, ...) {
    if (!(spectrumMsLevel %in% msLevel))
        return(x)
    bins <- MsCoreUtils::bin(x[, 2], x[, 1], size = binSize, breaks = breaks,
                             FUN = agg_fun, returnMids = FALSE)
    if (zero.rm) {
        keep <- which(bins != 0)
        cbind(mz = mids[keep], intensity = bins[keep])
    } else cbind(mz = mids, intensity = bins)
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
#' - `joinPeaks()`: maps peaks from two spectra allowing to specify the type of
#'   *join* that should be performed: `type = "outer"` each peak in `x` will be
#'   matched with each peak in `y`, for peaks that do not match any peak in the
#'   other spectra an `NA` intensity is returned. With `type = "left"` all peaks
#'   from the left spectrum (`x`) will be matched with peaks in `y`. Peaks in
#'   `y` that do not match any peak in `x` are omitted. `type = "right"` is the
#'   same as `type = "left"` only for `y`. Only peaks that can be matched
#'   between `x` and `y` are returned by `type = "inner"`, i.e. only
#'   peaks present in both spectra are reported.
#'
#' - `joinPeaksGnps()`: matches/maps peaks between spectra with the same
#'   approach used in GNPS: peaks are considered matching if a) the
#'   difference in their m/z values is smaller than defined by `tolerance`
#'   and `ppm` (this is the same as `joinPeaks`) **and** b) the difference of
#'   their m/z *adjusted* for the difference of the spectras' precursor is
#'   smaller than defined by `tolerance` and `ppm`. Based on this definition,
#'   peaks in `x` can match up to two peaks in `y` hence peaks in the returned
#'   matrices might be reported multiple times. Note that if one of
#'   `xPrecursorMz` or `yPrecursorMz` are `NA` or if both are the same, the
#'   results are the same as with [joinPeaks()]. To calculate GNPS similarity
#'   scores, [MsCoreUtils::gnps()] should be called on the aligned peak
#'   matrices (i.e. `compareSpectra` should be called with
#'   `MAPFUN = joinPeaksGnps` and `FUN = MsCoreUtils::gnps`).
#'
#' - `joinPeaksNone()`: does not perform any peak matching but simply returns
#'   the peak matrices in a `list`. This function should be used with the
#'   `MAPFUN` parameter of [compareSpectra()] if the spectra similarity
#'   function used (parameter `FUN` of `compareSpectra()`) performs its own
#'   peak matching and does hence not expect matched peak matrices as an input.
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
#' @param type For `joinPeaks()` and `joinPeaksGnps()`: `character(1)`
#'     specifying the type of join that should be performed. See function
#'     description for details.
#'
#' @param x `matrix` with two columns `"mz"` and `"intensity"` containing the
#'     m/z and intensity values of the mass peaks of a spectrum.
#'
#' @param xPrecursorMz for `joinPeaksGnps()`: `numeric(1)` with the precursor
#'     m/z of the spectrum `x`.
#'
#' @param y `matrix` with two columns `"mz"` and `"intensity"` containing the
#'     m/z and intensity values of the mass peaks of a spectrum.
#'
#' @param yPrecursorMz for `joinPeaksGnps()`: `numeric(1)` with the precursor
#'     m/z of the spectrum `y`.
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
#' @seealso
#'
#' - [compareSpectra()] for the function to calculate similarities between
#'   spectra.
#'
#' - [MsCoreUtils::gnps()] in the *MsCoreUtils* package for more information
#'   on the GNPS similarity score.
#'
#' @importFrom MsCoreUtils join ppm
#'
#' @export ppm
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
#' @importFrom MsCoreUtils join_gnps
#'
#' @rdname joinPeaks
joinPeaksGnps <- function(x, y, xPrecursorMz = NA_real_,
                          yPrecursorMz = NA_real_, tolerance = 0,
                          ppm = 0, type = "outer", ...) {
    map <- join_gnps(x[, 1L], y[, 1L], xPrecursorMz = xPrecursorMz,
                     yPrecursorMz = yPrecursorMz, tolerance = tolerance,
                     ppm = ppm, type = type, ...)
    list(x = x[map[[1L]], , drop = FALSE], y = y[map[[2L]], , drop = FALSE])
}

#' @export
#'
#' @rdname joinPeaks
joinPeaksNone <- function(x, y, ...) {
    list(x = x, y = y)
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

    n <- MsCoreUtils::noise(x[, 1L], x[, 2L], method = method, ...)

    l <- MsCoreUtils::localMaxima(x[, 2L], hws = halfWindowSize)

    p <- which(l & x[, 2L] > (snr * n))

    if (k > 0L) {
        cbind(mz = MsCoreUtils::refineCentroids(x = x[, 1L], y = x[, 2L], p = p,
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

#' @description
#'
#' Removes (Orbitrap) FFT artefacts (peaks) from a spectrum keeping peaks
#' potentially representing [13]C isotopes.
#'
#' @param x peaks matrix
#'
#' @param halfWindowSize `numeric(1)` defining the (m/z) window left and right
#'     of a peak to look for artefacts.
#'
#' @param threshold `numeric(1)` defining the relative intensity to the actual
#'     peak below which peaks are considered artefacts. Defaults to
#'     `threshold = 0.2` hence removing all peaks with an intensity below 0.2 *
#'     the peak intensity.
#'
#' @param keepIsotopes `logical(1)` whether it should be screened for potential
#'     isotope peaks. Defaults to `keepIsotopes = TRUE` thus candidate artefact
#'     peaks potentially being isotopes will **not** be removed.
#'
#' @param maxCharge `integer(1)`
#'
#' @param isotopeTolerance `numeric(1)`
#'
#' @author Jan Stanstrup
#'
#' @importFrom MsCoreUtils common
#'
#' @noRd
.peaks_remove_fft_artifact <- function(x, halfWindowSize = 0.05,
                                       threshold = 0.2,
                                       keepIsotopes = TRUE,
                                       maxCharge = 5,
                                       isotopeTolerance = 0.005, ...) {
    neutron   <- 1.0033548378 # really C12, C13 difference
    iso_dist  <- neutron / seq(from = 1, by = 1, to = maxCharge)
    ## just calculate isotopes that are in the halfWindowSize
    iso_dist <- iso_dist[iso_dist < halfWindowSize]
    find_isotopes <- keepIsotopes & length(iso_dist)
    if (find_isotopes)
        iso_dist <- c(-iso_dist, rev(iso_dist))

    mz <- x[, "mz"]
    int <- x[, "intensity"]

    ## left boundary
    lb <- findInterval(mz - halfWindowSize, mz) + 1L
    ## right boundary
    rb <- findInterval(mz + halfWindowSize, mz)

    ## region of interest (we just need to test if the window spans more than 1
    ## index)
    roi <- rb > lb

    ## test from the largest intensity
    idx <- seq_along(int)[roi][order(int[roi], decreasing = TRUE)]
    keep <- rep(TRUE, length(mz))

    for (i in idx) {
        if (keep[i]) {
            rem_candidate <- seq.int(lb[i], rb[i], by = 1L)
            rem_candidate <-
                rem_candidate[int[rem_candidate] / int[i] < threshold]

            if (find_isotopes) {
                cmm <- common(
                    mz[rem_candidate],
                    mz[i] + iso_dist,
                    tolerance = isotopeTolerance,
                    duplicates = "keep",
                    .check = FALSE
                )
                rem_candidate <- rem_candidate[!cmm]
            }

            keep[rem_candidate] <- FALSE
        }
    }
    x[keep, , drop = FALSE]
}

#' Function to keep only the monoisotopic peak for groups of (potential)
#' isotopologues. The peak with the lowest m/z is considered the monoisotopic
#' peak for each group of isotopologues. Isotope peaks are predicted using
#' the `MetaboCoreUtils::isotopologues` function. See tge respective
#' documentation for more information on parameters.
#'
#' @importFrom MetaboCoreUtils isotopologues isotopicSubstitutionMatrix
#'
#' @author Nir Shahaf, Johannes Rainer
#'
#' @param x peak `matrix` with m/z and intensity values.
#'
#' @noRd
.peaks_deisotope <-
    function(x, substDefinition = isotopicSubstitutionMatrix("HMDB_NEUTRAL"),
             tolerance = 0, ppm = 10, charge = 1, ...) {
        iso_grps <- MetaboCoreUtils::isotopologues(
                                         x, substDefinition = substDefinition,
                                         tolerance = tolerance, ppm = ppm,
                                         charge = charge)
        if (length(iso_grps)) {
            rem <- unique(unlist(lapply(iso_grps, `[`, -1), use.names = FALSE))
            x[-rem, , drop = FALSE]
        } else x
}

#' Function to *reduce* spectra keeping for each group of peaks with highly
#' similar m/z values (based on provided `tolerance` and `mz` the one with
#' the highest intensity. The *groups* of mass peaks are defined using the
#' `MsCoreUtils::group` function.
#'
#' @param x peaks `matrix` with m/z and instensity values.
#'
#' @author Nir Shahaf, Johannes Rainer
#'
#' @noRd
.peaks_reduce <- function(x, tolerance = 0, ppm = 10, ...) {
    grps <- group(x[, "mz"], tolerance = tolerance, ppm = ppm)
    l <- length(grps)
    if (l == grps[l])
        x
    else {
        il <- split(x[, "intensity"], as.factor(grps))
        idx <- vapply(il, which.max, integer(1), USE.NAMES = FALSE) +
            c(0, cumsum(lengths(il, use.names = FALSE))[-length(il)])
        x[idx, , drop = FALSE]
    }
}

.peaks_scale_intensities <- function(x, by = sum, spectrumMsLevel,
                                     msLevel, ...) {
    if (!spectrumMsLevel %in% msLevel || !nrow(x))
        return(x)
    ints <- x[, "intensity"]
    x[, "intensity"] <- ints / by(ints[!is.na(ints)])
    x
}

#' Combine peaks within each peak matrix if the difference of their m/z is
#' smaller than `tolerance` and `ppm`. Peak grouping is performed with the
#' `MsCoreUtils::group` function.
#'
#' Additional peak variables for combined peaks are dropped (replaced by `NA`)
#' if their values differ between the combined peaks.
#'
#' Note that `weighted = TRUE` overrides `mzFun` using the `weighted.mean`
#' to calculate the aggregated m/z values.
#'
#' @author Johannes Rainer
#'
#' @importFrom stats weighted.mean
#'
#' @noRd
.peaks_combine <- function(x, ppm = 20, tolerance = 0,
                           intensityFun = base::mean, mzFun = base::mean,
                           weighted = TRUE, spectrumMsLevel,
                           msLevel = spectrumMsLevel, ...) {
    if (!spectrumMsLevel %in% msLevel || !length(x))
        return(x)
    grps <- MsCoreUtils::group(x[, "mz"], tolerance = tolerance, ppm = ppm)
    lg <- length(grps)
    if (grps[lg] == lg)
        return(x)
    mzs <- split(x[, "mz"], grps)
    ints <- split(x[, "intensity"], grps)
    if (weighted)
        mzs <- unlist(mapply(
            mzs, ints, FUN = function(m, i) weighted.mean(m, i + 1),
            SIMPLIFY = FALSE, USE.NAMES = FALSE), use.names = FALSE)
    else mzs <- vapply1d(mzs, mzFun)
    ints <- vapply1d(ints, intensityFun)
    if (ncol(x) > 2) {
        lst <- lapply(x[!colnames(x) %in% c("mz", "intensity")], function(z) {
            z <- lapply(split(z, grps), unique)
            z[lengths(z) > 1] <- NA
            unlist(z, use.names = FALSE, recursive = FALSE)
        })
        do.call(cbind.data.frame,
                c(list(mz = mzs, intensity = ints), lst))
    } else
        cbind(mz = mzs, intensity = ints)
}

#' Keeps all peaks except those matching the precursor m/z with `tolerance`
#' and `ppm`.
#'
#' @author Johannes Rainer
#'
#' @noRd
.peaks_filter_precursor_ne <- function(x, ppm = 20, tolerance = 0,
                                       precursorMz, spectrumMsLevel,
                                       msLevel = spectrumMsLevel, ...) {
    if (!spectrumMsLevel %in% msLevel || !nrow(x))
        return(x)
    keep <- is.na(MsCoreUtils::closest(x[, "mz"], precursorMz, ppm = ppm,
                          tolerance = tolerance, duplicates = "keep",
                          .check = FALSE))
    x[keep, , drop = FALSE]
}

#' Keeps peaks with an m/z smaller than the precursor m/z. `tolerance` and
#' `ppm` are subtracted from the precursor m/z to ensure that also the
#' precursor m/z is removed.
#'
#' @author Johannes Rainer
#'
#' @noRd
.peaks_filter_precursor_keep_below <- function(x, ppm = 20, tolerance = 0,
                                               precursorMz, spectrumMsLevel,
                                               msLevel = spectrumMsLevel, ...) {
    if (!spectrumMsLevel %in% msLevel || !nrow(x))
        return(x)
    pmz <- precursorMz - tolerance - ppm(precursorMz, ppm = ppm)
    x[x[, "mz"] < pmz, , drop = FALSE]
}

#' filter a peak matrix `x` by (arbitrary) numeric ranges for spectra and/or
#' peaks variables. ranges for spectra and peaks variables are combined using
#' a logical AND, rows in the provided range matrices with a logical OR.
#'
#' Used by `filterPeaksRanges()` function for `Spectra`.
#'
#' @param svars `character` with the spectra variables for which filter ranges
#'     where provided.
#'
#' @param pvars `character` with the peaks variables for which filter ranges
#'     where provided.
#'
#' @param ranges `list` with `numeric` two-column matrices with the
#'     user-provided ranges. The number of rows of all matrices is expected
#'     to match.
#'
#' @param spectrumMsLevel `integer(1)` with the MS level of the peak matrix'
#'     spectrum.
#'
#' @param keep `logical(1)` whether mass peaks that match the filters should be
#'     kept or removed.
#'
#' @param ... values for all spectra variables defined in `svars` are expected
#'     to be passed through `...` as `name = value` pairs.
#'
#' @author Johannes Rainer
#'
#' @noRd
.peaks_filter_ranges <- function(x, svars = character(),
                                 pvars = character(),
                                 ranges, spectrumMsLevel,
                                 keep = TRUE, ...) {
    svalue <- list(..., msLevel = spectrumMsLevel)
    nx <- nrow(x)
    sel <- rep(FALSE, nx)
    for (i in seq_len(nrow(ranges[[1L]]))) {
        ## check ranges for spectra variables
        svars_ok <- vapply(svars, function(z)
            MsCoreUtils::between(svalue[[z]], ranges[[z]][i, ]), TRUE,
            USE.NAMES = FALSE)
        if (!anyNA(svars_ok) && all(svars_ok)) {
            if (length(pvars)) {
                ## check ranges for peaks variables
                tmp <- rowSums(do.call(cbind, lapply(pvars, function(z) {
                    MsCoreUtils::between(x[, z], ranges[[z]][i, ])
                }))) == length(pvars)
                tmp[is.na(tmp)] <- FALSE
                sel <- sel | tmp
            } else {
                ## No need to check further, because we have a match
                if (keep) return(x)
                else return(x[logical(), , drop = FALSE])
            }
        }
    }
    if (keep) x[sel, , drop = FALSE]
    else x[!sel, , drop = FALSE]
}

#' Check for presence of peaks defined by their m/z value. Note that this
#' function does **not** return a peak matrix, but only a logical of length 1!
#'
#' @return `logical(1)`
#' @noRd
.peaks_contain_mz <- function(x, mz = numeric(), tolerance = 0, ppm = 20,
                              condFun = any, ...) {
    condFun(common(mz, x[, "mz"], tolerance = tolerance, ppm = ppm))
}
