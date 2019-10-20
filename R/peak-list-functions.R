#' @title Function dealing with peak lists
#'
#' @description
#'
#' These functions are supposed to take a `list` of peak matrices (i.e. two
#' column matrices with columns `"mz"` and `"intensity"` as input and return
#' a single `matrix` combining the individual input peak matrices.
#'
#' @noRd
NULL

#' @title Combine peaks with similar m/z across spectra
#'
#' @description
#'
#' `combinePeaks` aggregates provided peak matrices into a single peak matrix.
#' Peaks are grouped by their m/z values with the `group()` function from the
#' `MsCoreUtils` package. All peaks in all provided spectra are first ordered
#' by their m/z and consecutively grouped into one group if the (pairwise)
#' difference between them is smaller than specified with parameter `tolerance`
#' and `ppm` (see [group()] for grouping details and examples).
#'
#' The m/z and intensity values for the resulting peak matrix are calculated
#' using the `mzFun` and `intensityFun` on the grouped m/z and intensity values.
#'
#' Parameters `main` and `unionPeaks` deal with special cases if e.g. only
#' peaks should be reported that are present in one of the provided peak
#' matrices. If for example `main = 2L` and `unionPeaks = FALSE`, only peaks
#' are returned wich are present in the second peak matrix.
#'
#' Setting `timeDomain` to `TRUE` causes grouping to be performed on the square
#' root of the m/z values.
#'
#' @details
#'
#' For general merging of spectra, the `tolerance` and/or `ppm` should be
#' manually specified based on the precision of the MS instrument. Peaks
#' from spectra with a difference in their m/z being smaller than `tolerance`
#' or smaller than `ppm` of their m/z are grouped into the same final peak.
#'
#'
#' Some details for the combination of consecutive spectra of an LC-MS run:
#'
#' The m/z values of the same ion in consecutive scans (spectra) of a LC-MS run
#' will not be identical. Assuming that this random variation is much smaller
#' than the resolution of the MS instrument (i.e. the difference between
#' m/z values within each single spectrum), m/z value groups are defined
#' across the spectra and those containing m/z values of the `main` spectrum
#' are retained.
#' Intensities and m/z values falling within each of these m/z groups are
#' aggregated using the `intensityFun` and `mzFun`, respectively. It is
#' highly likely that all QTOF profile data is collected with a timing circuit
#' that collects data points with regular intervals of time that are then later
#' converted into m/z values based on the relationship `t = k * sqrt(m/z)`. The
#' m/z scale is thus non-linear and the m/z scattering (which is in fact caused
#' by small variations in the time circuit) will thus be different in the lower
#' and upper m/z scale. m/z-intensity pairs from consecutive scans to be
#' combined are therefore defined by default on the square root of the m/z
#' values. With `timeDomain = FALSE`, the actual m/z values will be used.
#'
#' @param x `list` of peak matrices.
#'
#' @param intensityFun `function` to be used to combine intensity values for
#'     matching peaks. By default the mean intensity value is returned.
#'
#' @param mzFun `function` to be used to combine m/z values for matching peaks.
#'     By default the mean m/z value is returned.
#'
#' @param weighted `logical(1)` defining whether m/z values for matching peaks
#'     should be calculated by an intensity-weighted average of the individuak
#'     m/z values. This overrides parameter `mzFun`.
#'
#' @param tolerance `numeric(1)`
#'
#' @param ppm `numeric(1)`
#'
#' @param timeDomain `logical(1)`
#'
#' @param main `integer(1)`
#'
#' @param unionPeaks `logical(1)`
#'
#' @param ... additional parameters to the `mzFun` and `intensityFun` functions.
#'
#' @family peak matrix combining functions
#'
#' @return
#'
#' Peaks `matrix` with m/z and intensity values representing the aggregated
#' values across the provided peak matrices.
#'
#' @author Johannes Rainer
#'
#' @importFrom MsCoreUtils group
#'
#' @importFrom stats weighted.mean
#'
#' @noRd
#'
#' @examples
#'
#' set.seed(123)
#' mzs <- seq(1, 20, 0.1)
#' ints1 <- abs(rnorm(length(mzs), 10))
#' ints1[11:20] <- c(15, 30, 90, 200, 500, 300, 100, 70, 40, 20) # add peak
#' ints2 <- abs(rnorm(length(mzs), 10))
#' ints2[11:20] <- c(15, 30, 60, 120, 300, 200, 90, 60, 30, 23)
#' ints3 <- abs(rnorm(length(mzs), 10))
#' ints3[11:20] <- c(13, 20, 50, 100, 200, 100, 80, 40, 30, 20)
#'
#' ## Create the peaks matrices
#' p1 <- cbind(mz = mzs + rnorm(length(mzs), sd = 0.01),
#'     intensity = ints1)
#' p2 <- cbind(mz = mzs + rnorm(length(mzs), sd = 0.01),
#'     intensity = ints2)
#' p3 <- cbind(mz = mzs + rnorm(length(mzs), sd = 0.009),
#'     intensity = ints3)
#'
#' ## Combine the spectra. With `tolerance = 0` and `ppm = 0` only peaks with
#' ## **identical** m/z are combined
#' p <- combinePeaks(list(p1, p2, p3))
#'
#' ## Plot the spectra before and after combining
#' par(mfrow = c(2, 1), mar = c(4.3, 4, 1, 1))
#' plot(p1[, 1], p1[, 2], xlim = range(mzs[5:25]), type = "h", col = "red")
#' points(p2[, 1], p2[, 2], type = "h", col = "green")
#' points(p3[, 1], p3[, 2], type = "h", col = "blue")
#'
#' plot(p[, 1], p[, 2], xlim = range(mzs[5:25]), type = "h",
#'     col = "black")
#'
#' ## Combine spectra with `tolerance = 0.05`. This will merge all triplets.
#' p <- combinePeaks(list(p1, p2, p3), tolerance = 0.05)
#'
#' ## Plot the spectra before and after combining
#' par(mfrow = c(2, 1), mar = c(4.3, 4, 1, 1))
#' plot(p1[, 1], p1[, 2], xlim = range(mzs[5:25]), type = "h", col = "red")
#' points(p2[, 1], p2[, 2], type = "h", col = "green")
#' points(p3[, 1], p3[, 2], type = "h", col = "blue")
#'
#' plot(p[, 1], p[, 2], xlim = range(mzs[5:25]), type = "h",
#'     col = "black")
#'
#' ## With `intensityFun = max` the maximal intensity per peak is reported.
#' p <- combinePeaks(list(p1, p2, p3), tolerance = 0.05,
#'     intensityFun = max)
combinePeaks <- function(x, intensityFun = base::mean,
                         mzFun = base::mean, weighted = FALSE,
                         tolerance = 0, ppm = 0, timeDomain = FALSE,
                         main = 1L, unionPeaks = TRUE, ...) {
    if (length(x) == 1)
        return(x[[1]])
    mzs <- lapply(x, "[", , y = 1)
    mzs_lens <- lengths(mzs)
    mzs <- unlist(mzs, use.names = FALSE)
    mz_order <- order(mzs)
    mzs <- mzs[mz_order]
    if (timeDomain)
        mz_groups <- group(sqrt(mzs), tolerance = tolerance, ppm = ppm)
    else
        mz_groups <- group(mzs, tolerance = tolerance, ppm = ppm)
    ints <- unlist(lapply(x, "[", , y = 2), use.names = FALSE)[mz_order]
    if (!unionPeaks) {
        ## Want to keep only those groups with a m/z from the main spectrum.
        ## vectorized version from @sgibb
        is_in_main <- rep.int(seq_along(mzs_lens), mzs_lens)[mz_order] == main
        keep <- mz_groups %in% mz_groups[is_in_main]
        ## Keep only values for which a m/z in main is present.
        mz_groups <- mz_groups[keep]
        mzs <- mzs[keep]
        ints <- ints[keep]
    }
    if (weighted) {
        intsp <- split(ints, mz_groups)
        res <- cbind(mz = mapply(split(mzs, mz_groups), intsp,
                                 FUN = function(mz_vals, w)
                                     weighted.mean(mz_vals, w + 1,
                                                   na.rm = TRUE),
                                 USE.NAMES = FALSE, SIMPLIFY = TRUE),
                     intensity = vapply1d(intsp, FUN = intensityFun, ...))
    } else {
        res <- cbind(mz = vapply1d(split(mzs, mz_groups), FUN = mzFun, ...),
                     intensity = vapply1d(split(ints, mz_groups),
                                          FUN = intensityFun, ...))
    }
    if (is.unsorted(res[, 1]))
        stop("m/z values of combined spectrum are not ordered")
    res
}
