#' @title Calculating Precursor Purity for MS2 spectra
#'
#' @description
#'
#' MS instruments generally collect precursor ions in a discrete *m/z*
#' *isolation window* before fragmenting them and recording the respective
#' fragment (MS2) spectrum. Ideally, only a single ion species is fragmented,
#' depending also on the size of the isolation window, different ions (with
#' slightly different *m/z*) might be fragmented. The resulting MS2 spectrum
#' might thus contain fragments from different ions and hence be less *pure*.
#'
#' The `precursorPurity()` function calculates the **precursor purity** of MS2
#' (fragment) spectra expressed as the ratio between the itensity of the highest
#' signal in the isolation window to the sum of intensities of all MS1 peaks in
#' the isolation window. This is similar to the calculation performed in the
#' [*msPurity*](https://www.bioconductor.org/packages/release/bioc/html/msPurity.html)
#' Bioconductor package.
#'
#' The peak intensities within the isolation window is extracted from the last
#' MS1 spectrum before the respective MS2 spectrum. The spectra are thus
#' expected to be ordered by retention time. For the isolation window either the
#' isolation window reported in the `Spectra` object is used, or it is
#' calculated based on the MS2 spectra's precursor m/z. By default, the
#' isolation window is calculated based on the precursor m/z and parameters
#' `tolerance` and `ppm`: precursorMz +/- (`tolerance` +
#' `ppm(precursorMz, ppm)`). If the actually used precursor isolation window is
#' defined and available in the `Spectra` object, it can be used instead by
#' setting `useReportedIsolationWindow = TRUE` (default is
#' `useReportedIsolationWindow = FALSE`). Note that parameters `tolerance`
#' and `ppm` are ignored for `useReportedIsolationWindow = TRUE`.
#'
#' @note
#'
#' This approach is applicable only when fragment spectra are obtained through
#' data-dependent acquisition (DDA), as it assumes that the peak with the
#' highest intensity within the given isolation m/z window (from the previous
#' MS1 spectrum) corresponds to the precursor ion.
#'
#' The spectra in `object` have to be ordered by their retention time.
#'
#' @param object [Spectra()] object with LC-MS/MS data.
#'
#' @param tolerance `numeric(1)` defining an absolute value (in Da) to be used
#'     to define the isolation window. For the precursor purity calculation of
#'     an MS2 spectrum, all MS1 peaks from the previous MS1 scan with an m/z
#'     between the fragment spectrum's precursorMz +/- (tolerance +
#'     ppm(precursorMz, ppm)) are considered.
#'
#' @param ppm `numeric(1)` defining the m/z dependent acceptable difference in
#'     m/z. See documentation of parameter `tolerance` for more information.
#'
#' @param useReportedIsolationWindow `logical(1)` whether the reported
#'     isolation window, defined by spectra variables `isolationWindowLowerMz`
#'     and `isolationWindowUpperMz` in the input [Spectra] object, should be
#'     used instead of calculating the isolation window from the reported
#'     precursor m/z and parameters `tolerance` and `ppm`. Only few
#'     manufacturers report the isolation window with the spectra variables
#'     `isolationWindowLowerMz` and `isolationWindowTargetMz`, thus the default
#'     for this parameter is `FALSE`.
#'
#' @param BPPARAM parallel processing setup. Defaults to
#'     `BPPARAM = SerialParam()`. See [BiocParallel::SerialParam()] for
#'     more information.
#'
#' @return `numeric` vector of length equal to the number of spectra in
#'     `object`,  with values representing the calculated precursor purity
#'     for each spectrum. For MS1 spectra, `NA_real_` is returned. For MS2
#'     spectra, the purity is defined as the proportion of maximum signal to
#'     the total ion current within the isolation  window that is attributable
#'     to the selected precursor ion. If no matching  MS1 scan is found or
#'     the precursor peak is missing, `NA_real_` is returned.
#'
#' @author Ahlam Mentag, Johannes Rainer
#'
#' @seealso
#'
#' [addProcessing()] for other data analysis and manipulation functions.
#'
#' @export
#'
#' @examples
#'
#' ## Load a test DDA file
#' library(msdata)
#' fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML",
#'                  package = "msdata")
#' sps_dda <- Spectra(fl)
#'
#' ## Define the isolation window based on the MS2 spectra's precursor m/z
#' ## and parameter `tolerance`: isolation window with size 1Da:
#' pp <- precursorPurity(sps_dda, tolerance = 0.5)
#'
#' ## values for MS1 spectra are NA
#' head(pp[msLevel(sps_dda) == 1])
#'
#' head(pp[msLevel(sps_dda) == 2])
#'
#' ## Use the reported isolation window (if defined in the `Spectra`):
#' filterMsLevel(sps_dda, 2L) |>
#'     isolationWindowLowerMz() |>
#'     head()
#' filterMsLevel(sps_dda, 2L) |>
#'     isolationWindowUpperMz() |>
#'     head()
#'
#' pp_2 <- precursorPurity(sps_dda, useReportedIsolationWindow = TRUE)
#'
#' head(pp_2[msLevel(sps_dda) == 2])
precursorPurity <- function(object, tolerance = 0.05, ppm = 0,
                            useReportedIsolationWindow = FALSE,
                            BPPARAM = SerialParam()) {
    if (!inherits(object, "Spectra"))
        stop("'object' needs to be a 'Spectra' object.")
    f <- factor(dataOrigin(object))
    BPPARAM <- backendBpparam(object)
    res <- bplapply(split(object, f), .precursorPurity,
                    tolerance = tolerance, ppm = ppm,
                    useIsolationWindow = useReportedIsolationWindow,
                    BPPARAM = BPPARAM)
    unsplit(res, f = f)
}

#' @param x `Spectra` with MS1 and MS2 DDA data from a single sample/file (!)
#'
#' @param tolerance `numeric(1)` defining the half window size of m/z values
#'     to extract mass peaks from the MS1 spectrum around the MS2 spectrum's
#'     precursor m/z. The full m/z window is defined as 2 * tolerance, and
#'     peaks within this window (adjusted by `ppm`) are used to estimate purity.
#'
#' @param ppm `numeric(1)` defining an additional m/z-relative
#'     parts-per-million value to extend the extraction m/z window on both sides
#'     of the precursor m/z.
#'
#' @param useIsolationWindow `logical(1)` whether the actually reported
#'     isolation window lower and upper m/z from the input spectra object
#'     should be used **instead** of `tolerance` and `ppm`.
#'
#' @return `numeric` vector of length equal to `length(x)`, with estimated
#'     precursor purity values for MS2 spectra and `NA_real_` for MS1 spectra.
#'     The purity is calculated as the ratio of the intensity of the most
#'     intense peak to the total intensity within the isolation window in
#'     the preceding MS1 spectrum. If no MS1 spectrum or peaks are found
#'     within the window, `NA_real_` is returned.
#'
#' @author Ahlam Mentag, Johannes Rainer
#'
#' @importFrom MsCoreUtils ppm
#'
#' @importFrom MsCoreUtils between
#'
#' @noRd
.precursorPurity <- function(x, tolerance = 0.3, ppm = 0,
                             useIsolationWindow = FALSE) {
    if (is.unsorted(rtime(x)))
        stop("Spectra with data origin ", dataOrigin(x[1L]),
             " are not increasingly sorted by retention time.")
    ## Get pairs of MS1 and MS2 spectra
    ms2_idx <- which(msLevel(x) == 2L)
    ms1_all <- which(msLevel(x) == 1L)
    ms1_idx <- vapply(ms2_idx, function(i) max(ms1_all[ms1_all < i]), NA_real_)
    ms2_idx <- ms2_idx[is.finite(ms1_idx)]
    ms1_idx <- ms1_idx[is.finite(ms1_idx)]

    ratios <- rep(NA_real_, length(x))
    if (length(ms1_idx)) {
        if (useIsolationWindow) {
            l <- isolationWindowLowerMz(x)[ms2_idx]
            u <- isolationWindowUpperMz(x)[ms2_idx]
        } else {
            pmzs <- precursorMz(x[ms2_idx])
            if (ppm != 0)
                ppms <- ppm(pmzs, ppm)
            else ppms <- 0
            l <- pmzs - (tolerance + ppms)
            u <- pmzs + (tolerance + ppms)
        }
        pks <- peaksData(x[ms1_idx], c("mz", "intensity"), return.type = "list")
        for (i in seq_along(ms1_idx)) {
            p <- pks[[i]]
            p <- p[between(p[, 1L], c(l[i], u[i])), , drop = FALSE]
            if (nrow(p)) {
                intensities <- p[, 2L]
                ratio <- max(intensities) / sum(intensities)
            } else ratio <- NA_real_
            ratios[ms2_idx[i]] <- ratio
        }
    }
    ratios
}
