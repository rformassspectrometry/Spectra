#' @title Estimating precursor m/z valus for DDA data
#'
#' @description
#'
#' MS data from Waters instruments are calibrated through the *Lock Mass*, but,
#' while all m/z values of mass peaks in each spectrum will be calibrated by
#' this method, the reported precursor m/z might not. The precursor m/z in the
#' converted mzML file will have m/z values from quadrupole isolation windows
#' instead of accurate m/z values. See also the
#' [GNPS documentation](https://ccms-ucsd.github.io/GNPSDocumentation/fileconversion_waters/)
#' for more information.
#'
#' The `estimatePrecursorMz()` function estimates/adjusts the reported precursor
#' m/z of a fragment spectrum using the following approach: in data dependent
#' acquisition (DDA) mode, the MS instrument will select ions with the highest
#' intensities in one MS scan for fragmentation. Thus, for each fragment
#' spectrum, this method identifies in the previous MS1 spectrum the peak with
#' the highest intensity and an m/z value similar to the fragment spectrum's
#' reported precursor m/z (given parameters `tolerance` and `ppm`). This m/z
#' value is then reported. Since the fragment spectrum's potential MS1 mass
#' peak is selected based on its intensity, this method should **only be used
#' for DDA data**.
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
#' Users of this function should evaluate and compare the estimated precursor
#' m/z values with the originally reported one and only consider adjusted values
#' they feel comfortable with.
#'
#' @param object [Spectra()] object with **DDA** data.
#'
#' @param tolerance `numeric(1)` defining an absolute acceptable difference in
#'     m/z between the fragment spectra's reported precursor m/z and the
#'     MS1 peaks considered as the precursor peak. All MS1 peaks from the
#'     previous MS1 scan with an m/z between the fragment spectrum's
#'     precursorMz +/- (tolerance + ppm(precursorMz, ppm)) are considered.
#'
#' @param ppm `numeric(1)` defining the m/z dependent acceptable difference in
#'     m/z. See documentation of parameter `tolerance` for more information.
#'
#' @param BPPARAM parallel processing setup. Defaults to
#'     `BPPARAM = SerialParam()`. See [BiocParallel::SerialParam()] for
#'     more information.
#'
#' @return `numeric` of length equal to the number of spectra in `object` with
#'     the fragment spectra's estimated precursor m/z values. For MS1 spectra
#'     `NA_real_` values are returned. The original precursor m/z is reported
#'     for MS2 spectra for which no matching MS1 peak was found.
#'
#' @author Mar Garcia-Aloy, Johannes Rainer
#'
#' @seealso
#'
#' [addProcessing()] for other data analysis and manipulation functions.
#'
#' @export
#'
#' @examples
#'
#' ## Load a DDA test data set. For the present data set no large differences
#' ## between the reported and the *actual* precursor m/z are expected.
#' fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML", package = "msdata")
#' s <- Spectra(fl)
#'
#' pmz <- estimatePrecursorMz(s)
#'
#' ## plot the reported and estimated precursor m/z values against each other
#' plot(precursorMz(s), pmz)
#' abline(0, 1)
#'
#' ## They seem highly similar, but are they identical?
#' identical(precursorMz(s), pmz)
#' all.equal(precursorMz(s), pmz)
#'
#' ## Plot also the difference of m/z values against the m/z value
#' plot(precursorMz(s), precursorMz(s) - pmz, xlab = "precursor m/z",
#'     ylab = "difference reported - estimated precursor m/z")
#'
#' ## we could then replace the reported precursor m/z values
#' s$precursorMz <- pmz
estimatePrecursorMz <- function(object, tolerance = 0.3, ppm = 10,
                                BPPARAM = SerialParam()) {
    if (!inherits(object, "Spectra"))
        stop("'object' needs to be a 'Spectra' object.")
    f <- factor(dataOrigin(object))
    BPPARAM <- backendBpparam(object, BPPARAM = BPPARAM)
    res <- bplapply(split(object, f), .adjust_dda_precursor_mz,
                    tolerance = tolerance, ppm = ppm, BPPARAM = BPPARAM)
    unsplit(res, f = f)
}


#' @param x `Spectra` with MS1 and MS2 DDA data from a single sample/file (!)
#'
#' @param tolerance `numeric(1)` defining the half window size of m/z values
#'     to extract mass peaks from the MS1 around the MS2 spectrum's precursor
#'     m/z. The size of the final extraction m/z window is then 2 * tolerance,
#'     hence all mass peaks from the respective MS1 spectrum with an m/z value
#'     that is between the reported MS2 spectrum's precursor m/z value
#'     +/- tolerance will be considered for the precursor m/z adjustment.
#'
#' @param ppm `numeric(1)` defining an additional m/z-relative
#'     parts-per-million value of the MS2 spectrum's precursor m/z to increase
#'     the extraction m/z window defined by `tolerance` on both sides.
#'
#' @return `numeric` of length equal to `length(x)` with adjusted precursor m/z
#'     values or `NA_real_` for MS1 spectra. In case no MS1 spectrum was
#'     available for a MS2 spectrum or no mass peaks were present in that MS1
#'     spectrum within the selected m/z window an `NA_real_` value is reported.
#'
#' @author Mar Garcia-Aloy, Johannes Rainer
#'
#' @note this function is slowed down by repeatedly accessing peaksData.
#'     Changing the backend to a `MsBackendMemory` would significantly
#'     increase speed.
#'
#' @noRd
.adjust_dda_precursor_mz <- function(x, tolerance = 0.3, ppm = 10) {
    if (is.unsorted(rtime(x)))
        stop("Spectra with data origin ", dataOrigin(x[1L]),
             " are not increasingly sorted by retention time.")
    ms2_idx <- which(msLevel(x) == 2L)
    ms1_idx <- which(msLevel(x) == 1L)
    pmzs <- precursorMz(x)
    for (i in ms2_idx) {
        pmz <- NA_real_
        j <- ms1_idx[ms1_idx < i]
        if (length(j)) {
            j <- j[length(j)]
            pmz <- pmzs[i]
            s1 <- filterMzRange(
                x[j], c(pmzs[i] + c(-1, 1) * (tolerance + ppm(pmzs[i], ppm))))
            pks <- peaksData(s1, c("mz", "intensity"))[[1L]]
            if (nrow(pks))
                pmz <- pks[which.max(pks[, 2L]), 1L]
        }
        pmzs[i] <- pmz
    }
    pmzs
}
