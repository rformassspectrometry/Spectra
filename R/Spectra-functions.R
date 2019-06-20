#' @include hidden_aliases.R
NULL

.valid_processing_queue <- function(x) {
    if (length(x))
        if (!all(vapply(x, inherits, logical(1), "ProcessingStep")))
            return("'processingQueue' should only contain ProcessingStep objects.")
    NULL
}

#' @export addProcessing
#'
#' @rdname Spectra
addProcessing <- function(object, FUN, ...) {
    if (missing(FUN))
        return(object)
    object@processingQueue <- c(object@processingQueue,
                                list(ProcessingStep(FUN, ARGS = list(...))))
    validObject(object)
    object
}

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
.remove_peaks <- function(x, spectrumMsLevel, centroided = NA, t = "min",
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
#' @inheritParams .remove_peaks
#'
#' @return `matrix` with columns `"mz"` and `"intensity"`.
#'
#' @noRd
.clean_peaks <- function(x, spectrumMsLevel, all = FALSE,
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

#' @description
#'
#' Wrapper to allow use of the `isCentroided` function with the `.peaksapply`
#' function.
#'
#' @noRd
.is_centroided_peaks <- function(x, spectrumMsLevel, centroided = NA, ...) {
    .isCentroided(x, ...)
}

#' @description
#'
#' Apply the lazy evaluation processing queue to the peaks (i.e. m/z, intensity
#' matrix) and return the result.
#'
#' @param x `list` of peak `matrix`.
#'
#' @param msLevel `integer` with the MS level of each spectrum. Length has to
#'     match `length(x)`.
#'
#' @param centroided `logical` whether spectra are centroided. Length has to
#'     match `length(x)`.
#'
#' @param queue `list` of `ProcessStep` elements.
#'
#' @return `list` of peak `matrix` elements.
#'
#' @author Johannes Rainer
#'
#' @noRd
.apply_processing_queue <- function(x, msLevel, centroided, queue = NULL) {
    if (length(queue)) {
        for (i in seq_along(x)) {
            for (pStep in queue) {
                x[[i]] <- executeProcessingStep(pStep, x[[i]],
                                                spectrumMsLevel = msLevel[i],
                                                centroided = centroided[i])
            }
        }
    }
    x
}

#' @title Apply arbitrary functions and processing queue to peaks matrices
#'
#' @description
#'
#' This function applies the processing queue and an arbitrary function to
#' the peaks matrix of each spectrum of the `Spectra` object `object`.
#'
#' @param object `Spectra` object.
#'
#' @param FUN optional function to be applied to the peaks matrix. The peaks
#'     matrix will be passed to the first argument of the function. The function
#'     should also have arguments `...`.
#'
#' @param ... optional additional arguments to `FUN`.
#'
#' @param f `factor` or `vector` that can be coerced to one defining how the
#'     data should be split for parallel processing.
#'
#' @param BPPARAM parallel processing setup.
#'
#' @return `list` of `matrix` with the peaks.
#'
#' @author Johannes Rainer
#'
#' @importMethodsFrom BiocParallel bplapply
#'
#' @noRd
.peaksapply <- function(object, FUN = NULL, ..., f = dataStorage(object),
                        BPPARAM = bpparam()) {
    len <- length(object)
    if (!len)
        return(list())
    if (length(f) != len)
        stop("Length of 'f' has to match 'length(object)' (", len, ")")
    if (!is.factor(f))
        f <- factor(f)
    pqueue <- object@processingQueue
    if (!is.null(FUN))
        pqueue <- c(pqueue, ProcessingStep(FUN, ARGS = list(...)))
    ## Question whether we would use a slim version of the backend, i.e.
    ## reduce it to certain columns/spectra variables.
    res <- bplapply(split(object@backend, f), function(z, queue) {
        .apply_processing_queue(peaks(z), msLevel(z), centroided(z),
                                queue = queue)
    }, queue = pqueue, BPPARAM = BPPARAM)
    unsplit(res, f = f, drop = TRUE)
    ## unlist(res, recursive = FALSE, use.names = FALSE)
}

#' @description
#'
#' Simple helper function to test parameter msLevel. Returns `TRUE` if parameter
#' is OK, `FALSE` if a warning is thrown and throws an error if it is not
#' a numeric.
#'
#' @noRd
.check_ms_level <- function(object, msLevel) {
    if (!length(object))
        return(TRUE)
    if (!is.numeric(msLevel))
        stop("'msLevel' must be numeric")
    if (!any(msLevel(object) %in% msLevel)) {
        warning("Specified MS levels ", paste0(msLevel, collapse = ","),
                " not available in 'object'")
        FALSE
    } else TRUE
}
