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

#' @export applyProcessing
#'
#' @rdname Spectra
applyProcessing <- function(object, f = dataStorage(object),
                            BPPARAM = bpparam(), ...) {
    if (!length(object@processingQueue))
        return(object)
    if (isReadOnly(object@backend))
        stop(class(object@backend), " is read-only. 'applyProcessing' works ",
             "only with backends that support writing data.")
    if (!is.factor(f))
        f <- factor(f, levels = unique(f))
    if (length(f) != length(object))
        stop("length 'f' has to be equal to the length of 'object' (",
             length(object), ")")
    bknds <- bplapply(split(object@backend, f = f), function(z, queue) {
        peaks(z) <- .apply_processing_queue(peaks(z), msLevel(z),
                                            centroided(z), queue)
        z
    }, queue = object@processingQueue, BPPARAM = BPPARAM)
    bknds <- backendMerge(bknds)
    if (is.unsorted(f))
        bknds <- bknds[order(unlist(split(seq_along(bknds), f),
                                    use.names = FALSE))]
    object@backend <- bknds
    object@processing <- .logging(object@processing,
                                  "Applied processing queue with ",
                                  length(object@processingQueue),
                                  " steps")
    object@processingQueue <- list()
    object
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
