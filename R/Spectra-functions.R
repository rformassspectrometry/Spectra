#' @include hidden_aliases.R
NULL

.valid_processing_queue <- function(x) {
    if (length(x))
        if (!all(vapply1l(x, inherits, "ProcessingStep")))
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

#' @title Apply a function to a Spectra of length 1
#'
#' @description
#'
#' Utility function to apply any given function to a `Spectra` of length 1. If
#' `length(x) > 1` the function splits `x` and calls itself with in a `lapply`
#' call. The function `FUN` can expect to get a `Spectra` object of length 1
#' and can access all of it's variables through the accessor methods or the `$`
#' operator.
#'
#' @note
#'
#' This function supports parallel processing which is however only suggested
#' if the function applied to the `Spectra` of length 1 is computationally
#' intense.
#'
#' For faster processing, especially if only peak data is involved, consider
#' using `.peaksapply` instead, as it takes care of proper splitting the data
#' for parallel processing.
#'
#' @param x `Spectra` object.
#'
#' @param FUN `function` to be applied to each spectrum.
#'
#' @param ... Additional parameters for `FUN`.
#'
#' @param BPPARAM Optional settings for `BiocParallel`-based parallel
#'     processing. See [bpparam()] for more informations and options.
#'
#' @return A `list` with the result of `FUN`.
#'
#' @author Johannes Rainer
#'
#' @importFrom BiocParallel SerialParam
#'
#' @noRd
.lapply <- function(x, FUN = identity, ..., BPPARAM = SerialParam()) {
    res <- bplapply(split(x, seq_along(x)), FUN = FUN, ..., BPPARAM = BPPARAM)
    names(res) <- spectraNames(x)
    res
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

#' @title Spectrum comparison
#'
#' @description
#'
#' Compare each spectrum in `x` with each spectrum in `y` with function `FUN`.
#' Mapping between the peaks in both spectra can be defined with the `MAPFUN`
#' function.
#'
#' @note
#'
#' Results might be slightly different if `x` and `y` are switched, because of
#' the matching of the peaks by their m/z between spectra. Matching is performed
#' by applying a `tolerance` and `ppm` to the second spectrum (i.e. from `y`)
#' and, especially for `ppm`, the maximal accepted difference depend thus on
#' the m/z values from the second spectrum.
#'
#' @param x [Spectra()]
#'
#' @param y [Spectra()]
#'
#' @param FUN `function` to compare the (m/z matched) intensities from each
#'     spectrum in `x` with those in `y`.
#'
#' @param MAPFUN `function` to map peaks between the two compared spectra. See
#'     [joinPeaks()].
#'
#' @param tolerance `numeric(1)` allowing to define a constant maximal accepted
#'     difference between m/z values for peaks to be matched.
#'
#' @param ppm `numeric(1)` allowing to define a relative, m/z-dependent,
#'     maximal accepted difference between m/z values for peaks to be matched.
#'
#' @param ... additional parameters passed to `FUN` and `MAPFUN`.
#'
#' @return
#'
#' `matrix` with number of rows equal to length of `x` and number of columns
#' equal `length(y)`.
#'
#' @importFrom stats cor
#'
#' @examples
#'
#' library(Spectra)
#' fl <- system.file("TripleTOF-SWATH/PestMix1_SWATH.mzML", package = "msdata")
#' sps <- Spectra(fl, source = MsBackendMzR())
#'
#' sps <- sps[1:10]
#' res <- .compare_spectra(sps, sps, ppm = 20, FUN = cor,
#'     use = "pairwise.complete.obs")
#' res <- .compare_spectra(sps, sps, FUN = cor,
#'     use = "pairwise.complete.obs", method = "spearman", ppm = 40)
#'
#' res <- .compare_spectra(sps[1], sps, FUN = cor,
#'     use = "pairwise.complete.obs")
#' res <- .compare_spectra(sps, sps[1], FUN = cor,
#'     use = "pairwise.complete.obs")
#' @noRd
.compare_spectra <- function(x, y = NULL, MAPFUN = joinPeaks, tolerance = 0,
                             ppm = 20, FUN = cor, ...) {
    x_idx <- seq_along(x)
    y_idx <- seq_along(y)

    nx <- length(x_idx)
    ny <- length(y_idx)

    mat <- matrix(NA_real_, nrow = nx, ncol = ny,
                  dimnames = list(spectraNames(x), spectraNames(y)))

    ## This code duplication may be overengineering.
    if (nx >= ny) {
        for (i in x_idx) {
            px <- peaks(x[i])[[1L]]
            for (j in y_idx) {
                peak_map <- MAPFUN(px, peaks(y[j])[[1L]],
                                   tolerance = tolerance, ppm = ppm, ...)
                mat[i, j] <- FUN(peak_map[[1L]][, 2L], peak_map[[2L]][, 2L],
                                 ...)
            }
        }
    } else {
        for (j in y_idx) {
            py <- peaks(y[j])[[1L]]
            for (i in x_idx) {
                peak_map <- MAPFUN(peaks(x[i])[[1]], py,
                                   tolerance = tolerance, ppm = ppm, ...)
                mat[i, j] <- FUN(peak_map[[1L]][, 2L], peak_map[[2L]][, 2L],
                                 ...)
            }
        }
    }
    mat
}

#' @description
#'
#' Compare each spectrum in `x` with each other spectrum in `x`. Makes use of
#' `combn` to avoid calculating combinations twice.
#'
#' @inheritParams .compare_spectra
#'
#' @importFrom utils combn
#'
#' @importFrom MsCoreUtils vapply1d
#' @examples
#'
#' res <- .compare_spectra(sps, sps)
#' res_2 <- .compare_spectra_self(sps)
#' all.equal(res, res_2)    # comparison x[1], x[2] != x[2], x[1]
#' all.equal(res[!lower.tri(res)], res_2[!lower.tri(res_2)])
#'
#' @noRd
.compare_spectra_self <- function(x, MAPFUN = joinPeaks, tolerance = 0,
                                  ppm = 20, FUN = cor, ...) {
    nx <- length(x)
    m <- matrix(NA_real_, nrow = nx, ncol = nx,
                  dimnames = list(spectraNames(x), spectraNames(x)))

    cb <- which(lower.tri(m, diag = TRUE), arr.ind = TRUE)
    for (i in seq_len(nrow(cb))) {
        cur <- cb[i, 2L]
        if (i == 1L || cb[i - 1L, 2L] != cur)
            py <- px <- peaks(x[cur])[[1L]]
        else
            py <- peaks(x[cb[i, 1L]])[[1L]]
        map <- MAPFUN(px, py, tolerance = tolerance, ppm = ppm, ...)
        m[cb[i, 1L], cur] <- m[cur, cb[i, 1L]] <-
            FUN(map[[1L]][, 2L], map[[2L]][, 2L], ...)
    }
    m
}
