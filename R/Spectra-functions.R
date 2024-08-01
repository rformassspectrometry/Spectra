#' @include hidden_aliases.R
NULL

.valid_processing_queue <- function(x) {
    if (length(x) && !all(vapply1l(x, inherits, "ProcessingStep")))
        return("'processingQueue' should only contain ProcessingStep objects.")
    NULL
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
#' @param spectraData `data.frame` with additional spectra variables as
#'     requested with `addProcessing`.
#'
#' @return `list` of peak `matrix` elements.
#'
#' @author Johannes Rainer
#'
#' @importFrom ProtGenerics executeProcessingStep
#'
#' @noRd
.apply_processing_queue <- function(x, spectraData = data.frame(),
                                    queue = NULL) {
    x <- unname(x)
    if (length(queue)) {
        args <- list()
        ls <- length(spectraData)
        if (ls)
            colnames(spectraData) <- sub("^msLevel$", "spectrumMsLevel",
                                         colnames(spectraData))
        for (i in seq_along(x)) {
            if (ls)
                args <- as.list(spectraData[i, ])
            for (pStep in queue) {
                x[[i]] <- do.call(pStep@FUN, args = c(x[i], args, pStep@ARGS))
            }
        }
    }
    x
}

.processingQueueVariables <- function(x) {
    if (.hasSlot(x, "processingQueueVariables"))
        x@processingQueueVariables
    else c("msLevel", "centroided")
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
#' @param spectraVariables `character` defining spectra variables that should
#'     be passed to the function.
#'
#' @param ... optional additional arguments to `FUN`.
#'
#' @param f `factor` or `vector` that can be coerced to one defining how the
#'     data should be split for parallel processing. Set to `NULL` or
#'     `factor()` to disable splitting and parallel processing.
#'
#' @param columns `character` defining the columns that should be returned.
#'     This will be passed to the backend's `peaksData` function.
#'
#' @param BPPARAM parallel processing setup.
#'
#' @return `list` of `matrix` with the peaks.
#'
#' @author Johannes Rainer
#'
#' @importMethodsFrom BiocParallel bplapply
#'
#' @note
#'
#' Who's calling `.peaksapply`:
#'
#' In *Spectra.R*:
#'
#' - `peaksData`: would pass optional parameter `f` down.
#' - `as,Spectra,list`.
#' - `intensity`: would pass optional parameter `f` down.
#' - `ionCount`.
#' - `isCentroided`.
#' - `isEmpty`.
#' - `mz`: would pass optional parameter `f` down.
#' - `lengths`.
#' - `$`.
#' - `bin`.
#'
#' In *Spectra-functions.R*:
#'
#' - `.compare_spectra_chunk`: has its own `chunkSize` parameter and ignores
#'   any other splitting. Does also some weird internal chunking by
#'   `dataStorage` that we might want to skip/disable (or use
#'   `backendParallelFactor` instead).
#' - `.compare_spectra_self`: does not do any splitting.
#' - `.combine_spectra`: uses `f = x$dataStorage` to split. Needs to be fixed.
#' - `.estimate_precursor_intensity`: disables parallel/chunk wise processing.
#'   But it's invoked splitting based on `dataOrigin`, so, should be OK.
#'
#' @noRd
.peaksapply <- function(object, FUN = NULL, ...,
                        spectraVariables = .processingQueueVariables(object),
                        f = .parallel_processing_factor(object),
                        columns = c("mz", "intensity"),
                        BPPARAM = backendBpparam(object)) {
    len <- length(object)
    lf <- length(f)
    if (!len)
        return(list())
    if (lf && lf != len)
        stop("Length of 'f' has to match 'length(object)' (", len, ")")
    if (!is.factor(f))
        f <- factor(f)
    pqueue <- object@processingQueue
    if (!is.null(FUN))
        pqueue <- c(pqueue, ProcessingStep(FUN, ARGS = list(...)))
    ll <- length(levels(f))
    if (length(pqueue) || ll > 1) {
        if (!all(c("mz", "intensity") %in% columns)) {
            scols <- columns
            columns <- union(c("mz", "intensity"), columns)
            pqueue <- c(
                pqueue, ProcessingStep(
                            FUN = function(x, scols, ...) {
                                x[, scols, drop = FALSE]
                            }, ARGS = list(scols = scols)))
        }
        .local <- function(z, queue, svars) {
            if (length(svars))
                spd <- as.data.frame(spectraData(z, columns = svars))
            else spd <- NULL
            .apply_processing_queue(peaksData(z, columns = columns), spd,
                                    queue = queue)
        }
        if (ll > 1) {
            res <- bplapply(split(object@backend, f), .local,
                            queue = pqueue, svars = spectraVariables,
                            BPPARAM = BPPARAM)
            unsplit(res, f = f, drop = TRUE)
        } else
            .local(object@backend, pqueue, spectraVariables)
    } else peaksData(object@backend, columns = columns)
}

#' @title Apply a function to subsets of Spectra
#'
#' @description
#'
#' Applying a function to a Spectra eventually splitting it into chunks for
#' parallel/serial processing. By default `x` is splitted by `dataStorage`
#' applying the function `FUN` to each chunk. This has the advantage of
#' possibly running calculations in parallel **and** requiring, depending on
#' the backend, less memory (e.g. in case of an on-disk backend only the m/z
#' and intensity values for the current chunk has to be loaded into memory).
#'
#' While being partially redundant to the `.peaksapply` function, this function
#' allows the function `FUN` to get access to ALL spectra variables.
#'
#' Note that results are returned as generated by `split`, i.e. in this order
#' and as a `list` of values. Eventually use `unlist` later.
#'
#' @param x `Spectra` object.
#'
#' @param FUN `function` to be applied to each spectrum.
#'
#' @param f factor to split `x` into chunks for parallel processing. Defaults
#'     to splitting `x` by `dataStorage`.
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
.lapply <- function(x, FUN, f = as.factor(dataStorage(x)), ...,
                    BPPARAM = SerialParam()) {
    bplapply(split(x, f), FUN = FUN, ..., BPPARAM = BPPARAM)
}

#' @export applyProcessing
#'
#' @rdname Spectra
applyProcessing <- function(object, f = processingChunkFactor(object),
                            BPPARAM = bpparam(), ...) {
    queue <- object@processingQueue
    if (!length(queue))
        return(object)
    if (isReadOnly(object@backend))
        stop(class(object@backend), " is read-only. 'applyProcessing' works ",
             "only with backends that support writing data.")
    BPPARAM <- backendBpparam(object@backend, BPPARAM)
    svars <- .processingQueueVariables(object)
    pv <- peaksVariables(object)
    if (length(f)) {
        if (!is.factor(f))
            f <- factor(f, levels = unique(f))
        if (length(f) != length(object))
            stop("length 'f' has to be equal to the length of 'object' (",
                 length(object), ")")
        bknds <- bplapply(
            split(object@backend, f = f), function(z, queue, pv, svars) {
                if (length(svars))
                    spd <- as.data.frame(spectraData(z, columns = svars))
                else spd <- NULL
                peaksData(z) <- .apply_processing_queue(
                    peaksData(z, columns = pv), spd, queue)
                z
            }, queue = queue, pv = pv, svars = svars, BPPARAM = BPPARAM)
        bknds <- backendMerge(bknds)
        if (is.unsorted(f))
            bknds <- bknds[order(unlist(split(seq_along(bknds), f),
                                        use.names = FALSE))]
        object@backend <- bknds
    } else {
        if (length(svars))
            spd <- as.data.frame(spectraData(object@backend, columns = svars))
        else spd <- NULL
        peaksData(object@backend) <- .apply_processing_queue(
            peaksData(object@backend, columns = pv), spd, queue)
    }
    object@processing <- .logging(object@processing,
                                  "Applied processing queue with ",
                                  length(object@processingQueue),
                                  " steps")
    object@processingQueue <- list()
    if (!.hasSlot(object, "processingQueueVariables"))
        object <- updateObject(object, check = FALSE)
    object@processingQueueVariables <- character()
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
    if (!any(uniqueMsLevels(object) %in% msLevel)) {
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
#' function. This function *realizes* the peak matrices in chunks to balance
#' between performance (high `chunkSize`) and low memory demand
#' (small `chunkSize`).
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
#' @param chunkSize `integer(1)` defining the number of peak matrices that
#'     should be realized (i.e. loaded into memory) in each iteration. Large
#'     `chunkSize` will result in faster processing, small `chunkSize` in
#'     lower memory demand. Defaults to `chunkSize = 10000`.
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
#' res <- .compare_spectra_chunk(sps, sps, ppm = 20)
#'
#' res <- .compare_spectra_chunk(sps[1], sps, FUN = correlateSpectra)
#' res <- .compare_spectra_chunk(sps, sps[1], FUN = correlateSpectra)
#'
#' @noRd
.compare_spectra_chunk <- function(x, y = NULL, MAPFUN = joinPeaks,
                                   tolerance = 0, ppm = 20,
                                   FUN = ndotproduct,
                                   chunkSize = 10000, ...,
                                   BPPARAM = bpparam()) {
    nx <- length(x)
    ny <- length(y)
    bppx <- backendBpparam(x@backend, BPPARAM)
    if (is(bppx, "SerialParam"))
        fx <- NULL
    else fx <- backendParallelFactor(x@backend)
    pqvx <- .processingQueueVariables(x)
    bppy <- backendBpparam(y@backend, BPPARAM)
    if (is(bppy, "SerialParam"))
        fy <- NULL
    else fy <- backendParallelFactor(y@backend)
    pqvy <- .processingQueueVariables(y)
    if (nx <= chunkSize && ny <= chunkSize) {
        mat <- .peaks_compare(
            .peaksapply(x, f = fx, spectraVariables = pqvx, BPPARAM = bppx),
            .peaksapply(y, f = fy, spectraVariables = pqvy, BPPARAM = bppy),
            MAPFUN = MAPFUN, tolerance = tolerance, ppm = ppm,
            FUN = FUN, xPrecursorMz = precursorMz(x),
            yPrecursorMz = precursorMz(y), ...)
        dimnames(mat) <- list(spectraNames(x), spectraNames(y))
        return(mat)
    }

    x_idx <- seq_len(nx)
    x_chunks <- split(x_idx, ceiling(x_idx / chunkSize))
    rm(x_idx)

    y_idx <- seq_len(ny)
    y_chunks <- split(y_idx, ceiling(y_idx / chunkSize))
    rm(y_idx)

    mat <- matrix(NA_real_, nrow = nx, ncol = ny,
                  dimnames = list(spectraNames(x), spectraNames(y)))
    for (x_chunk in x_chunks) {
        x_buff <- x[x_chunk]
        if (length(fx))
            fxc <- backendParallelFactor(x_buff@backend)
        else fxc <- NULL
        x_peaks <- .peaksapply(x_buff, f = fxc, spectraVariables = pqvx,
                               BPPARAM = bppx)
        for (y_chunk in y_chunks) {
            y_buff <- y[y_chunk]
            if (length(fy))
                fyc <- backendParallelFactor(y_buff@backend)
            else fyc <- NULL
            mat[x_chunk, y_chunk] <- .peaks_compare(
                x_peaks, .peaksapply(y_buff, f = fyc, spectraVariables = pqvy,
                                     BPPARAM = bppy),
                MAPFUN = MAPFUN, tolerance = tolerance, ppm = ppm,
                FUN = FUN, xPrecursorMz = precursorMz(x_buff),
                yPrecursorMz = precursorMz(y_buff), ...)
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
                                  ppm = 20, FUN = ndotproduct, ...) {
    nx <- length(x)
    m <- matrix(NA_real_, nrow = nx, ncol = nx,
                  dimnames = list(spectraNames(x), spectraNames(x)))

    sv <- .processingQueueVariables(x)
    cb <- which(lower.tri(m, diag = TRUE), arr.ind = TRUE)
    pmz <- precursorMz(x)
    bpp <- SerialParam()
    for (i in seq_len(nrow(cb))) {
        cur <- cb[i, 2L]
        if (i == 1L || cb[i - 1L, 2L] != cur) {
            py <- px <- .peaksapply(x[cur], spectraVariables = sv, f = NULL,
                                    BPPARAM = bpp)[[1L]]
            pmzx <- pmzy <- pmz[cur]
        } else {
            py <- .peaksapply(x[cb[i, 1L]], spectraVariables = sv, f = NULL,
                              BPPARAM = bpp)[[1L]]
            pmzy <- pmz[cb[i, 1L]]
        }
        map <- MAPFUN(px, py, tolerance = tolerance, ppm = ppm,
                      xPrecursorMz = pmzx, yPrecursorMz = pmzy,
                      .check = FALSE,...)
        m[cb[i, 1L], cur] <- m[cur, cb[i, 1L]] <-
            FUN(map[[1L]], map[[2L]], xPrecursorMz = pmzx,
                yPrecursorMz = pmzy, ppm = ppm, tolerance = tolerance, ...)
    }
    m
}

#' @description
#'
#' Combine a `Spectra` of length `n` to a `Spectra` of length `1`.
#' Takes all *spectra variables* from the first spectrum and combines the
#' peaks into a single peaks matrix.
#'
#' @note
#'
#' If the backend of the input `Spectra` is *read-only* we're first converting
#' that to a `MsBackendMemory` to support replacing peak data.
#'
#' @param x `Spectra`, ideally from a single `$dataStorage`.
#'
#' @param f `factor` to group spectra in `x`. It is supposed that input checks
#'     have been performed beforehand.
#'
#' @param FUN `function` to be applied to the list of peak matrices.
#'
#' @return `Spectra` of length 1. If the input `Spectra` uses a read-only
#'     backend that is changed to a `MsBackendMemory`.
#'
#' @author Johannes Rainer
#'
#' @noRd
#'
#' @examples
#'
#' spd <- DataFrame(msLevel = c(2L, 2L, 2L), rtime = c(1, 2, 3))
#' spd$mz <- list(c(12, 14, 45, 56), c(14.1, 34, 56.1), c(12.1, 14.15, 34.1))
#' spd$intensity <- list(c(10, 20, 30, 40), c(11, 21, 31), c(12, 22, 32))
#'
#' sps <- Spectra(spd)
#'
#' res <- .combine_spectra(sps)
#' res$mz
#' res$intensity
#'
#' res <- .combine_spectra(sps, FUN = combinePeaksData, tolerance = 0.1)
#' res$mz
#' res$intensity
.combine_spectra <- function(x, f = x$dataStorage,
                             FUN = combinePeaksData, ...) {
    if (!is.factor(f))
        f <- factor(f)
    else f <- droplevels(f)
    x_new <- x[match(levels(f), f)]
    if (isReadOnly(x_new@backend))
        x_new <- setBackend(x_new, MsBackendMemory())
    bpp <- SerialParam()
    peaksData(x_new@backend) <- lapply(
        split(.peaksapply(x, f = factor(), BPPARAM = bpp), f = f),
        FUN = FUN, ...)
    x_new@processingQueue <- list()
    validObject(x_new)
    x_new
}

#' Utility to concatenate a `list` of `Spectra` objects into a single `Spectra`.
#' This function is called in `c`.
#'
#' @param x `list` of `Spectra` objects.
#'
#' @return `Spectra`.
#'
#' @author Johannes Rainer
#'
#' @noRd
.concatenate_spectra <- function(x) {
    cls <- vapply1c(x, class)
    if (any(cls != "Spectra"))
        stop("Can only concatenate 'Spectra' objects")
    pqs <- lapply(x, function(z) z@processingQueue)
    ## For now we stop if there is any of the processingQueues not empty. Later
    ## we could even test if they are similar, and if so, merge.
    if (any(lengths(pqs)))
        stop("Can not concatenate 'Spectra' objects with non-empty ",
             "processing queue. Consider calling 'applyProcessing' before.")
    metad <- do.call(c, lapply(x, function(z) z@metadata))
    procs <- unique(unlist(lapply(x, function(z) z@processing)))
    object <- new(
        "Spectra", metadata = metad,
        backend = backendMerge(lapply(x, function(z) z@backend)),
        processing = c(procs, paste0("Merge ", length(x),
                                     " Spectra into one [", date(), "]"))
    )
    validObject(object)
    object
}

#' @export concatenateSpectra
#'
#' @rdname Spectra
concatenateSpectra <- function(x, ...) {
    .concatenate_spectra(unlist(unname(list(unname(x), ...))))
}

#' @export combineSpectra
#'
#' @rdname Spectra
combineSpectra <- function(x, f = x$dataStorage, p = x$dataStorage,
                           FUN = combinePeaksData, ..., BPPARAM = bpparam()) {
    if (!is.factor(f))
        f <- factor(f, levels = unique(f))
    if (!is.factor(p))
        p <- factor(p, levels = unique(p))
    BPPARAM <- backendBpparam(x@backend, BPPARAM)
    if (length(f) != length(x) || length(p) != length(x))
        stop("length of 'f' and 'p' have to match length of 'x'")
    if (isReadOnly(x@backend))
        message("Backend of the input object is read-only, will change that",
                " to an 'MsBackendMemory'")
    if (nlevels(p) > 1) {
        ## We split the workload by storage file. This ensures memory efficiency
        ## for file-based backends.
        res <- bpmapply(FUN = .combine_spectra, split(x, p), split(f, p),
                        MoreArgs = list(FUN = FUN, ...), BPPARAM = BPPARAM)
        .concatenate_spectra(res)
    } else
        .combine_spectra(x, f = f, FUN = FUN, ...)
}

#' @description
#'
#' Internal function to check if any (or all) of the provided `mz` values are
#' in the spectras' m/z.
#'
#' @param x `Spectra` object
#'
#' @param mz `numeric` of m/z value(s) to check in each spectrum of `x`.
#'
#' @param tolarance `numeric(1)` with the tolerance.
#'
#' @param ppm `numeric(1)` with the ppm.
#'
#' @param condFun `function` such as `any` or `all`.
#'
#' @param parallel `BiocParallel` parameter object.
#'
#' @return `logical` same length than `x`.
#'
#' @author Johannes Rainer
#'
#' @importFrom MsCoreUtils common
#'
#' @noRd
.has_mz <- function(x, mz = numeric(), tolerance = 0, ppm = 20, condFun = any,
                    parallel = SerialParam()) {
    mzs <- mz(x, BPPARAM = parallel)
    vapply(mzs, FUN = function(z)
        condFun(common(mz, z, tolerance = tolerance, ppm = ppm)), logical(1))
}

#' @description
#'
#' Same as `.has_mz` only that a different `mz` is used for each spectrum in
#' `x`. Length of `mz` is thus expected to be equal to length of `x`.
#'
#' @param mz `numeric` **same length as `x`**.
#'
#' @author Johannes Rainer
#'
#' @noRd
.has_mz_each <- function(x, mz = numeric(), tolerance = 0, ppm = 20,
                         parallel = SerialParam()) {
    mzs <- mz(x, BPPARAM = parallel)
    mapply(mzs, mz, FUN = function(z, y) {
        if (is.na(y))
            NA # if m/z is NA is better to return NA instead of FALSE
        else common(y, z, tolerance = tolerance, ppm = ppm)
    })
}


#' @export joinSpectraData
#'
#' @rdname Spectra
joinSpectraData <- function(x, y,
                            by.x = "spectrumId",
                            by.y,
                            suffix.y = ".y") {

    stopifnot(inherits(x, "Spectra"))
    stopifnot(inherits(y, "DataFrame"))
    if (missing(by.y))
        by.y <- by.x
    x_vars <- spectraVariables(x)
    y_vars <- names(y)
    if (length(by.x) != 1 | length(by.y) != 1)
        stop("'by.x' and 'by.y' must be of length 1.")
    if (!is.character(by.x) | !is.character(by.y))
        stop("'by.x' and 'by.y' must be characters.")
    if (any(duplicated(x_vars)))
        stop("Duplicated spectra variables in 'x'.")
    if (any(duplicated(y_vars)))
        stop("Duplicated names in 'y'.")
    if (!by.x %in% x_vars)
        stop("'by.x' not found in spectra variables.")
    if (!by.y %in% y_vars)
        stop("'by.y' not found.")
    ## Keep only rows that (1) have a non-NA by.y and (2) that are in
    ## x[[by.x]] (given that we do a left join here).
    y <- y[!is.na(y[[by.y]]), ]
    y <- y[y[[by.y]] %in% spectraData(x)[[by.x]], ]
    k <- match(y[[by.y]], spectraData(x)[[by.x]])
    n <- length(x)
    ## check for by.x and by.y keys duplicates
    if (anyDuplicated(spectraData(x)[[by.x]]))
        stop("Duplicates found in the 'x' key!")
    if (anyDuplicated(y[[by.y]]))
        warning("Duplicates found in the 'y' key. Only last instance will be kept!")
    ## Don't need by.y anymore
    y_vars <- y_vars[-grep(by.y, y_vars)]
    y <- y[, y_vars, drop = FALSE]
    ## Check if there are any shared x_vars and y_vars. If so, the
    ## y_vars are appended suffix.y.
    if (length(xy_vars <- intersect(x_vars, y_vars))) {
        y_vars[y_vars %in% xy_vars] <- paste0(y_vars[y_vars %in% xy_vars], suffix.y[1])
        names(y) <- y_vars
    }
    for (y_var in y_vars) {
        ## Initialise the new column as a list, vector or List (in
        ## that order!)
        if (is.list(y[[y_var]])) {
            x[[y_var]] <- vector("list", length = n)
        } else if (is.vector(y[[y_var]])) {
            x[[y_var]] <- rep(NA, n)
        } else if (inherits(y[[y_var]], "List")) {
            x[[y_var]] <- as(vector("list", length = n), class(y[[y_var]]))
        }
        ## Populate the new column
        x[[y_var]][k] <- y[[y_var]]
    }
    x
}


#' @export
#'
#' @rdname Spectra
processingLog <- function(x) {
    x@processing
}

#' estimate precursor intensities based on MS1 peak intensity. This function
#' assumes that `x` is a `Spectra` with data **from a single file/sample**.
#'
#' @author Johannes Rainer
#'
#' @importFrom stats approx
#'
#' @noRd
.estimate_precursor_intensity <- function(x, ppm = 20, tolerance = 0,
                                          method = c("previous",
                                                     "interpolation"),
                                          msLevel = 2L) {
    method <- match.arg(method)
    if (any(msLevel > 2))
        warning("Estimation of precursor intensities for MS levels > 2 is ",
                "not yet validated")
    if (!length(x))
        return(numeric())
    if (is.unsorted(rtime(x))) {
        idx <- order(rtime(x))
        x <- x[idx]
    } else idx <- integer()
    pmz <- precursorMz(x)
    pmi <- rep(NA_real_, length(pmz))
    idx_ms2 <- which(!is.na(pmz) & msLevel(x) == msLevel)
    x_ms1 <- filterMsLevel(x, msLevel = msLevel - 1L)
    pks <- .peaksapply(x_ms1, f = NULL, BPPARAM = SerialParam())
    for (i in idx_ms2) {
        rt_i <- rtime(x[i])
        before_idx <- which(rtime(x_ms1) < rt_i)
        before_int <- NA_real_
        before_rt <- NA_real_
        if (length(before_idx)) {
            before_idx <- before_idx[length(before_idx)]
            before_rt <- rtime(x_ms1)[before_idx]
            peaks <- pks[[before_idx]]
            peak_idx <- closest(pmz[i], peaks[, "mz"],
                                tolerance = tolerance,
                                ppm = ppm, duplicates = "closest",
                                .check = FALSE)
            before_int <- unname(peaks[peak_idx, "intensity"])
        }
        if (method == "interpolation") {
            after_idx <- which(rtime(x_ms1) > rt_i)
            after_int <- NA_real_
            after_rt <- NA_real_
            if (length(after_idx)) {
                after_idx <- after_idx[1L]
                after_rt <- rtime(x_ms1)[after_idx]
                peaks <- pks[[after_idx]]
                peak_idx <- closest(pmz[i], peaks[, "mz"],
                                    tolerance = tolerance,
                                    ppm = ppm, duplicates = "closest",
                                    .check = FALSE)
                after_int <- unname(peaks[peak_idx, "intensity"])
            }
            if (is.na(before_int)) {
                before_int <- after_int
            } else {
                if (!is.na(after_int))
                    before_int <- approx(c(before_rt, after_rt),
                                         c(before_int, after_int),
                                         xout = rt_i)$y
            }
        }
        pmi[i] <- before_int
    }
    if (length(idx))
        pmi <- pmi[order(idx)]
    pmi
}

#' @title Apply a function stepwise to chunks of data
#'
#' @description
#'
#' `chunkapply()` splits `x` into chunks and applies the function `FUN` stepwise
#' to each of these chunks. Depending on the object it is called, this
#' function might reduce memory demand considerably, if for example only the
#' full data for a single chunk needs to be loaded into memory at a time (e.g.,
#' for `Spectra` objects with on-disk or similar backends).
#'
#' @param x object to which `FUN` should be applied. Can be any object that
#'     supports `split`.
#'
#' @param FUN the function to apply to `x`.
#'
#' @param ... additional parameters to `FUN`.
#'
#' @param chunkSize `integer(1)` defining the size of each chunk into which `x`
#'     should be splitted.
#'
#' @param chunks optional `factor` or length equal to `length(x)` defining the
#'     chunks into which `x` should be splitted.
#'
#' @return Depending on `FUN`, but in most cases a vector/result object of
#'     length equal to `length(x)`.
#'
#' @export
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Apply a function (`sqrt`) to each element in `x`, processed in chunks of
#' ## size 200.
#' x <- rnorm(n = 1000, mean = 500)
#' res <- chunkapply(x, sqrt, chunkSize = 200)
#' length(res)
#' head(res)
#'
#' ## For such a calculation the vectorized `sqrt` would however be recommended
#' system.time(sqrt(x))
#' system.time(chunkapply(x, sqrt, chunkSize = 200))
#'
#' ## Simple example splitting a numeric vector into chunks of 200 and
#' ## aggregating the values within the chunk using the `mean`. Due to the
#' ## `unsplit` the result has the same length than the input with the mean
#' ## value repeated.
#' x <- 1:1000
#' res <- chunkapply(x, mean, chunkSize = 200)
#' length(res)
#' head(res)
chunkapply <- function(x, FUN, ..., chunkSize = 1000L, chunks = factor()) {
    lx <- length(x)
    lc <- length(chunks)
    if (lx <= length(levels(chunks)) || (!lc & lx <= chunkSize))
        return(FUN(x, ...))
    if (lc > 0) {
        if (lc != lx)
            stop("Length of 'chunks' does not match length of 'x'")
    } else chunks <- .chunk_factor(lx, chunkSize = chunkSize)
    unsplit(lapply(split(x, chunks), FUN, ...), f = chunks)
}

.chunk_factor <- function(len, chunkSize = 1000L) {
    if (len <= chunkSize)
        return(as.factor(rep(1L, len)))
    as.factor(rep(1:ceiling(len / chunkSize), each = chunkSize)[seq_len(len)])
}

#' @rdname Spectra
#'
#' @author Nir Shahaf, Johannes Rainer
#'
#' @export
deisotopeSpectra <-
    function(x, substDefinition = isotopicSubstitutionMatrix("HMDB_NEUTRAL"),
             tolerance = 0, ppm = 20, charge = 1) {
        im <- force(substDefinition)
        x@processing <- .logging(x@processing, "Deisotope spectra.")
        addProcessing(x, .peaks_deisotope, tolerance = tolerance, ppm = ppm,
                      substDefinition = im, charge = charge)
    }

#' @rdname Spectra
#'
#' @author Nir Shahaf, Johannes Rainer
#'
#' @export
reduceSpectra <- function(x, tolerance = 0, ppm = 20) {
    x@processing <- .logging(x@processing, "For groups of peaks with similar ",
                             "m/z keep the one with the highest intensity.")
    addProcessing(x, .peaks_reduce, tolerance = tolerance, ppm = ppm)
}

#' @rdname Spectra
#'
#' @author Nir Shahaf
#'
#' @export
filterPrecursorMaxIntensity <- function(x, tolerance = 0, ppm = 20) {
    pmz <- precursorMz(x)
    pmi <- precursorIntensity(x)
    if (any(!is.na(pmz)) && all(is.na(pmi)))
        message("Most/all precursor intensities are 'NA'. Please add ",
                "(estimated) precursor intensities to 'x' (see ",
                "'?estimatePrecursorIntensity' for more information and ",
                "examples.")
    idx <- order(pmz, na.last = NA)
    mz_grps <- group(pmz[idx], tolerance = tolerance, ppm = ppm)
    if (any(duplicated(mz_grps))) {
        keep <- is.na(pmz)
        keep[vapply(split(idx, as.factor(mz_grps)),
                    function(z) {
                        if (length(z) == 1L) z
                        else {
                            ## returns the first if all intensities are NA.
                            if (length(res <- z[which.max(pmi[z])])) res
                            else z[1L]
                        }
                    }, integer(1L),
                    USE.NAMES = FALSE)] <- TRUE
        x <- x[keep]
    }
    x@processing <- .logging(
        x@processing, "Filter: for groups of spectra with similar precursor ",
        "m/z, keep the one with the highest precursor intensity")
    x
}

#' @rdname Spectra
#'
#' @author Nir Shahaf
#'
#' @export
filterPrecursorIsotopes <-
    function(x, tolerance = 0, ppm = 20,
             substDefinition = isotopicSubstitutionMatrix("HMDB_NEUTRAL")) {
        pmz <- precursorMz(x)
        pmi <- precursorIntensity(x)
        if (any(!is.na(pmz)) && all(is.na(pmi)))
            message("Most/all precursor intensities are 'NA'. Please add ",
                    "(estimated) precursor intensities to 'x' (see ",
                    "'?estimatePrecursorIntensity' for more information and ",
                    "examples.")
        idx <- order(pmz, na.last = NA)
        if (length(idx)) {
            keep <- rep(TRUE, length(pmz))
            iso_grps <- isotopologues(
                cbind(pmz[idx], pmi[idx]),
                tolerance = tolerance, ppm = ppm,
                substDefinition = substDefinition)
            rem <- unlist(lapply(iso_grps, function(z) z[-1]),
                          use.names = FALSE)
            if (length(rem))
                keep[idx[rem]] <- FALSE
            x <- x[keep]
        }
        x@processing <- .logging(
            x@processing, "Filter: for groups of spectra for precursors ",
            "representing potential isotopes keep only the spectrum of the ",
            "monoisotopic precursor.")
        x
}

#' @rdname Spectra
#'
#' @author Johannes Rainer
#'
#' @export
scalePeaks <- function(x, by = sum, msLevel. = uniqueMsLevels(x)) {
    msl <- force(msLevel.)
    x <- addProcessing(x, .peaks_scale_intensities, msLevel = msl,
                       by = by, spectraVariables = "msLevel")
    x@processing <- .logging(
        x@processing, "Scale peak intensities in spectra of MS level(s) ",
        paste0(msLevel., collapse = ", "), ".")
    x
}

#' @rdname Spectra
#'
#' @export
filterPrecursorPeaks <- function(object, tolerance = 0, ppm = 20,
                                 mz = c("==", ">="),
                                 msLevel. = uniqueMsLevels(object)) {
    if (!inherits(object, "Spectra"))
        stop("'object' is expected to be an instance of class 'Spectra'.")
    if (!.check_ms_level(object, msLevel.))
        return(object)
    if (length(tolerance) != 1)
        stop("'tolerance' should be of length 1")
    if (length(ppm) != 1)
        stop("'ppm' should be of length 1")
    mz <- match.arg(mz)
    if (mz == "==") {
        FUN <- .peaks_filter_precursor_ne
        msg <- "matching precursor m/z."
    } else {
        FUN <- .peaks_filter_precursor_keep_below
        msg <- "above precursor."
    }
    object <- addProcessing(
        object, FUN, tolerance = tolerance, ppm = ppm,
        msLevel = msLevel.,
        spectraVariables = c("msLevel", "precursorMz"))
    object@processing <- .logging(
        object@processing, "Filter: remove peaks ", msg)
    object
}

#' @title Define a factor for parallel access to peaks data
#'
#' @description
#'
#' Define a factor to split a `Spectra` for parallel processing based on
#' parallel processing setup,
#'
#' If `processingChunkSize(x)` is defined and its value is smaller then
#' `length(x)`, a factor depending on this chunk size is returned.
#' Otherwise `backendParallelFactor(x)` is returned which is the optimal
#' splitting for parallel processing suggested by the backend.
#'
#' Properties/considerations:
#'
#' - in-memory backends: don't split if not requested by the user (e.g. by
#'   specifying `f` or `chunkSize`.
#' - on-disk backends: file-based backends, such as `MsBackendMzR`: perform
#'   per file parallel processing if `f` or `chunkSize` is not defined.
#'   Other on-disk backends: only if requested by the user.
#'
#' @param x `Spectra` object.
#'
#' @param chunkSize `integer` defining the size of chunks into which `x` should
#'     be split.
#'
#' @return `factor` with the desired way how to split `x` for chunk-wise
#'     processing (or `factor()` to disable).
#'
#' @author Johannes Rainer
#'
#' @noRd
.parallel_processing_factor <- function(x, chunkSize = processingChunkSize(x)) {
    lx <- length(x)
    if (is.finite(chunkSize) && chunkSize < lx)
        .chunk_factor(lx, chunkSize = chunkSize)
    else backendParallelFactor(x@backend)
}

#' @title Parallel and chunk-wise processing of `Spectra`
#'
#' @description
#'
#' Many operations on `Spectra` objects, specifically those working with
#' the actual MS data (peaks data), allow a chunk-wise processing in which
#' the `Spectra` is splitted into smaller parts (chunks) that are
#' iteratively processed. This enables parallel processing of the data (by
#' data chunk) and also reduces the memory demand  since only the MS data
#' of the currently processed subset is loaded into memory and processed.
#' This chunk-wise processing, which is by default disabled, can be enabled
#' by setting the processing chunk size of a `Spectra` with the
#' `processingChunkSize()` function to a value which is smaller than the
#' length of the `Spectra` object. Setting `processingChunkSize(sps) <- 1000`
#' will cause any data manipulation operation on the `sps`, such as
#' `filterIntensity()` or `bin()`, to be performed eventually in parallel for
#' sets of 1000 spectra in each iteration.
#'
#' Such chunk-wise processing is specifically useful for `Spectra` objects
#' using an *on-disk* backend or for very large experiments. For small data
#' sets or `Spectra` using an in-memory backend, a direct processing might
#' however be more efficient. Setting the chunk size to `Inf` will disable
#' the chunk-wise processing.
#'
#' For some backends a certain type of splitting and chunk-wise processing
#' might be preferable. The `MsBackendMzR` backend for example needs to load
#' the MS data from the original (mzML) files, hence chunk-wise processing
#' on a per-file basis would be ideal. The [backendParallelFactor()] function
#' for `MsBackend` allows backends to suggest a preferred splitting of the
#' data by returning a `factor` defining the respective data chunks. The
#' `MsBackendMzR` returns for example a `factor` based on the *dataStorage*
#' spectra variable. A `factor` of length 0 is returned if no particular
#' preferred splitting should be performed. The suggested chunk definition
#' will be used if no finite `processingChunkSize()` is defined. Setting
#' the `processingChunkSize` overrides `backendParallelFactor`.
#'
#' See the *Large-scale data handling and processing with Spectra* for more
#' information and examples.
#'
#' Functions to configure parallel or chunk-wise processing:
#'
#' - `processingChunkSize()`: allows to get or set the size of the chunks for
#'   parallel processing or chunk-wise processing of a `Spectra` in general.
#'   With a value of `Inf` (the default) no chunk-wise processing will be
#'   performed.
#'
#' - `processingChunkFactor()`: returns a `factor` defining the chunks into
#'   which a `Spectra` will be split for chunk-wise (parallel) processing.
#'   A `factor` of length 0 indicates that no chunk-wise processing will be
#'   performed.
#'
#' @note
#'
#' Some backends might not support parallel processing at all.
#' For these, the `backendBpparam()` function will always return a
#' `SerialParam()` independently on how parallel processing was defined.
#'
#' @param x `Spectra`.
#'
#' @param value `integer(1)` defining the chunk size.
#'
#' @return `processingChunkSize()` returns the currently defined processing
#'     chunk size (or `Inf` if it is not defined). `processingChunkFactor()`
#'     returns a `factor` defining the chunks into which `x` will be split
#'     for (parallel) chunk-wise processing or a `factor` of length 0 if
#'     no splitting is defined.
#'
#' @author Johannes Rainer
#'
#' @export
processingChunkSize <- function(x) {
    if (.hasSlot(x, "processingChunkSize"))
        x@processingChunkSize
    else Inf
}

#' @rdname processingChunkSize
#'
#' @export
`processingChunkSize<-` <- function(x, value) {
    if (length(value) != 1L)
        stop("'value' has to be of length 1")
    if (!.hasSlot(x, "processingChunkSize"))
        x <- updateObject(x)
    x@processingChunkSize <- value
    x
}

#' @rdname processingChunkSize
#'
#' @export
processingChunkFactor <- function(x) {
    if (!inherits(x, "Spectra"))
        stop("'x' is supposed to be a 'Spectra' object")
    .parallel_processing_factor(x)
}

#' @title Filter peaks based on spectra and peaks variable ranges
#'
#' @description
#'
#' The `filterPeaksRanges()` function allows to filter the peaks matrices of a
#' [Spectra] object using any set of range-based filters on numeric spectra
#' variables or peaks variables. These ranges can be passed to the function
#' using the `...` as `<variable name> = <range>` pairs. `<variable name>`
#' has to be an available spectra or peaks variable. `<range>` can be a
#' `numeric` of length 2 defining the lower and upper boundary, or a `numeric`
#' two-column matrix (multi-row matrices are also supported, see further
#' below). `filterPeaksRanges(s, mz = c(200, 300))` would for example reduce
#' the peaks matrices of the `Spectra` object `s` to mass peaks with an m/z
#' value between 200 and 300. `filterPeaksRanges()` returns the original
#' `Spectra` object with the filter operation added to the processing queue.
#' Thus, the filter gets **only** applied when the peaks data gets extracted
#' with `mz()`, `intensity()` or `peaksData()`. If ranges for both spectra
#' **and** peaks variables are defined, the function evaluates first whether
#' the spectra variable value for a spectrum is within the provided range and,
#' if so, applies also the peaks variable-based filter (otherwise an empty
#' peaks matrix is returned).
#'
#' If more than one spectra variable and/or peaks variable are defined, their
#' filter results are combined with a logical AND: a peak matrix is only
#' returned for a spectrum if all values of spectra variables are within the
#' provided (respective) ranges for spectra variables, and this matrix is
#' further filtered to contain only those peaks which values are within the
#' provided peaks variable ranges.
#'
#' **Filtering with multiple ranges** per spectra and peaks variables is also
#' supported: ranges can also be provided as multi-row numeric (two-column)
#' matrices. In this case, the above described procedure is applied for each
#' row separately and their results are combined with a logical OR, i.e.
#' peaks matrices are returned that match any of the conditions/filters
#' of a row. The number of rows of the provided ranges (being it for spectra
#' or peaks variables) have to match.
#'
#' **Missing value handling**: any comparison which involves a missing value
#' (being it a spectra variable value, a peaks variable value or a value
#' in one of the provided ranges) is treated as a logical `FALSE`. For
#' example, if the retention time of a spectrum is `NA` and the data is
#' filtered using a retention time range, an empty peaks matrix is returned
#' (for `keep = TRUE`, for `keep = FALSE` the full peaks matrix is returned).
#'
#' @note
#'
#' In contrast to some other *filter* functions, this function does not provide
#' a `msLevel` parameter that allows to define the MS level of spectra on which
#' the filter should be applied. The filter(s) will always be applied to
#' **all** spectra (irrespectively of their MS level). Through combination of
#' multiple filter ranges it is however possible to apply MS level-dependent
#' filters (see examples below for details).
#'
#' The filter will not be applied immediately to the data but only executed when
#' the mass peak data is accessed (through `peaksData()`, `mz()` or
#' `intensity()`) or by calling `applyProcessing()`.
#'
#' @param object A [Spectra] object.
#'
#' @param ... the ranges for the spectra and/or peaks variables. Has to be
#'     provided as `<name> = <range>` pairs with `<name>` being the name of a
#'     spectra or peaks variable (of numeric data type) and `<range>` being
#'     either a `numeric` of length 2 or a `numeric` two column matrix (see
#'     function desription above for details),
#'
#' @param keep `logical(1)` whether to keep (default) or remove peaks that
#'     match the provided range(s).
#'
#' @author Johannes Rainer
#'
#' @name filterPeaksRanges
#'
#' @export
#'
#' @examples
#'
#' ## Define a test Spectra
#' d <- data.frame(rtime = c(123.2, 134.2), msLevel = c(1L, 2L))
#' d$mz <- list(c(100.1, 100.2, 100.3, 200.1, 200.2, 300.3),
#'     c(100.3, 100.4, 200.2, 400.3, 400.4))
#' ## Use the index of the mass peak within the spectrum as index for
#' ## better illustration of filtering results
#' d$intensity <- list(c(1:6), 1:5)
#' s <- Spectra(d)
#' s
#'
#' ## Filter peaks removing all mass peaks with an m/z between 200 and 300
#' res <- filterPeaksRanges(s, mz = c(200, 300), keep = FALSE)
#' res
#'
#' ## The Spectra object has still the same length and spectra variables
#' length(res)
#' res$rtime
#'
#' ## The filter gets applied when mass peak data gets extracted, using either
#' ## `mz()`, `intensity()` or `peaksData()`. The filtered peaks data does
#' ## not contain any mass peaks with m/z values between 200 and 300:
#' peaksData(res)[[1L]]
#' peaksData(res)[[2L]]
#'
#' ## We next combine spectra and filter variables. We want to keep only mass
#' ## peaks of MS2 spectra that have an m/z between 100 and 110.
#' res <- filterPeaksRanges(s, mz = c(100, 110), msLevel = c(2, 2))
#' res
#' length(res)
#'
#' ## Only data for peaks are returned for which the spectra's MS level is
#' ## between 2 and 2 and with an m/z between 100 and 110. The peaks data for
#' ## the first spectrum, that has MS level 1, is thus empty:
#' peaksData(res)[[1L]]
#'
#' ## While the peaks matrix for the second spectrum (with MS level 2) contains
#' ## the mass peaks with m/z between 100 and 110.
#' peaksData(res)[[2L]]
#'
#' ## To keep also the peaks data for the first spectrum, we need to define
#' ## an additional set of ranges, which we define using a second row in each
#' ## ranges matrix. We use the same filter as above, i.e. keeping only mass
#' ## peaks with an m/z between 100 and 110 for spectra with MS level 2, but
#' ## add an additional row for MS level 1 spectra keeping mass peaks with an
#' ## m/z between 0 and 2000. Filter results of different rows are combined
#' ## using a logical OR, i.e. peaks matrices with mass peaks are returned
#' ## matching either the first, or the second row.
#' res <- filterPeaksRanges(s, mz = rbind(c(100, 110), c(0, 1000)),
#'     msLevel = rbind(c(2, 2), c(1, 1)))
#'
#' ## The results for the MS level 2 spectrum are the same as before, but with
#' ## the additional row we keep the full peaks matrix of the MS1 spectrum:
#' peaksData(res)[[1L]]
#' peaksData(res)[[2L]]
#'
#' ## As a last example we define a filter that keeps all mass peaks with an
#' ## m/z either between 100 and 200, or between 300 and 400.
#' res <- filterPeaksRanges(s, mz = rbind(c(100, 200), c(300, 400)))
#' peaksData(res)[[1L]]
#' peaksData(res)[[2L]]
#'
#' ## Such filters could thus be defined to restrict/filter the MS data to
#' ## specific e.g. retention time and m/z ranges.
filterPeaksRanges <- function(object, ..., keep = TRUE) {
    if (!inherits(object, "Spectra"))
        stop("'object' is expected to be a 'Spectra' object.")
    dots <- list(...)
    variables <- names(dots)
    if (!length(variables))
        return(object)
    ## check that:
    ## - variables are in spectraVariables
    pvars <- peaksVariables(object)
    svars <- spectraVariables(object)
    if (!all(variables %in% c(svars, pvars)))
        stop("Provided filter variable(s): ",
             paste0("\"", variables[!variables %in% c(svars, pvars)], "\"",
                    collapse = ", "), " are not valid spectra variables. ",
             "Use 'spectraVariables(object)' and 'peaksVariables()' to list ",
             "available variables.")
    ## - range parameters are defined correctly
    err <- paste0("Range parameters have to be either a 'numeric' of length ",
                  "2 or a 'numeric' matrix with two columns.")
    dots <- lapply(dots, function(z) {
        if (is.null(nrow(z))) {
            if (length(z) != 2)
                stop(err)
            z <- matrix(z, ncol = 2)
        }
        if (!is.matrix(z) | !is.numeric(z)) stop(err)
        z
    })
    ## - number for rows of matrices matches.
    nr <- unlist(lapply(dots, nrow), use.names = FALSE)
    if (any(nr != nr[1L]))
        stop("Number of rows of the range matrices have to match.")
    ## OK, now proceed to split by svar and pvar and pass to the peaks function.
    pvars <- intersect(variables, pvars)
    svars <- intersect(variables, svars)
    object <- addProcessing(object, .peaks_filter_ranges, ranges = dots,
                            svars = svars, pvars = pvars,
                            spectraVariables = c(svars, "msLevel"), keep = keep)
    if (keep) keep_or_remove <- "select"
    else keep_or_remove <- "remove"
    object@processing <- .logging(
        object@processing, "Filter: ", keep_or_remove, " peaks based on ",
        "user-provided ranges for ", length(variables), " variables")
    object
}
