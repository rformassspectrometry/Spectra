#' @include hidden_aliases.R
NULL

#' @title The Spectra class to manage and access MS data
#'
#' @aliases Spectra-class
#'
#' @name Spectra
#'
#' @description
#'
#' The `Spectra` class encapsules spectral mass spectrometry data and
#' related metadata.
#'
#' It supports multiple data backends, e.g. in-memory ([MsBackendDataFrame()]),
#' on-disk as mzML ([MsBackendMzR()]).
#'
#' @details
#'
#' The `Spectra` class uses by default a lazy data manipulation strategy,
#' i.e. data manipulations such as performed with `removePeaks` are not applied
#' immediately to the data, but applied on-the-fly to the spectrum data once it
#' is retrieved.
#'
#' @section Creation of objects, conversion and changing the backend:
#'
#' `Spectra` classes can be created with the `Spectra` constructor function
#' which supports the following formats:
#'
#' - parameter `object` is a `DataFrame` containing the spectrum data. The
#'   provided `backend` (by default a [MsBackendDataFrame-class]) will be
#'   initialized with that data.
#' - parameter `object` is missing, in which case it is supposed that the data
#'   is provided by the [MsBackend-class] class passed along with the `backend`
#'   argument.
#'
#' `Spectra` classes are usually created with the `readSpectra`
#' function that reads general spectrum metadata information from the  mass
#' spectrometry data files.
#'
#' @section Accessing data:
#'
#' - `fromFile`: get the file/sample assignment of each spectrum. Returns an
#'   integer vector of length equal to the number of spectra.
#'
#' - `length`: get the number of spectra in the object.
#'
#' - `msLevel`: get the spectra's MS level. Returns an integer vector (names
#'   being spectrum names, length equal to the number of spectra) with the MS
#'   level for each spectrum.
#'
#' @section Data manipulation and analysis methods:
#'
#' Many data manipulation operations, such as those listed in this section, are
#' not applied immediately to the spectra, but added to a
#' *lazy processinq queue*. Operations stored in this queue are applied
#' on-the-fly to spectra data each time it is accessed. This lazy
#' execution guarantees the same functionality for `Spectra` objects with
#' any backend, i.e. backends supporting to save changes to spectrum data
#' ([MsBackendDataFrame()] as well as read-only backends (such
#' as the [MsBackendMzR()]).
#'
#' - `clean`: remove 0-intensity data points. For `all = FALSE` (the default)
#'   0-intensity peaks next to non-zero intensity peaks are retained while with
#'   `all = TRUE` all 0-intensity peaks are removed.
#'
#' - `removePeaks`: *remove* peaks lower or equal to a threshold intensity
#'   value `t` by setting their intensity to `0`. With the default `t = "min"`
#'   all peaks with an intensity smaller or equal to the minimal non-zero
#'   intensity is set to `0`. If the spectrum is in profile mode, ranges of
#'   successive non-0 peaks <= `t` are set to 0. If the spectrum is centroided,
#'   then individual peaks <= `t` are set to 0. Note that the number of peaks
#'   is not changed unless `clean` is called after `removePeaks`.
#'
#' @return See individual method description for the return value.
#'
#' @param all for `clean`: `logical(1)` whether all 0 intensity peaks should be
#'     removed (`TRUE`) or whether 0-intensity peaks directly adjacent to a
#'     non-zero intensity peak should be kept (`FALSE`).
#'
#' @param backend For `Spectra`: [MsBackend-class] to be used as backend. See
#'     section on creation of `Spectra` objects for details.
#'
#' @param BPPARAM Parallel setup configuration. See [bpparam()] for more
#'     information. This is passed directly to the [backendInitialize()] method
#'     of the [MsBackend-class].
#'
#' @param metadata For `Spectra`: optional `list` with metadata information.
#'
#' @param msLevel. `integer` defining the MS level(s) of the spectra to which
#'     the function should be applied. For `filterMsLevel`: the MS level to
#'     which `object` should be subsetted.
#'
#' @param object For `Spectra`: either a `DataFrame` or `missing`. See section
#'     on creation of `Spectra` objects for details. For all other methods a
#'     `Spectra` object.
#'
#' @param processingQueue For `Spectra`: optional `list` of
#'     [ProcessingStep-class] objects.
#'
#' @param t for `removePeaks`: a `numeric(1)` defining the threshold or `"min"`.
#'
#' @param x A `Spectra` object.
#'
#' @param ... Additional arguments.
#'
#' @author Sebastian Gibb, Johannes Rainer
#'
#' @md
#'
#' @exportClass Spectra
#'
#' @exportMethod Spectra
#'
#' @examples
#'
#' ## Create a Spectra providing a `DataFrame` containing the spectrum data.
#'
#' spd <- DataFrame(msLevel = c(1L, 2L), rtime = c(1.1, 1.2))
#' spd$mz <- list(c(100, 103.2, 104.3, 106.5), c(45.6, 120.4, 190.2))
#' spd$intensity <- list(c(200, 400, 34.2, 17), c(12.3, 15.2, 6.8))
#'
#' data <- Spectra(spd)
#' data
#'
#' msLevel(data)
NULL

#' The Spectra class
#'
#' The [Spectra-class] encapsulates data and meta-data for mass
#' spectrometry experiments.
#'
#'
#' @slot backend A derivate of [MsBackend-class] holding/controlling the spectra
#' data.
#' @slot processingQueue `list` of `ProcessingStep` objects.
#' @slot processing A `character` storing logging information.
#' @slot metadata A `list` storing experiment metadata.
#' @slot version A `characher(1)` containing the class version.
#'
#' @name Spectra-class
#' @docType class
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#'
#' @importClassesFrom S4Vectors DataFrame
#'
#' @importMethodsFrom S4Vectors lapply
#'
#' @importFrom S4Vectors DataFrame
#'
#' @noRd
setClass(
    "Spectra",
    slots = c(
        backend = "MsBackend",
        processingQueue = "list",
        ## logging
        processing = "character",
        ## metadata
        metadata = "list",
        version = "character"
    ),
    prototype = prototype(version = "0.1")
)

setValidity("Spectra", function(object) {
    msg <- .valid_processing_queue(object@processingQueue)
    if (length(msg)) msg
    else TRUE
})

#' @rdname hidden_aliases
#'
#' @importMethodsFrom methods show
#'
#' @importFrom utils capture.output
#'
#' @exportMethod show
setMethod("show", "Spectra",
    function(object) {
        cat("MSn data (", class(object)[1L], ") with ",
            length(object@backend), " spectra in a ", class(object@backend),
            " backend:\n", sep = "")
        if (length(object@backend)) {
            txt <- capture.output(show(object@backend))
            cat(txt[-1], sep = "\n")
        }
        if (length(object@processingQueue))
            cat("Lazy evaluation queue:", length(object@processingQueue),
                "processing step(s)\n")
        cat("Processing:\n", paste(object@processing, collapse="\n"), "\n")
    })

#' @rdname Spectra
setMethod("Spectra", "DataFrame", function(object, processingQueue = list(),
                                           metadata = list(), ...,
                                           backend = MsBackendDataFrame(),
                                           BPPARAM = bpparam()) {
    object$fromFile <- rep(1L, nrow(object))
    new("Spectra", metadata = metadata, processingQueue = processingQueue,
        backend = backendInitialize(
            backend, files = if (nrow(object)) NA_character_ else character(),
            object, BPPARAM = BPPARAM))
})
#' @rdname Spectra
setMethod("Spectra", "missing", function(object, processingQueue = list(),
                                         metadata = list(), ...,
                                         backend = MsBackendDataFrame(),
                                         BPPARAM = bpparam()) {
    new("Spectra", metadata = metadata, processingQueue = processingQueue,
        backend = backend)
})

## ACCESSOR METHODS

#' @rdname Spectra
setMethod("fromFile", "Spectra", function(object) fromFile(object@backend))

#' @rdname Spectra
#'
#' @exportMethod length
setMethod("length", "Spectra", function(x) length(x@backend))

#' @rdname Spectra
setMethod("msLevel", "Spectra", function(object) msLevel(object@backend))


## DATA MANIPULATION METHODS

#' @rdname Spectra
#'
#' @exportMethod removePeaks
setMethod("removePeaks", "Spectra", function(object, t = "min", msLevel.) {
    if (!is.numeric(t) & t != "min")
        stop("Argument 't' has to be either numeric of 'min'.")
    if (missing(msLevel.))
        msLevel. <- base::sort(unique(msLevel(object)))
    if (!is.numeric(msLevel.))
        stop("'msLevel' must be numeric.")
    object <- addProcessingStep(object, .remove_peaks, t = t,
                                msLevel. = msLevel.)
    object@processing <- c(object@processing,
                           paste0("Signal <= ", t, " in MS level(s) ",
                                  paste0(msLevel., collapse = ", "),
                                  " set to 0 [", date(), "]"))
    object
})

#' @rdname Spectra
#'
#' @exportMethod clean
setMethod("clean", "Spectra", function(object, all = FALSE, msLevel.) {
    if (!is.logical(all) || length(all) != 1)
        stop("Argument 'all' must be a logical of length 1")
    if (missing(msLevel.))
        msLevel. <- base::sort(unique(msLevel(object)))
    if (!is.numeric(msLevel.))
        stop("'msLevel' must be numeric.")
    object <- addProcessingStep(object, .clean_peaks, all = all,
                                msLevel. = msLevel.)
    object@processing <- c(object@processing,
                           paste0("Spectra of MS level(s) ",
                                  paste0(msLevel., collapse = ", "),
                                  " cleaned [", date(), "]"))
    object
})

## #' @rdname Spectra
## #'
## #' @importMethodsFrom BiocParallel bplapply
## #'
## #' @importFrom BiocParallel bpparam
## #'
## #' @export
## readSpectra <- function(file, sampleData, backend = BackendMzR(),
##                               smoothed = NA, metadata = list(), ...,
##                               BPPARAM = bpparam()) {
##     ## if (missing(backend) || !inherits(backend))
##     if (missing(file) || length(file) == 0)
##         stop("Parameter 'file' is required")
##     if (!all(file.exists(file)))
##         stop("Input file(s) can not be found")
##     file <- normalizePath(file)
##     if (!missing(sampleData)) {
##         if (is.data.frame(sampleData))
##             sampleData <- DataFrame(sampleData)
##     } else {
##         sampleData <- DataFrame(sampleIdx = seq_along(file))
##     }
##     if (!is.logical(smoothed))
##         stop("smoothed should be a logical")
##     .read_file <- function(z, files, smoothed) {
##         file_number <- match(z, files)
##         suppressPackageStartupMessages(
##             require("MSnbase", quietly = TRUE, character.only = TRUE))
##         msd <- mzR::openMSfile(z)
##         on.exit(mzR::close(msd))
##         hdr <- mzR::header(msd)
##         sp_idx <- seq_len(nrow(hdr))
##         rownames(hdr) <- formatFileSpectrumNames(fileIds = file_number,
##                                                  spectrumIds = seq_along(sp_idx),
##                                                  nSpectra = length(sp_idx),
##                                                  nFiles = length(files))
##         ## rename totIonCurrent and peaksCount, as detailed in
##         ## https://github.com/lgatto/MSnbase/issues/105#issuecomment-229503816
##         names(hdr) <- sub("peaksCount", "originalPeaksCount", names(hdr))
##         ## Add also:
##         ## o fileIdx -> links to fileNames property
##         ## o spIdx -> the index of the spectrum in the file.
##         hdr$fileIdx <- file_number
##         hdr$spIdx <- sp_idx
##         hdr$smoothed <- smoothed
##         if (isCdfFile(z)) {
##             if (!any(colnames(hdr) == "polarity"))
##                 hdr$polarity <- NA
##         }
##         ## Order the fdData by acquisitionNum to force use of acquisitionNum
##         ## as unique ID for the spectrum (issue #103). That way we can use
##         ## the spIdx (is the index of the spectrum within the file) for
##         ## subsetting and extracting.
##         if (!all(sort(hdr$acquisitionNum) == hdr$acquisitionNum))
##             warning(paste("Unexpected acquisition number order detected.",
##                           "Please contact the maintainers or open an issue",
##                           "on https://github.com/lgatto/MSnbase.",
##                           sep = "\n")) ## see issue #160
##         hdr[order(hdr$acquisitionNum), ]
##     }
##     spectraData <- DataFrame(
##         do.call(rbind, bplapply(file, .read_file, files = file,
##                                 smoothed = smoothed, BPPARAM = BPPARAM)))
##     msnexp <- new("Spectra",
##         backend = backendInitialize(BackendMzR(), file),
##         sampleData = sampleData,
##         spectraData = spectraData,
##         processingQueue = list(),
##         metadata = metadata,
##         processing = paste0("Data loaded [", date(), "]")
##     )

##     if (!inherits(backend, "BackendMzR"))
##         msnexp <- setBackend(msnexp, backend, ..., BPPARAM = BPPARAM)

##     msnexp
## }

## #' @rdname Spectra
## setGeneric("setBackend", function(object, backend, ..., BPPARAM = bpparam())
##     standardGeneric("setBackend"))
## #' @rdname hidden_aliases
## #'
## #' @importMethodsFrom BiocParallel bpmapply
## setMethod(
##     "setBackend",
##     c("Spectra", "MsBackend"),
##     function(object, backend, ..., BPPARAM = bpparam()) {
##     backend <- backendInitialize(backend, fileNames(object), object@spectraData,
##                                  ...)
##     ## update fileIdx, useful to split src backends across cores
##     spd <- object@spectraData
##     spd$fileIdx <- 1L
##     spd <- split(spd, object@spectraData$fileIdx)

##     ## keep current modCount
##     backend@modCount <- object@backend@modCount

##     backendSplitByFile(backend, object@spectraData) <-
##         bpmapply(function(dst, src, spd, queue) {
##             backendWriteSpectra(
##                 dst, backendReadSpectra(src, spd), spd, updateModCount=FALSE
##             )
##         },
##         dst = backendSplitByFile(backend, object@spectraData),
##         src = backendSplitByFile(object@backend, object@spectraData),
##         spd = spd,
##         SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM)

##     object@backend <- backend
##     object@processing <- c(object@processing,
##                            paste0("Backend set to '", class(backend),
##                                   "' [", date(), "]"))
##     validObject(object)
##     object
## })

## #' @rdname hidden_aliases
## setMethod(
##     "setBackend",
##     c("Spectra", "BackendMzR"),
##     function(object, backend, ..., BPPARAM = bpparam()) {
##     if (any(object@backend@modCount))
##         stop("Can not change backend to 'BackendMzR' because the ",
##              "data was changed.")
##     object@backend <- backendInitialize(backend, fileNames(object),
##                                         object@spectraData, ...)
##     object@processing <- c(object@processing,
##                            paste0("Backend set to 'BackendMzR' [", date(), "]"))
##     validObject(object)
##     object
## })
