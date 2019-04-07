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
#' It supports multiple data backends, e.g. in-memory ([MsBackend-class]),
#' on-disk as mzML ([MsBackendMzR-class]) or HDF5 ([MsBackendHdf5-class]).
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
#' `Spectra` classes are usually created with the `readSpectra`
#' function that reads general spectrum metadata information from the  mass
#' spectrometry data files.
#'
#' Alternatively it is possible to create a new object from a list of `Spectrum`
#' objects using the `Spectra` function. Additional spectrum metadata
#' columns can be provided with the `spectraData` argument, sample annotations
#' with the `sampleData` argument and arbitrary metadata with the `metadata`
#' argument. Note that objects created with the `Spectra` constructor
#' function can not use the `BackendMzR` as backend.
#'
#' `Spectra` objects can be converted to a `list` or
#' [S4Vectors::List-class] of `Spectrum` objects with the `as(object, "list")`
#' and `as(object, "List")` function, respectively.
#'
#' The [MsBackend-class] can be changed with the `setBackend` function by
#' specifying the new [MsBackend-class] with the `backend` parameter. See
#' examples for more details.
#'
#' @return See individual method description for the return value.
#'
#' @author Sebastian Gibb, Johannes Rainer
#'
#' @md
#'
#' @exportClass Spectra
#'
#' @examples
#'
#' ## Create an Spectra from a list of Spectrum objects
NULL

#' The Spectra class
#'
#' The [Spectra-class] encapsulates data and meta-data for mass
#' spectrometry experiments.
#'
#'
#' @slot backend A derivate of [Backend-class] holding/controlling the spectra
#' data.
#' @slot spectraData A [S4Vectors::DataFrame-class] storing spectra metadata.
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
    prototype = prototype(version = "0.1"),
    validity = .valid_spectra
)

## #' @rdname Spectra

## Constructor:
## provide spectraData, backend, metadata, spectra, ...
#' @param spectraData `DataFrame` with spectrum properties. Check content for
#'     specific columns.
#'
#' @param backend `Backend` class. Either `spectra` or `Backend` have to be
#'     provided.
#'
#' @param spectra `list` of `data.frame`s with columns named `"mz"` and
#'     `"intensity"`.
#'
#' @param metadata `list` of optional additional metadata.
Spectra <- function(spectraData = DataFrame(), backend = MsBackend()) {
}

## Spectra <- function(x, spectraData, metadata, ...) {
##     if (missing(x) || any(vapply(x, function(z) !inherits(z, "Spectrum"),
##                                  logical(1))))
##         stop("'x' has to be a list of 'Spectrum' objects")
##     fdata <- DataFrame(do.call(rbind, lapply(x, .spectrum_header)))
##     if (!missing(spectraData))
##         spectraData <- .combine_featureData(spectraData, fdata)
##     else spectraData <- fdata
##     if (all(is.na(spectraData$fileIdx))) {
##         spectraData$fileIdx <- 1L
##         for (i in seq_along(x))
##             x[[i]]@fromFile <- 1L
##     }
##     file <- unique(spectraData$fileIdx)
##     if (!is.null(names(x)))
##         rownames(spectraData) <- names(x)
##     else {
##         ## Use acquisitionNum to order the spectra per file. If that is NA
##         ## spectra get named according to their
##         from_file <- factor(spectraData$fileIdx, unique(spectraData$fileIdx))
##         names(x) <- formatFileSpectrumNames(
##             fileIds = spectraData$fileIdx,
##             spectrumIds = unsplit(lapply(
##                 split(spectraData$acquisitionNum, from_file),
##                 order), from_file),
##             nFiles = length(file))
##         rownames(spectraData) <- names(x)
##     }
##     if (missing(sampleData))
##         sampleData <- DataFrame(sampleIdx = unique(fdata$fileIdx))
##     backend <- BackendMemory()
##     backend <- backendInitialize(backend, as.character(file), spectraData, ...)
##     backend <- backendWriteSpectra(backend, x, spectraData)
##     new("Spectra",
##         backend = backend,
##         sampleData = sampleData,
##         spectraData = spectraData,
##         processingQueue = list(),
##         processing = paste0("Data loaded [", date(), "]")
##         )
## }

#' @description
#'
#' Combine two featureData `DataFrame`s. The resulting `DataFrame` contains all
#' columns from both `x` and `y`. For columns present in both `DataFrame`s those
#' in `x` will be used. Also, the resulting `DataFrame` uses the row names of
#' `x` unless `x` has no row names.
#'
#' @param x `DataFrame`
#'
#' @param y `DataFrame`
#'
#' @return `DataFrame` with the merged columns.
#'
#' @author Johannes Rainer
#'
#' @importMethodsFrom S4Vectors nrow rownames colnames cbind
#'
#' @noRd
.combine_featureData <- function(x, y) {
    if (nrow(x) != nrow(y))
        stop("'x' and 'y' have to have the same number of rows")
    if (is.null(rownames(x)) & !is.null(rownames(y)))
        rownames(x) <- rownames(y)
    cols_y <- !(colnames(y) %in% colnames(x))
    if (any(cols_y))
        x <- cbind(x, y[, cols_y, drop = FALSE])
    x
}

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
            nrow(object@spectraData), " spectra:\n", sep = "")
        if (nrow(object@spectraData)) {
            txt <- capture.output(
                object@spectraData[, c("msLevel", "retentionTime", "totIonCurrent")])
            cat(txt[-1], sep = "\n")
        }
        show(object@backend)
        if (length(object@processingQueue))
            cat("Lazy evaluation queue:", length(object@processingQueue),
                "processing step(s)\n")
        cat("Processing:\n", paste(object@processing, collapse="\n"), "\n")
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
