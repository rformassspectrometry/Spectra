#' @include hidden_aliases.R
NULL

################################################################################
##
## Spectra class, creation, data representation, export
##
################################################################################

#' @title The Spectra class to manage and access MS data
#'
#' @name Spectra-class
#'
#' @description
#'
#' The `Spectra` class encapsules spectral mass spectrometry (MS) data and
#' related metadata. The MS data is represented by a *backend* extending the
#' virual [MsBackend] class which provides the data to the `Spectra` object.
#' The `Spectra` class implements only data accessor, filtering and analysis
#' methods for the MS data and relies on its *backend* to provide the MS data.
#' This allows to change data representations of a `Spectra` object depending
#' on the user's needs and properties of the data. Different backends and
#' their properties are explained in the [MsBackend] documentation.
#'
#' Documentation on other topics and functionality of `Spectra`can be found in:
#'
#' LLLLLLL add links to individual documentations.
#' - [processingChunkSize()] for information on parallel and chunk-wise data
#'   processing.
#' - [plotSpectra()] for visualization of `Spectra`.
#'
#'
#' @section Creation of objects:
#'
#' `Spectra` classes can be created with the `Spectra()` constructor function
#' which supports the following formats:
#'
#' - parameter `object` is a `data.frame` or `DataFrame` containing the
#'   full spectrum data (spectra variables in columns as well as columns
#'   with the individual MS peak data, *m/z* and intensity). The provided
#'   `backend` (by default a [MsBackendMemory-class]) will be initialized
#'   with that data.
#'
#' - parameter `object` is a [MsBackend-class] (assumed to be already
#'   initialized).
#'
#' - parameter `object` is missing, in which case it is supposed that the data
#'   is provided by the [MsBackend-class] class passed along with the `backend`
#'   argument.
#'
#' - parameter `object` is of type `character` and is expected to be the file
#'   names(s) from which spectra should be imported. Parameter `source` allows
#'   to define a [MsBackend-class] that is able to import the data from the
#'   provided source files. The default value for `source` is [MsBackendMzR()]
#'   which allows to import spectra data from mzML, mzXML or CDF files.
#'
#' With `...` additional arguments can be passed to the backend's
#' [backendInitialize()] method. Parameter `backend` allows to specify which
#' [MsBackend-class] should be used for data representation and storage.
#'
#'
#' @section Data representation of a `Spectra`:
#'
#' The MS data which can be accessed through the `Spectra` object is
#' *represented* by its backend, which means that this backend defines how
#' and where the data is stored (e.g. in memory or on disk). The `Specrta`
#' object relies on the backend to provide the MS data whenever it needs it
#' for data processing.
#' Different backends with different properties, such as minimal memory
#' requirement or fast data access, are defined in the *Spectra* package or
#' one of the MsBackend* packages. More information on backends and their
#' properties is provided in the documentation of [MsBackend].
#'
#' On-disk backends keep only a limited amount of data in memory retrieving
#' most of the data (usually the MS peak data) upon request on-the-fly from
#' their on-disk data representations. Moving the on-disk data storage of such
#' a backend or a serialized object to a different location in the file
#' system will cause data corruption. The `dataStorageBasePath()` and
#' `dataStorageBasePath<-` functions allow in such cases (and if thebackend
#' classes support this operation), to get or change the *base*
#' path to the directory of the backend's data storage. In-memory backends
#' such as [MsBackendMemory] or [MsBackendDataFrame] keeping all MS data in
#' memory don't support, and need, this function, but for [MsBackendMzR] this
#' function can be used to update/adapt the path to the directory containing
#' the original data files. Thus, for `Spectra` objects (using this backend)
#' that were moved to another file system or computer, these functions allow to
#' adjust/adapt the base file path.
#'
#'
#' @section Changing data representation of a `Spectra`:
#'
#' The data representation, i.e. the backend of a `Spectra` object can be
#' changed with the `setBackend()` method that takes an instance of the new
#' backend as second parameter `backend`. A call to
#' `setBackend(sps, backend = MsBackendDataFrame())`
#' would for example change the backend of `sps` to the *in-memory*
#' `MsBackendDataFrame`. Changing to a backend is only supported if that
#' backend has a `data` parameter in its `backendInitialize()` method and if
#' `supportsSetBackend()` returns `TRUE` for that backend. `setBackend()` will
#' transfer the full spectra data from the originating backend as a `DataFrame`
#' to the new backend.
#'
#' Generally, it is not possible to change **to** a read-only backend such as
#' the [MsBackendMzR()] backend.
#'
#' The definition of the function is:
#' `setBackend(object, backend, ..., f = dataStorage(object),
#'     BPPARAM = bpparam())` and its parameters are:
#'
#' - `object`: the `Spectra` object.
#'
#' - `backend`: an instance of the new backend, e.g. `[MsBackendMemory()]`.
#'
#' - `f`: factor allowing to parallelize the change of the backends. By
#'   default the process of copying the spectra data from the original to the
#'   new backend is performed separately (and in parallel) for each file. Users
#'   are advised to use the default setting.
#'
#' - `...`: optional additional arguments passed to the [backendInitialize()]
#'   method of the new `backend`.
#'
#' - `BPPARAM`: setup for the parallel processing. See [bpparam()] for
#'   details.
#'
#'
#' @section Exporting data from a `Spectra` object:
#'
#' Data from a `Spectra` object can be **exported** to a file with the
#' `export()` function. The actual export of the data is performed by
#' the `export` method of the [MsBackend] class defined with the mandatory
#' parameter `backend` which defines also the format in which the data
#' is exported. Note however that not all backend classes support
#' export of data. From the `MsBackend` classes in the `Spectra` package
#' currently only the `MsBackendMzR` backend supports data export (to
#' mzML/mzXML file(s)); see the help page of the [MsBackend-class] for
#' information on its arguments or the examples below or the vignette
#' for examples.
#'
#' The definition of the function is
#' `export(object, backend,  ...)` and its
#' parameters are:
#'
#' - `object`: the `Spectra` object to be exported.
#'
#' - `backend`: instance of a class extending [MsBackend] which supports export
#'   of the data (i.e. which has a defined `export` method).
#'
#' - `...`: additional parameters specific for the `MsBackend` passed with
#'   parameter `backend`.
#'
#'
#' @details
#'
#' The `Spectra` class uses by default a lazy data manipulation strategy,
#' i.e. data manipulations such as performed with `replaceIntensitiesBelow()`
#' are not applied immediately to the data, but applied on-the-fly to the
#' spectrum data once it is retrieved. For some backends that allow to write
#' data back to the data storage (such as the [MsBackendMemory()],
#' [MsBackendDataFrame()] and [MsBackendHdf5Peaks()]) it is possible to apply
#' to queue with the `applyProcessing` function. See the *Data manipulation and
#' analysis *methods* section below for more details.
#'
#' Clarifications regarding scan/acquisition numbers and indices:
#'
#' - A `spectrumId` (or `spectrumID`) is a vendor specific field in
#'   the mzML file that contains some information about the
#'   run/spectrum, e.g.: `controllerType=0 controllerNumber=1
#'   scan=5281 file=2`
#'
#' - `acquisitionNum` is a more a less sanitize spectrum id generated
#'   from the `spectrumId` field by `mzR` (see
#'   [here](https://github.com/sneumann/mzR/blob/master/src/pwiz/data/msdata/MSData.cpp#L552-L580)).
#'
#' - `scanIndex` is the `mzR` generated sequence number of the
#'   spectrum in the raw file (which doesn't have to be the same as
#'   the `acquisitionNum`)
#'
#' See also [this issue](https://github.com/lgatto/MSnbase/issues/525).
#'
#' @md
#'
#' @author Sebastian Gibb, Johannes Rainer, Laurent Gatto, Philippine Louail
#'
#' @exportClass Spectra
#'
#' @exportMethod Spectra
#'
#' @examples
#'
#' ## ---- CREATION OF SPECTRA OBJECTS ----
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
#' ## Create a Spectra from mzML files and use the `MsBackendMzR` on-disk
#' ## backend.
#' sciex_file <- dir(system.file("sciex", package = "msdata"),
#'     full.names = TRUE)
#' sciex <- Spectra(sciex_file, backend = MsBackendMzR())
#' sciex
#'
#'
#' ## ---- CHANGING DATA REPRESENTATIONS ----
#'
#' ## The MS data is on disk and will be read into memory on-demand. We can
#' ## however change the backend to a MsBackendMemory backend which will
#' ## keep all of the data in memory.
#' sciex_im <- setBackend(sciex, MsBackendMemory())
#' sciex_im
#'
#' ## The `MsBackendMemory()` supports the `setBackend()` method:
#' supportsSetBackend(MsBackendMemory())
#'
#' ## Thus, it is possible to change to that backend with `setBackend()`. Most
#' ## read-only backends however don't support that, such as the
#' ## `MsBackendMzR` and `setBackend()` would fail to change to that backend.
#' supportsSetBackend(MsBackendMzR())
#'
#' ## The on-disk object `sciex` is light-weight, because it does not keep the
#' ## MS peak data in memory. The `sciex_im` object in contrast keeps all the
#' ## data in memory and its size is thus much larger.
#' object.size(sciex)
#' object.size(sciex_im)
#'
#' ## The spectra variable `dataStorage` returns for each spectrum the location
#' ## where the data is stored. For in-memory objects:
#' head(dataStorage(sciex_im))
#'
#' ## While objects that use an on-disk backend will list the files where the
#' ## data is stored.
#' head(dataStorage(sciex))
#'
#' ## The spectra variable `dataOrigin` returns for each spectrum the *origin*
#' ## of the data. If the data is read from e.g. mzML files, this will be the
#' ## original mzML file name:
#' head(dataOrigin(sciex))
#' head(dataOrigin(sciex_im))
#'
#'
#' ## ---- DATA EXPORT ----
#'
#' ## Some `MsBackend` classes provide an `export()` method to export the data
#' ## to the file format supported by the backend.
#' ## The `MsBackendMzR` for example allows to export MS data to mzML or
#' ## mzXML file(s), the `MsBackendMgf` (defined in the MsBackendMgf R package)
#' ## would allow to export the data in mgf file format.
#' ## Below we export the MS data in `data`. We call the `export()` method on
#' ## this object, specify the backend that should be used to export the data
#' ## (and which also defines the output format) and provide a file name.
#' fl <- tempfile()
#' export(data, MsBackendMzR(), file = fl)
#'
#' ## This exported our data in mzML format. Below we read the first 6 lines
#' ## from that file.
#' readLines(fl, n = 6)
#'
#' ## If only a single file name is provided, all spectra are exported to that
#' ## file. To export data with the `MsBackendMzR` backend to different files, a
#' ## file name for each individual spectrum has to be provided.
#' ## Below we export each spectrum to its own file.
#' fls <- c(tempfile(), tempfile())
#' export(data, MsBackendMzR(), file = fls)
#'
#' ## Reading the data from the first file
#' res <- Spectra(backendInitialize(MsBackendMzR(), fls[1]))
#'
#' mz(res)
#' mz(data)

#' The Spectra class
#'
#' The [Spectra-class] encapsulates data and meta-data for mass
#' spectrometry experiments.
#'
#' @slot backend A derivate of [MsBackend-class] holding/controlling the spectra
#'     data.
#'
#' @slot processingQueue `list` of `ProcessingStep` objects.
#'
#' @slot processingQueueVariables `character` of spectraVariables that should
#'     be passed to the processing step function.
#'
#' @slot processing A `character` storing logging information.
#'
#' @slot metadata A `list` storing experiment metadata.
#'
#' @slot version A `character(1)` containing the class version.
#'
#' @docType class
#'
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
        processingQueueVariables = "character",
        ## logging
        processing = "character",
        ## metadata
        metadata = "list",
        processingChunkSize = "numeric",
        version = "character"
    ),
    prototype = prototype(version = "0.3",
                          processingChunkSize = Inf)
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
        lp <- length(object@processing)
        if (lp) {
            lps <- object@processing
            if (lp > 3) {
                lps <- lps[1:3]
                lps <- c(lps, paste0("...", lp - 3, " more processings. ",
                                     "Use 'processingLog' to list all."))
            }
            cat("Processing:\n", paste(lps, collapse="\n "), "\n")
        }
    })

#' @rdname Spectra-class
setMethod("Spectra", "missing", function(object, processingQueue = list(),
                                         metadata = list(), ...,
                                         backend = MsBackendMemory(),
                                         BPPARAM = bpparam()) {
    new("Spectra", metadata = metadata, processingQueue = processingQueue,
        backend = backend)
})

#' @rdname Spectra-class
setMethod("Spectra", "MsBackend", function(object, processingQueue = list(),
                                           metadata = list(), ...,
                                           BPPARAM = bpparam()) {
    new("Spectra", metadata = metadata, processingQueue = processingQueue,
        backend = object)
})

#' @rdname Spectra-class
#'
#' @importFrom methods callNextMethod
setMethod("Spectra", "character", function(object, processingQueue = list(),
                                           metadata = list(),
                                           source = MsBackendMzR(),
                                           backend = source,
                                           ..., BPPARAM = bpparam()) {
    if (!length(object))
        Spectra(backend, metadata = metadata,
                processingQueue = processingQueue)
    else
        callNextMethod(object = object, processingQueue = processingQueue,
                       metadata = metadata, source = source, backend = backend,
                       ..., BPPARAM = BPPARAM)
})

#' @rdname Spectra-class
setMethod("Spectra", "ANY", function(object, processingQueue = list(),
                                     metadata = list(),
                                     source = MsBackendMemory(),
                                     backend = source,
                                     ..., BPPARAM = bpparam()) {
    sp <- new("Spectra", metadata = metadata, processingQueue = processingQueue,
              backend = backendInitialize(
                  source, object, ...,
                  BPPARAM = backendBpparam(source, BPPARAM)))
    if (class(source)[1L] != class(backend)[1L])
        setBackend(sp, backend, ..., BPPARAM = backendBpparam(backend, BPPARAM))
    else sp
})

#' @rdname Spectra-class
#'
#' @importMethodsFrom ProtGenerics setBackend
#'
#' @exportMethod setBackend
setMethod(
    "setBackend", c("Spectra", "MsBackend"),
    function(object, backend, f = processingChunkFactor(object), ...,
             BPPARAM = bpparam()) {
        backend_class <- class(object@backend)[1L]
        BPPARAM <- backendBpparam(object@backend, BPPARAM)
        BPPARAM <- backendBpparam(backend, BPPARAM)
        if (!supportsSetBackend(backend))
            stop(class(backend), " does not support 'setBackend'")
        if (!length(object)) {
            bknds <- backendInitialize(
                backend, data = spectraData(object@backend), ...)
        } else {
            if (!is.factor(f))
                f <- force(factor(f, levels = unique(f)))
            if (length(f) && (length(levels(f)) > 1)) {
                if (length(f) != length(object))
                    stop("length of 'f' has to match the length of 'object'")
                bknds <- bplapply(
                    split(object@backend, f = f),
                    function(z, ...) {
                        backendInitialize(backend,
                                          data = spectraData(z), ...,
                                          BPPARAM = SerialParam())
                    }, ..., BPPARAM = BPPARAM)
                bknds <- backendMerge(bknds)
                ## That below ensures the backend is returned in its original
                ## order - unsplit does unfortunately not work.
                if (is.unsorted(f))
                    bknds <- bknds[order(unlist(split(seq_along(bknds), f),
                                                use.names = FALSE))]
            } else {
                bknds <- backendInitialize(
                    backend, data = spectraData(object@backend), ...)
            }
        }
        object@backend <- bknds
        object@processing <- .logging(object@processing,
                                      "Switch backend from ",
                                      backend_class, " to ",
                                      class(object@backend))
        object
    })

#' @rdname Spectra-class
#'
#' @export
setMethod("export", "Spectra",
          function(object, backend, ...) {
              if (missing(backend))
                  stop("Parameter 'backend' is required.")
              export(backend, object, ...)
          })

#' @rdname Spectra-class
setMethod("dataStorageBasePath", "Spectra", function(object) {
    dataStorageBasePath(object@backend)
})

#' @rdname Spectra-class
setReplaceMethod("dataStorageBasePath", "Spectra", function(object, value) {
    dataStorageBasePath(object@backend) <- value
    object
}
)

## CONTINUNE HERE:
## - check if some additional methods/functions need to be moved up.

################################################################################
##
## Merging, splitting and aggregating Spectra: length of Spectra is changed
##
################################################################################

#' @title Merging, splitting and aggregating Spectra
#'
#' @aliases [,Spectra-method

################################################################################
##
## Filtering, subsetting Spectra: subsetting Spectra and its data content.
##
################################################################################

#' @title Filtering and subsetting Spectra objects
#'

################################################################################
##
## Accessing and adding/setting/changing MS data.
##
################################################################################

#' @title Accessing mass spectrometry data
#'
#'

################################################################################
##
## Data manipulation and analysis operations (lazy processing)
##
################################################################################

#' @title Data manipulation and analysis methods
#'

################################################################################
##
## Spectra similarity calculations
##
################################################################################

#' @title Spectra similarity calculations


#' @rdname Spectra
#'
#' @importFrom MsCoreUtils vapply1c
#'
#' @exportMethod c
setMethod("c", "Spectra", function(x, ...) {
    .concatenate_spectra(unname(list(unname(x), ...)))
})

#' @rdname Spectra
setMethod("split", "Spectra", function(x, f, drop = FALSE, ...) {
    bcknds <- split(x@backend, f, ...)
    lapply(bcknds, function(b) {
        slot(x, "backend", check = FALSE) <- b
        x
    })
})


#### ---------------------------------------------------------------------------
##
##                          ACCESSOR METHODS
##
#### ---------------------------------------------------------------------------

#' @rdname Spectra
setMethod("acquisitionNum", "Spectra", function(object)
    acquisitionNum(object@backend))

#' @rdname Spectra
setMethod(
    "peaksData", "Spectra",
    function(object, columns = c("mz", "intensity"),
             f = processingChunkFactor(object), ..., BPPARAM = bpparam()) {
        if (length(object@processingQueue) || length(f))
            SimpleList(.peaksapply(object, columns = columns, f = f))
        else SimpleList(peaksData(object@backend, columns = columns))
    })

#' @rdname Spectra
setMethod("peaksVariables", "Spectra", function(object)
    peaksVariables(object@backend))

#' @importFrom methods setAs
setAs("Spectra", "list", function(from, to) {
    .peaksapply(from)
})

setAs("Spectra", "SimpleList", function(from, to) {
    peaksData(from)
})

#' @rdname Spectra
setMethod("centroided", "Spectra", function(object) {
    centroided(object@backend)
})

#' @rdname Spectra
setReplaceMethod("centroided", "Spectra", function(object, value) {
    centroided(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("collisionEnergy", "Spectra", function(object) {
    collisionEnergy(object@backend)
})

#' @rdname Spectra
setReplaceMethod("collisionEnergy", "Spectra", function(object, value) {
    collisionEnergy(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("dataOrigin", "Spectra", function(object) dataOrigin(object@backend))

#' @rdname Spectra
setReplaceMethod("dataOrigin", "Spectra", function(object, value) {
    dataOrigin(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("dataStorage", "Spectra",
          function(object) dataStorage(object@backend))

#' @rdname Spectra
setMethod("dropNaSpectraVariables", "Spectra", function(object) {
    object@backend <- dropNaSpectraVariables(object@backend)
    object
})

#' @rdname Spectra
setMethod("intensity", "Spectra", function(object,
                                           f = processingChunkFactor(object),
                                           ...) {
    if (length(object@processingQueue) || length(f))
        NumericList(.peaksapply(object, FUN = function(z, ...) z[, 2],
                                f = f, ...), compress = FALSE)
    else intensity(object@backend)
})

#' @rdname Spectra
setMethod("ionCount", "Spectra", function(object) {
    if (length(object))
        unlist(.peaksapply(
            object, FUN = function(pks, ...) sum(pks[, 2], na.rm = TRUE)),
            use.names = FALSE)
    else numeric()
})

#' @rdname Spectra
setMethod("isCentroided", "Spectra", function(object, ...) {
    if (length(object))
        unlist(.peaksapply(object, FUN = .peaks_is_centroided),
               use.names = FALSE)
    else logical()
})

#' @rdname Spectra
setMethod("isEmpty", "Spectra", function(x) {
    if (length(x))
        unlist(.peaksapply(x, FUN = function(pks, ...) nrow(pks) == 0),
               use.names = FALSE)
    else logical()
})

#' @rdname Spectra
setMethod("isolationWindowLowerMz", "Spectra", function(object) {
    isolationWindowLowerMz(object@backend)
})

#' @rdname Spectra
setReplaceMethod("isolationWindowLowerMz", "Spectra", function(object, value) {
    isolationWindowLowerMz(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("isolationWindowTargetMz", "Spectra", function(object) {
    isolationWindowTargetMz(object@backend)
})

#' @rdname Spectra
setReplaceMethod("isolationWindowTargetMz", "Spectra", function(object, value) {
    isolationWindowTargetMz(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("isolationWindowUpperMz", "Spectra", function(object) {
    isolationWindowUpperMz(object@backend)
})

#' @rdname Spectra
setReplaceMethod("isolationWindowUpperMz", "Spectra", function(object, value) {
    isolationWindowUpperMz(object@backend) <- value
    object
})

#' @rdname Spectra
#'
#' @exportMethod containsMz
setMethod("containsMz", "Spectra", function(object, mz = numeric(),
                                            tolerance = 0,
                                            ppm = 20, which = c("any", "all"),
                                            BPPARAM = bpparam()) {
    cond_fun <- match.fun(match.arg(which))
    if (all(is.na(mz)))
        return(rep(NA, length(object)))
    mz <- unique(sort(mz))
    BPPARAM <- backendBpparam(object@backend, BPPARAM)
    ## TODO: fix to use .peaksapply instead.
    if (is(BPPARAM, "SerialParam"))
        .has_mz(object, mz, tolerance = tolerance, ppm = ppm,
                condFun = cond_fun, parallel = BPPARAM)
    else {
        sp <- SerialParam()
        f <- as.factor(dataStorage(object))
        res <- .lapply(object, FUN = .has_mz, mz = mz, tolerance = tolerance,
                       condFun = cond_fun, parallel = sp, f = f,
                       BPPARAM = BPPARAM)
        unsplit(res, f = f)
    }
})

#' @rdname Spectra
#'
#' @exportMethod containsNeutralLoss
setMethod("containsNeutralLoss", "Spectra", function(object, neutralLoss = 0,
                                                     tolerance = 0, ppm = 20,
                                                     BPPARAM = bpparam()) {
    BPPARAM <- backendBpparam(object@backend, BPPARAM)
    ## TODO: FIX me to use chunk size.
    if (is(BPPARAM, "SerialParam")) {
        .has_mz_each(object, precursorMz(object) - neutralLoss,
                     tolerance = tolerance, ppm = ppm, parallel = BPPARAM)
    } else {
        sp <- SerialParam()
        f <- as.factor(dataStorage(object))
        res <- .lapply(object, FUN = function(obj, n, tol, ppm, par) {
            .has_mz_each(obj, precursorMz(obj) - n, tolerance = tol,
                         ppm = ppm, parallel = sp)
        }, n = neutralLoss, tol = tolerance, ppm = ppm, par = sp, f = f,
                       BPPARAM = BPPARAM)
        unsplit(res, f = f)
    }
})

#' @rdname Spectra
#'
#' @importMethodsFrom ProtGenerics spectrapply
#'
#' @exportMethod spectrapply
setMethod("spectrapply", "Spectra", function(object, FUN, ...,
                                             chunkSize = integer(),
                                             f = factor(),
                                             BPPARAM = SerialParam()) {
    if (missing(FUN))
        FUN <- identity
    if (length(chunkSize))
        return(chunkapply(object, FUN, ..., chunkSize = chunkSize))
    if (!length(f))
        f <- as.factor(seq_along(object))
    .lapply(object, FUN = FUN, f = f, ...,
            BPPARAM = backendBpparam(object@backend, BPPARAM))
})

#' @rdname Spectra
#'
#' @exportMethod length
setMethod("length", "Spectra", function(x) length(x@backend))

#' @rdname Spectra
setMethod("msLevel", "Spectra", function(object) msLevel(object@backend))

#' @rdname Spectra
setMethod("mz", "Spectra", function(object, f = processingChunkFactor(object),
                                    ...) {
    if (length(object@processingQueue) || length(f))
        NumericList(.peaksapply(object, FUN = function(z, ...) z[, 1],
                                f = f, ...), compress = FALSE)
    else mz(object@backend)
})

#' @rdname Spectra
#'
#' @exportMethod lengths
setMethod("lengths", "Spectra", function(x, use.names = FALSE) {
    f <- .parallel_processing_factor(x)
    if (length(x)) {
        if (length(x@processingQueue) || length(f))
            unlist(.peaksapply(x, FUN = function(pks, ...) nrow(pks)),
                   use.names = use.names)
        else lengths(x@backend, use.names = use.names)
    } else integer()
})

#' @rdname Spectra
setMethod("polarity", "Spectra", function(object) {
    polarity(object@backend)
})

#' @rdname Spectra
setReplaceMethod("polarity", "Spectra", function(object, value) {
    polarity(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("precScanNum", "Spectra", function(object) {
    precScanNum(object@backend)
})

#' @rdname Spectra
setMethod("precursorCharge", "Spectra", function(object) {
    precursorCharge(object@backend)
})

#' @rdname Spectra
setMethod("precursorIntensity", "Spectra", function(object) {
    precursorIntensity(object@backend)
})

#' @rdname Spectra
setMethod("precursorMz", "Spectra", function(object) {
    precursorMz(object@backend)
})

#' @rdname Spectra
setMethod("rtime", "Spectra", function(object) {
    rtime(object@backend)
})

#' @rdname Spectra
setReplaceMethod("rtime", "Spectra", function(object, value) {
    rtime(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("scanIndex", "Spectra", function(object) {
    scanIndex(object@backend)
})

#' @rdname Spectra
setMethod(
    "selectSpectraVariables", "Spectra",
    function(object, spectraVariables = union(spectraVariables(object),
                                              peaksVariables(object))) {
        spectraVariables <- union(spectraVariables, "dataStorage")
        object@backend <- selectSpectraVariables(
            object@backend, spectraVariables = spectraVariables)
        object
    })

#' @rdname Spectra
setMethod("smoothed", "Spectra", function(object) {
    smoothed(object@backend)
})

#' @rdname Spectra
setReplaceMethod("smoothed", "Spectra", function(object, value) {
    smoothed(object@backend) <- value
    object
})

#' @rdname Spectra
#'
#' @importMethodsFrom ProtGenerics spectraData
#'
#' @exportMethod spectraData
setMethod(
    "spectraData", "Spectra",
    function(object, columns = spectraVariables(object)) {
        if (length(object@processingQueue) &&
            length(pcns <- intersect(columns, peaksVariables(object)))) {
            ## If user requests peaks variables we need to ensure that the
            ## processing queue is executed.
            scns <- setdiff(columns, pcns)
            if (length(scns))
                spd <- spectraData(object@backend, columns = scns)
            else
                spd <- make_zero_col_DFrame(nrow = length(object))
            pkd <- peaksData(object, columns = pcns)
            ## Add individual peaks variables to the `DataFrame`.
            for (pcn in pcns) {
                vals <- lapply(pkd, `[`, , pcn)
                if (pcn %in% c("mz", "intensity"))
                    vals <- NumericList(vals, compress = FALSE)
                spd <- do.call(`[[<-`, list(spd, i = pcn, value = vals))
            }
            spd
        } else
            spectraData(object@backend, columns = columns)
    })

#' @rdname Spectra
#'
#' @importMethodsFrom ProtGenerics spectraData<-
#'
#' @exportMethod spectraData<-
setReplaceMethod("spectraData", "Spectra", function(object, value) {
    if (!inherits(value, "DataFrame"))
        stop("'spectraData<-' expects a 'DataFrame' as input.", call. = FALSE)
    pvs <- peaksVariables(object)
    if (length(object@processingQueue) &&
        any(colnames(value) %in% pvs))
        stop("Can not replace peaks variables with a non-empty processing ",
             "queue. Please use 'object <- applyProcessing(object)' to apply ",
             "and clear the processing queue. Note that 'applyProcessing' ",
             "requires a *writeable* backend. Use e.g. 'object <- ",
             "setBackend(object, MsBackendMemory())' if needed.")
    pvs <- setdiff(pvs, colnames(value))
    if (length(pvs)) {
        sd <- spectraData(object, pvs)
        for (pv in pvs) {
            value <- do.call("$<-", list(value, name = pv, sd[, pv]))
        }
        object@processingQueue <- list()
    }
    spectraData(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("spectraNames", "Spectra", function(object) {
    spectraNames(object@backend)
})

#' @rdname Spectra
setReplaceMethod("spectraNames", "Spectra", function(object, value) {
    spectraNames(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("spectraVariables", "Spectra", function(object) {
    setdiff(spectraVariables(object@backend), peaksVariables(object@backend))
})

#' @rdname Spectra
setMethod("tic", "Spectra", function(object, initial = TRUE) {
    if (!length(object))
        return(numeric())
    if (initial)
        tic(object@backend, initial = initial)
    else ionCount(object)
})

#' @rdname Spectra
#'
#' @importMethodsFrom S4Vectors $
#'
#' @export
setMethod("$", "Spectra", function(x, name) {
    if (!(name %in% c(spectraVariables(x@backend), peaksVariables(x@backend))))
        stop("No spectra variable '", name, "' available")
    if (name == "mz")
        mz(x)
    else if (name == "intensity")
        intensity(x)
    else {
        if (length(x@processingQueue) && name %in% peaksVariables(x))
            .peaksapply(x, FUN = function(z, ...) z[, name],
                        columns = c("mz", "intensity", name))
        else
            do.call("$", list(x@backend, name))
    }
})

#' @rdname Spectra
#'
#' @export
setReplaceMethod("$", "Spectra", function(x, name, value) {
    if (length(x@processingQueue) &&
        any(name %in% peaksVariables(x)))
        stop("Can not replace peaks variables with a non-empty processing ",
             "queue. Please use 'object <- applyProcessing(object)' to apply ",
             "and clear the processing queue. Note that 'applyProcessing' ",
             "requires a *writeable* backend. Use e.g. 'object <- ",
             "setBackend(object, MsBackendMemory())' if needed.")
    x@backend <- do.call("$<-", list(x@backend, name, value))
    x
})

#' @rdname Spectra
#'
#' @export
setMethod("[[", "Spectra", function(x, i, j, ...) {
    if (!is.character(i))
        stop("'i' is supposed to be a character defining the spectra ",
             "variable to access.")
    if (!missing(j))
        stop("'j' is not supported.")
    if (!(i %in% c(spectraVariables(x), "mz", "intensity")))
        stop("No spectra variable '", i, "' available")
    if (i == "mz")
        mz(x)
    else if (i == "intensity")
        intensity(x)
    else
        do.call("[[", list(x@backend, i))
})

#' @rdname Spectra
#'
#' @export
setReplaceMethod("[[", "Spectra", function(x, i, j, ..., value) {
    if (!is.character(i))
        stop("'i' is supposed to be a character defining the spectra ",
             "variable to replace or create.")
    if (!missing(j))
        stop("'j' is not supported.")
    x@backend <- do.call("[[<-", list(x@backend, i = i, value = value))
    x
})

#### ---------------------------------------------------------------------------
##
##                      FILTERING AND SUBSETTING
##
#### ---------------------------------------------------------------------------

#' @rdname Spectra
setMethod("[", "Spectra", function(x, i, j, ..., drop = FALSE) {
    if (!missing(j))
        stop("Subsetting 'Spectra' by columns is not (yet) supported")
    if (missing(i))
        return(x)
    slot(x, "backend", check = FALSE) <- x@backend[i = i]
    x
})

#' @rdname Spectra
setMethod("filterAcquisitionNum", "Spectra", function(object, n = integer(),
                                                      dataStorage = character(),
                                                      dataOrigin = character()) {
    if (length(dataStorage) && !is.character(dataStorage))
        stop("'dataStorage' is expected to be of type character")
    if (length(dataOrigin) && !is.character(dataOrigin))
        stop("'dataOrigin' is expected to be of type character")
    object@backend <- filterAcquisitionNum(object@backend, n,
                                           dataStorage, dataOrigin)
    object@processing <- .logging(object@processing,
                                  "Filter: select by: ", length(n),
                                  " acquisition number(s) in ",
                                  max(length(dataStorage), length(dataOrigin)),
                                  " file(s)")
    object
})

#' @rdname Spectra
setMethod("filterEmptySpectra", "Spectra", function(object) {
    object@backend <- object@backend[as.logical(lengths(object))]
    object@processing <- .logging(object@processing,
                                  "Filter: removed empty spectra.")
    object
})

#' @rdname Spectra
setMethod("filterDataOrigin", "Spectra", function(object,
                                                  dataOrigin = character()) {
    if (length(dataOrigin) && !is.character(dataOrigin))
        stop("'dataOrigin' is expected to be of type character")
    object@backend <- filterDataOrigin(object@backend, dataOrigin = dataOrigin)
    object@processing <- .logging(object@processing,
                                  "Filter: select data origin(s) ",
                                  paste0(dataOrigin, collapse = ", "))
    object
})

#' @rdname Spectra
setMethod("filterDataStorage", "Spectra", function(object,
                                                   dataStorage = character()) {
    if (length(dataStorage) && !is.character(dataStorage))
        stop("'dataStorage' is expected to be of type character")
    object@backend <- filterDataStorage(object@backend, dataStorage)
    object@processing <- .logging(object@processing,
                                  "Filter: select data storage(s) ",
                                  paste0(dataStorage, collapse = ", "))
    object
})

#' @rdname Spectra
#'
#' @exportMethod filterFourierTransformArtefacts
setMethod("filterFourierTransformArtefacts", "Spectra",
          function(object, halfWindowSize = 0.05, threshold = 0.2,
                   keepIsotopes = TRUE, maxCharge = 5,
                   isotopeTolerance = 0.005) {
              object <- addProcessing(object, .peaks_remove_fft_artifact,
                                      halfWindowSize = halfWindowSize,
                                      threshold = threshold,
                                      keepIsotopes = keepIsotopes,
                                      maxCharge = maxCharge,
                                      isotopeTolerance = isotopeTolerance)
              object@processing <- .logging(
                  object@processing, "Remove fast fourier artefacts.")
              object
          })

#' @rdname Spectra
#'
#' @importMethodsFrom ProtGenerics filterIntensity
#'
#' @exportMethod filterIntensity
setMethod("filterIntensity", "Spectra",
          function(object, intensity = c(0, Inf),
                   msLevel. = uniqueMsLevels(object), ...) {
              if (!.check_ms_level(object, msLevel.))
                  return(object)
              if (is.numeric(intensity)) {
                  if (length(intensity) == 1)
                      intensity <- c(intensity, Inf)
                  if (length(intensity) != 2)
                      stop("'intensity' should be of length specifying a ",
                           "lower intensity limit or of length two defining ",
                           "a lower and upper limit.")
                  object <- addProcessing(object, .peaks_filter_intensity,
                                          intensity = intensity,
                                          msLevel = msLevel.,
                                          spectraVariables = "msLevel")
                  object@processing <- .logging(
                      object@processing, "Remove peaks with intensities ",
                      "outside [", intensity[1], ", ", intensity[2],
                      "] in spectra of MS level(s) ",
                      paste0(msLevel., collapse = ", "), ".")
              } else {
                  if (is.function(intensity)) {
                      object <- addProcessing(
                          object, .peaks_filter_intensity_function,
                          intfun = intensity, msLevel = msLevel.,
                          args = list(...), spectraVariables = "msLevel")
                      object@processing <- .logging(
                          object@processing, "Remove peaks based on their ",
                          "intensities and a user-provided function ",
                          "in spectra of MS level(s) ",
                          paste0(msLevel., collapse = ", "), ".")
                  }
                  else stop("'intensity' has to be numeric or a function")
              }
              object
          })


#' @rdname Spectra
setMethod("filterIsolationWindow", "Spectra", function(object, mz = numeric()) {
    object@backend <- filterIsolationWindow(object@backend, mz = mz)
    object@processing <- .logging(object@processing,
                                  "Filter: select spectra containing m/z ",
                                  mz, " in their isolation window")
    object
})

#' @rdname Spectra
setMethod("filterMsLevel", "Spectra", function(object, msLevel. = integer()) {
    object@backend <- filterMsLevel(object@backend, msLevel = msLevel.)
    object@processing <- .logging(object@processing,
                                  "Filter: select MS level(s) ",
                                  paste0(unique(msLevel.), collapse = " "))
    object
})

#' @rdname Spectra
#'
#' @importMethodsFrom ProtGenerics filterMzRange
#'
#' @export
setMethod("filterMzRange", "Spectra",
          function(object, mz = numeric(), msLevel. = uniqueMsLevels(object),
                   keep = TRUE) {
              if (!.check_ms_level(object, msLevel.))
                  return(object)
              if (!length(mz)) mz <- c(-Inf, Inf)
              else mz <- range(mz)
              object <- addProcessing(object, .peaks_filter_mz_range, mz = mz,
                                      msLevel = msLevel., keep = keep,
                                      spectraVariables = "msLevel")
              if (keep) keep_or_remove <- "select"
              else keep_or_remove <- "remove"
              object@processing <- .logging(
                  object@processing, "Filter: ", keep_or_remove,
                  " peaks with an m/z within [", mz[1L], ", ", mz[2L], "]")
              object
          })

#' @rdname Spectra
#'
#' @importMethodsFrom ProtGenerics filterMzValues
#'
#' @export
setMethod("filterMzValues", "Spectra",
          function(object, mz = numeric(), tolerance = 0, ppm = 20,
                   msLevel. = uniqueMsLevels(object), keep = TRUE) {
              if (!.check_ms_level(object, msLevel.))
                  return(object)
              l <- length(mz)
              if (length(tolerance) != 1)
                  stop("'tolerance' should be of length 1")
              if (length(ppm) != 1)
                  stop("'ppm' should be of length 1")
              if (is.unsorted(mz)) {
                  idx <- order(mz)
                  mz <- mz[idx]
                  if (length(tolerance) == l)
                      tolerance <- tolerance[idx]
                  if (length(ppm) == l)
                      ppm <- ppm[idx]
              }
              object <- addProcessing(object, .peaks_filter_mz_value,
                                      mz = mz, tolerance = tolerance,
                                      ppm = ppm, msLevel = msLevel.,
                                      keep = keep, spectraVariables = "msLevel")
              if (length(mz) <= 3)
                  what <- paste0(format(mz, digits = 4), collapse = ", ")
              else what <- ""
              if (keep)
                  keep_or_remove <- "select"
              else keep_or_remove <- "remove"
              object@processing <- .logging(
                  object@processing, "Filter: ", keep_or_remove,
                  " peaks matching provided m/z values ", what)
              object
          })

#' @rdname Spectra
setMethod("filterPolarity", "Spectra", function(object, polarity = integer()) {
    object@backend <- filterPolarity(object@backend, polarity = polarity)
    object@processing <- .logging(object@processing,
                                  "Filter: select spectra with polarity ",
                                  paste0(polarity, collapse = " "))
    object
})

#' @rdname Spectra
#'
#' @export
setMethod("filterPrecursorMz", "Spectra",
          function(object, mz = numeric()) {
              .Deprecated(
                  msg = paste0("'filterPrecursorMz' is deprecated. Please use",
                               " 'filterPrecursorMzRange' instead."))
              object@backend <- filterPrecursorMzRange(object@backend, mz)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select spectra with a precursor m/z within [",
                  paste0(mz, collapse = ", "), "]")
              object
          })

#' @rdname Spectra
setMethod("filterPrecursorMzRange", "Spectra",
          function(object, mz = numeric()) {
              object@backend <- filterPrecursorMzRange(object@backend, mz)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select spectra with a precursor m/z within [",
                  paste0(mz, collapse = ", "), "]")
              object
          })

#' @rdname Spectra
setMethod("filterPrecursorMzValues", "Spectra",
          function(object, mz = numeric(), ppm = 20, tolerance = 0) {
              object@backend <- filterPrecursorMzValues(
                  object@backend, sort(mz), ppm = ppm, tolerance = tolerance)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select spectra with precursor m/z matching ",
                  paste0(mz, collapse = ", "), "")
              object
          })

#' @rdname Spectra
setMethod("filterPrecursorCharge", "Spectra",
          function(object, z = integer()) {
              z <- unique(z)
              object@backend <- filterPrecursorCharge(object@backend, z)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select spectra with a precursor charge ",
                  paste0(z, collapse = ", "))
              object
          })

#' @rdname Spectra
setMethod("filterPrecursorScan", "Spectra",
          function(object, acquisitionNum = integer(), f = dataOrigin(object)) {
              if (!all(f %in% unique(dataOrigin(object))))
                  stop("'f' must be in dataOrigin().")
              object@backend <- filterPrecursorScan(object@backend,
                                                    acquisitionNum,
                                                    f = dataOrigin(object))
              object@backend <- filterDataOrigin(object@backend, dataOrigin = f)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select parent/children scans for ",
                  paste0(acquisitionNum, collapse = " "))
              object
          })

#' @rdname Spectra
setMethod("filterRt", "Spectra",
          function(object, rt = numeric(), msLevel. = uniqueMsLevels(object)) {
              if (!is.numeric(msLevel.))
                  stop("Please provide a numeric MS level.")
              if (length(rt) != 2L || !is.numeric(rt) || rt[1] >= rt[2])
                  stop("Please provide a lower and upper numeric retention",
                       " time range.")
              if (length(rt))
                  rt <- range(rt)
              else rt <- c(-Inf, Inf)
              object@backend <- filterRt(object@backend, rt, msLevel.)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select retention time [", rt[1], "..", rt[2],
                  "] on MS level(s) ", paste0(msLevel., collapse = " "))
              object
          })

#' @rdname Spectra
setMethod("reset", "Spectra", function(object, ...) {
    object@backend <- reset(object@backend)
    object@processingQueue <- list()
    if (!.hasSlot(object, "processingQueueVariables"))
        object <- updateObject(object, check = FALSE)
    object@processingQueueVariables <- character()
    object@processing <- .logging(object@processing, "Reset object.")
    object
})

#' @rdname Spectra
setMethod("filterRanges", "Spectra",
          function(object, spectraVariables = character(), ranges = numeric(),
                   match = c("all", "any")){
              object@backend <- filterRanges(object@backend, spectraVariables,
                                             ranges, match)
              object@processing <- .logging(object@processing,
                                            "Filter: select spectra with a ",
                                            spectraVariables, " within: [",
                                            ranges[seq(ranges)%% 2 != 0], ", ",
                                            ranges[seq(ranges)%% 2 == 0], "]"
                                            )
              object
              })

#' @rdname Spectra
setMethod("filterValues", "Spectra",
          function(object, spectraVariables = character(), values = numeric(),
                   ppm = 0, tolerance = 0, match = c("all", "any")){
              object@backend <- filterValues(object@backend, spectraVariables,
                                             values, ppm, tolerance, match)
              object@processing <- .logging(object@processing,
                                            "Filter: select spectra with a ",
                                            spectraVariables, " similar to: ",
                                            values)
              object
              })

#### ---------------------------------------------------------------------------
##
##                      DATA MANIPULATION METHODS
##
#### ---------------------------------------------------------------------------

#' @rdname Spectra
#'
#' @importMethodsFrom ProtGenerics bin
#'
#' @exportMethod bin
setMethod("bin", "Spectra", function(x, binSize = 1L, breaks = NULL,
                                     msLevel. = uniqueMsLevels(x),
                                     FUN = sum, zero.rm = TRUE) {
    if (!.check_ms_level(x, msLevel.))
        return(x)
    if (!length(breaks)) {
        mzr <- range(.peaksapply(filterMsLevel(x, msLevel.),
                                 function(z, ...) z[c(1L, nrow(z))]
                                 ), na.rm = TRUE)
        breaks <- seq(floor(mzr[1]), ceiling(mzr[2]), by = binSize)
        breaks <- MsCoreUtils:::.fix_breaks(breaks, mzr)
    }
    mids <- (breaks[-length(breaks)] + breaks[-1L]) / 2
    x <- addProcessing(x, .peaks_bin, breaks = breaks, mids = mids,
                       agg_fun = FUN, msLevel = msLevel., zero.rm = zero.rm,
                       spectraVariables = "msLevel")
    x@processing <- .logging(x@processing,
                             "Spectra of MS level(s) ",
                             paste0(msLevel., collapse = ", "),
                             " binned.")
    x
})

#' @rdname Spectra
#'
#' @exportMethod compareSpectra
#'
#' @importFrom MsCoreUtils ndotproduct
#'
#' @importMethodsFrom ProtGenerics compareSpectra
#'
#' @exportMethod compareSpectra
setMethod("compareSpectra", signature(x = "Spectra", y = "Spectra"),
          function(x, y, MAPFUN = joinPeaks, tolerance = 0, ppm = 20,
                   FUN = ndotproduct, ..., SIMPLIFY = TRUE) {
              mat <- .compare_spectra_chunk(x, y, MAPFUN = MAPFUN,
                                            tolerance = tolerance,
                                            ppm = ppm, FUN = FUN, ...)
              if (SIMPLIFY && (length(x) == 1 || length(y) == 1))
                  mat <- as.vector(mat)
              mat
          })
#' @rdname Spectra
setMethod("compareSpectra", signature(x = "Spectra", y = "missing"),
          function(x, y = NULL, MAPFUN = joinPeaks, tolerance = 0, ppm = 20,
                   FUN = ndotproduct, ..., SIMPLIFY = TRUE) {
              if (length(x) == 1)
                  return(compareSpectra(x, x, MAPFUN = MAPFUN,
                                        tolerance = tolerance,
                                        ppm = ppm, FUN = FUN, ...,
                                        SIMPLIFY = SIMPLIFY))
              mat <- .compare_spectra_self(x, MAPFUN = MAPFUN, FUN = FUN,
                                           tolerance = tolerance, ppm = ppm,
                                           ...)
              if (SIMPLIFY && length(x) == 1)
                  mat <- as.vector(mat)
              mat
          })

## estimateMzResolution

## estimateNoise

## normalize

#' @rdname Spectra
#'
#' @exportMethod pickPeaks
setMethod("pickPeaks", "Spectra",
          function(object, halfWindowSize = 2L,
                   method = c("MAD", "SuperSmoother"), snr = 0, k = 0L,
                   descending = FALSE, threshold = 0,
                   msLevel. = uniqueMsLevels(object), ...) {
    if (!.check_ms_level(object, msLevel.))
        return(object)
    if (!is.integer(halfWindowSize) || length(halfWindowSize) != 1L ||
        halfWindowSize <= 0L)
        stop("Argument 'halfWindowSize' has to be an integer of length 1 ",
             "and > 0.")
    if (!is.numeric(snr) || length(snr) != 1L || snr < 0L)
        stop("Argument 'snr' has to be a numeric of length 1 that is >= 0.")
    if (!is.integer(k) || length(k) != 1L || k < 0L)
        stop("Argument 'k' has to be an integer of length 1 that is >= 0.")
    if (!is.logical(descending) || length(descending) != 1L ||
        is.na(descending))
        stop("Argument 'descending' has to be just TRUE or FALSE")
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        threshold < 0L || threshold > 1L)
        stop("Argument 'threshold' has to be a numeric of length 1 ",
             "that is >= 0 and <= 1.")

    method <- match.arg(method)

    object <- addProcessing(object, .peaks_pick,
                            halfWindowSize = halfWindowSize, method = method,
                            snr = snr, k = k, descending = descending,
                            threshold = threshold, msLevel = msLevel., ...,
                            spectraVariables = c("msLevel", "centroided"))
    object$centroided[msLevel(object) %in% msLevel.] <- TRUE
    object@processing <- .logging(object@processing,
                                  "Peak picking with ", method,
                                  " noise estimation, hws = ", halfWindowSize,
                                  ", snr = ", snr,
                                  if (k > 0) " and centroid refinement")
    object
})

## quantify

## removeReporters

#' @rdname Spectra
#'
#' @exportMethod replaceIntensitiesBelow
setMethod("replaceIntensitiesBelow", "Spectra",
          function(object, threshold = min, value = 0,
                   msLevel. = uniqueMsLevels(object)) {
              if (!is.numeric(threshold) && !is.function(threshold))
                  stop("Argument 'threshold' has to be either numeric or ",
                       "a function.")
              if (!.check_ms_level(object, msLevel.))
                  return(object)
              object <- addProcessing(
                  object, .peaks_replace_intensity, threshold = threshold,
                  value = value, msLevel = msLevel.,
                  spectraVariables = c("msLevel", "centroided"))
              msg <- ifelse(
                  is.function(threshold),
                  yes = "a threshold defined by a provided function",
                  no = threshold)
              object@processing <- .logging(object@processing,
                                            "Signal <= ", msg,
                                            " in MS level(s) ",
                                            paste0(msLevel., collapse = ", "),
                                            " set to 0")
              object
          })


#' @rdname Spectra
#'
#' @importFrom ProtGenerics smooth
#' @importFrom MsCoreUtils coefMA coefWMA coefSG
#' @exportMethod smooth
setMethod("smooth", "Spectra",
          function(x, halfWindowSize = 2L,
                   method = c("MovingAverage", "WeightedMovingAverage",
                              "SavitzkyGolay"),
                   msLevel. = uniqueMsLevels(x), ...) {
    if (!.check_ms_level(x, msLevel.))
        return(x)
    if (!is.integer(halfWindowSize) || length(halfWindowSize) != 1L ||
        halfWindowSize <= 0L)
        stop("Argument 'halfWindowSize' has to be an integer of length 1 ",
             "and > 0.")

    method <- match.arg(method)
    coef <- switch(method,
                   MovingAverage = coefMA(halfWindowSize),
                   WeightedMovingAverage = coefWMA(halfWindowSize),
                   SavitzkyGolay = coefSG(halfWindowSize, ...))

    x <- addProcessing(x, .peaks_smooth, halfWindowSize = halfWindowSize,
                       coef = coef, msLevel = msLevel., ...,
                       spectraVariables = "msLevel")
    x@processing <- .logging(x@processing, "Spectra smoothing with ", method,
                                           ", hws = ", halfWindowSize)
    x
})

#' @exportMethod addProcessing
#'
#' @importFrom ProtGenerics ProcessingStep
#'
#' @importMethodsFrom ProtGenerics addProcessing
#'
#' @importClassesFrom ProtGenerics ProcessingStep
#'
#' @importFrom methods .hasSlot
#'
#' @importFrom BiocGenerics updateObject
#'
#' @rdname Spectra
setMethod("addProcessing", "Spectra", function(object, FUN, ...,
                                               spectraVariables = character()) {
    if (missing(FUN))
        return(object)
    object@processingQueue <- c(object@processingQueue,
                                list(ProcessingStep(FUN, ARGS = list(...))))
    if (!.hasSlot(object, "processingQueueVariables"))
        object <- updateObject(object)
    object@processingQueueVariables <- union(object@processingQueueVariables,
                                             spectraVariables)
    validObject(object)
    object
})

#' @rdname Spectra
#'
#' @export
coreSpectraVariables <- function() .SPECTRA_DATA_COLUMNS

#' @rdname Spectra
setMethod("uniqueMsLevels", "Spectra", function(object, ...) {
    uniqueMsLevels(object@backend, ...)
})

#' @rdname Spectra
setMethod("backendBpparam", "Spectra", function(object, BPPARAM = bpparam()) {
    backendBpparam(object@backend, BPPARAM)
})

#' @rdname hidden_aliases
setMethod("combinePeaks", "list", function(object, ...) {
    .Deprecated("combinePeaksData", old = "combinePeaks",
                msg = paste0("'combinePeaks' for lists of peak matrices is ",
                             "deprecated; please use 'combinePeaksData' ",
                             "instead."))
    combinePeaksData(object, ...)
})

#' @rdname Spectra
#'
#' @exportMethod combinePeaks
setMethod("combinePeaks", "Spectra", function(object, tolerance = 0, ppm = 20,
                                              intensityFun = base::mean,
                                              mzFun = base::mean,
                                              weighted = TRUE,
                                              msLevel. = uniqueMsLevels(object),
                                              ...) {
    object <- addProcessing(
        object, .peaks_combine, ppm = ppm, tolerance = tolerance,
        intensityFun = intensityFun, mzFun = mzFun, weighted = weighted,
        msLevel = force(msLevel.), spectraVariables = "msLevel")
    object@processing <- .logging(
        object@processing, "Combining peaks within each spectrum with ppm = ",
        ppm, " and tolerance = ", tolerance, ".")
    object
})


#' @rdname Spectra
#'
#' @importFrom MsCoreUtils entropy nentropy
#'
#' @export
setMethod("entropy", "Spectra", function(object, normalized = TRUE) {
    if (length(object)) {
    if (normalized) entropy_fun <- nentropy
    else entropy_fun <- entropy
    unlist(.peaksapply(
        object, FUN = function(pks, ...) entropy_fun(pks[, "intensity"])),
        use.names = FALSE
        )
    } else numeric()
})
#' @rdname Spectra
setMethod("entropy", "ANY", function(object, ...) {
    MsCoreUtils::entropy(object)
})

#' @export
#' @rdname Spectra
#'
#' @param spectraVars `character()` indicating what spectra variables to add to
#'     the `DataFrame`. Default is `spectraVariables(object)`, i.e. all
#'     available variables.
#'
#' @examples
#'
#' ## Convert a subset of the Spectra object to a long DataFrame.
#' asDataFrame(sciex, i = 1:3, spectraVars = c("rtime", "msLevel"))
asDataFrame <- function(object, i = seq_along(object),
                        spectraVars = spectraVariables(object)) {
    stopifnot(inherits(object, "Spectra"))
    object <- object[i]
    n <- sapply(peaksData(object), nrow)
    v <- spectraData(object)[rep(seq_along(object), n), spectraVars]
    p <- do.call(rbind, as.list(peaksData(object)))
    cbind(p, v)
}

#' @title Estimate Precursor Intensities
#'
#' @description
#'
#' Some MS instrument manufacturers don't provide precursor intensities for
#' fragment spectra. These can however be estimated, given that also MS1
#' spectra are available. The `estimatePrecursorIntensity()` funtion defines the
#' precursor intensities for MS2 spectra using the intensity of the matching
#' MS1 peak from the closest MS1 spectrum (i.e. the last MS1 spectrum measured
#' before the respective MS2 spectrum). With `method = "interpolation"` it is
#' also possible to calculate the precursor intensity based on an interpolation
#' of intensity values (and retention times) of the matching MS1 peaks from the
#' previous and next MS1 spectrum. See below for an example.
#'
#' @param object `Spectra` with MS1 and MS2 spectra.
#'
#' @param ppm `numeric(1)` with the maximal allowed relative difference of m/z
#'     values between the precursor m/z of a spectrum and the m/z of the
#'     respective ion on the MS1 scan.
#'
#' @param tolerance `numeric(1)` with the maximal allowed difference of m/z
#'     values between the precursor m/z of a spectrum and the m/z of the
#'     respective ion on the MS1 scan.
#'
#' @param method `character(1)` defining whether the precursor intensity
#'     should be estimated on the previous MS1 spectrum (`method = "previous"`,
#'     the default) or based on an interpolation on the previous and next
#'     MS1 spectrum (`method = "interpolation"`).
#'
#' @param msLevel. `integer(1)` the MS level for which precursor intensities
#'     should be estimated. Defaults to `2L`.
#'
#' @param f `factor` (or vector to be coerced to `factor`) defining which
#'     spectra belong to the same original data file (sample).
#'     Defaults to `f = dataOrigin(x)`.
#'
#' @param BPPARAM Parallel setup configuration. See [bpparam()] for more
#'     information. This is passed directly to the [backendInitialize()] method
#'     of the [MsBackend-class].
#'
#' @author Johannes Rainer with feedback and suggestions from Corey Broeckling
#'
#' @importMethodsFrom ProtGenerics estimatePrecursorIntensity
#'
#' @exportMethod estimatePrecursorIntensity
#'
#' @rdname estimatePrecursorIntensity
#'
#' @examples
#'
#' #' ## Calculating the precursor intensity for MS2 spectra:
#' ##
#' ## Some MS instrument manufacturer don't report the precursor intensities
#' ## for MS2 spectra. The `estimatePrecursorIntensity` function can be used
#' ## in these cases to calculate the precursor intensity on MS1 data. Below
#' ## we load an mzML file from a vendor providing precursor intensities and
#' ## compare the estimated and reported precursor intensities.
#' tmt <- Spectra(msdata::proteomics(full.names = TRUE)[5],
#'     backend = MsBackendMzR())
#' pmi <- estimatePrecursorIntensity(tmt)
#' plot(pmi, precursorIntensity(tmt))
#'
#' ## We can also replace the original precursor intensity values with the
#' ## newly calculated ones
#' tmt$precursorIntensity <- pmi
setMethod(
    "estimatePrecursorIntensity", "Spectra",
    function(object, ppm = 20, tolerance = 0,
             method = c("previous", "interpolation"),
             msLevel. = 2L, f = dataOrigin(object), BPPARAM = bpparam()) {
        if (is.factor(f))
            f <- as.character(f)
        f <- factor(f, levels = unique(f))
        BPPARAM <- backendBpparam(object@backend, BPPARAM)
        unlist(bplapply(split(object, f),
                        FUN = .estimate_precursor_intensity, ppm = ppm,
                        tolerance = tolerance, method = method,
                        msLevel = msLevel., BPPARAM = BPPARAM),
               use.names = FALSE)
    })
