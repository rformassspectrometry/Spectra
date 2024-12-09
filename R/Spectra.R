#' @include hidden_aliases.R
NULL

################################################################################
##
## Spectra class, creation, data representation, export
##
################################################################################

#' @title The Spectra class to manage and access MS data
#'
#' @name Spectra
#'
#' @aliases Spectra-class
#' @aliases Spectra
#' @aliases setBackend
#' @aliases export
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
#' - [spectraData()] for accessing and using MS data through `Spectra` objects.
#' - [filterMsLevel()] to subset and filter `Spectra` objects.
#' - [plotSpectra()] for visualization of `Spectra` orbjects.
#' - [processingChunkSize()] for information on parallel and chunk-wise data
#'   processing.
#' - [combineSpectra()] for merging, aggregating and splitting of `Spectra`
#'   objects.
#' - [combinePeaks()] for merging and aggregating `Spectra`'s mass peaks data.
#' - [addProcessing()] for data analysis functions.
#' - [compareSpectra()] for spectra similarity calculations.
#'
#' @param backend For `Spectra()`: [MsBackend-class] to be used as backend. See
#'     section on creation of `Spectra` objects for details. For `setBackend()`:
#'     instance of [MsBackend-class] that supports `setBackend()` (i.e. for
#'     which `supportsSetBackend()` returns `TRUE`). Such backends have a
#'     parameter `data` in their `backendInitialize()` function that support
#'     passing the full spectra data to the initialize method. See section on
#'     creation of `Spectra` objects for details.
#'     For `export()`: [MsBackend-class] to be used to export the data.
#'
#' @param BPPARAM Parallel setup configuration. See [bpparam()] for more
#'     information. This is passed directly to the [backendInitialize()] method
#'     of the [MsBackend-class].
#'
#' @param f For `setBackend()`: factor defining how to split the data
#'     for parallelized copying of the spectra data to the new backend. For
#'     some backends changing this parameter can lead to errors. Defaults to
#'     [processingChunkFactor()].
#'
#' @param metadata For `Spectra()`: optional `list` with metadata information.
#'
#' @param object For `Spectra()`: an object to instantiate the `Spectra`
#'     object and initialize the with data.. See section on creation of
#'     `Spectra` objects for details. For all other methods a `Spectra` object.
#'
#' @param processingQueue For `Spectra()`: optional `list` of
#'     [ProcessingStep-class] objects.
#'
#' @param source For `Spectra()`: instance of [MsBackend-class] that can be
#'     used to import spectrum data from the provided files. See section
#'    *Creation of objects* for more details.
#'
#' @param value For `dataStorageBasePath()`: A `character` vector that defines
#'     the base directory where the data storage files can be found.
#'
#' @param ... Additional arguments.
#'
#' @section Data stored in a `Spectra` object:
#'
#' The `Spectra` object is a container for MS data that includes mass peak
#' data (*m/z* and related intensity values, also referred to as *peaks data*
#' in the context of `Spectra`) and metadata of individual spectra (so called
#' *spectra variables*). While a core set of spectra variables (the
#' `coreSpectraVariables()`) are guaranteed to be provided by a
#' `Spectra`, it is possible to add arbitrary additional spectra variables to
#' a `Spectra` object.
#'
#' The `Spectra` object is designed to contain MS data of a (large) set of mass
#' spectra. The data is organized *linearly* and can be thought of a list of
#' mass spectra, i.e. each element in the `Spectra` is one spectrum.
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
#' spectrum data once it is retrieved. This enables data manipulation
#' operations also for *read only* data representations. For some backends that
#' allow to write data back to the data storage (such as the
#' [MsBackendMemory()], [MsBackendDataFrame()] and [MsBackendHdf5Peaks()]) it
#' is possible to apply to queue with the [applyProcessing()] function (see
#' the [applyProcessing()] function for details).
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
#' ##  --------  CREATION OF SPECTRA OBJECTS  --------
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
#' ##  --------  CHANGING DATA REPRESENTATIONS  --------
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
#' ##  -------- DATA EXPORT  --------
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
NULL

#' The Spectra class
#'
#' The [Spectra] class encapsulates data and meta-data for mass
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

#' @rdname Spectra
setMethod("Spectra", "missing", function(object, processingQueue = list(),
                                         metadata = list(), ...,
                                         backend = MsBackendMemory(),
                                         BPPARAM = bpparam()) {
    if (length(backend))
        new("Spectra", metadata = metadata, processingQueue = processingQueue,
            backend = backend)
    else callNextMethod()
})

#' @rdname Spectra
setMethod("Spectra", "MsBackend", function(object, processingQueue = list(),
                                           metadata = list(), ...,
                                           BPPARAM = bpparam()) {
    new("Spectra", metadata = metadata, processingQueue = processingQueue,
        backend = object)
})

#' @rdname Spectra
#'
#' @importFrom methods callNextMethod
setMethod("Spectra", "character", function(object, processingQueue = list(),
                                           metadata = list(),
                                           source = MsBackendMzR(),
                                           backend = source,
                                           ..., BPPARAM = bpparam()) {
    sp <- .create_spectra(object, processingQueue = processingQueue,
                          metadata = metadata, backend = source,
                          ..., BPPARAM = BPPARAM)
    if (class(source)[1L] != class(backend)[1L])
        setBackend(sp, backend, ..., BPPARAM = backendBpparam(backend, BPPARAM))
    else sp
})

#' @rdname Spectra
setMethod("Spectra", "ANY", function(object, processingQueue = list(),
                                     metadata = list(),
                                     source = MsBackendMemory(),
                                     backend = source,
                                     ..., BPPARAM = bpparam()) {
    sp <- .create_spectra(object, processingQueue = processingQueue,
                          metadata = metadata, backend = source,
                          ..., BPPARAM = BPPARAM)
    if (class(source)[1L] != class(backend)[1L])
        setBackend(sp, backend, ..., BPPARAM = backendBpparam(backend, BPPARAM))
    else sp
})

.create_spectra <- function(object, processingQueue = list(), metadata = list(),
                            backend = MsBackendMemory(), ...,
                            BPPARAM = bpparam()) {
    if (missing(object))
        backend <- backendInitialize(
            backend, ..., BPPARAM = backendBpparam(backend, BPPARAM))
    else backend <- backendInitialize(
             backend, object, ..., BPPARAM = backendBpparam(backend, BPPARAM))
    new("Spectra", metadata = metadata, processingQueue = processingQueue,
        backend = backend)
}

#' @rdname Spectra
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
                    bknds <- extractByIndex(
                        bknds, order(unlist(split(seq_along(bknds), f),
                                            use.names = FALSE)))
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

#' @rdname Spectra
#'
#' @export
setMethod("export", "Spectra",
          function(object, backend, ...) {
              if (missing(backend))
                  stop("Parameter 'backend' is required.")
              export(backend, object, ...)
          })

#' @rdname Spectra
setMethod("dataStorageBasePath", "Spectra", function(object) {
    dataStorageBasePath(object@backend)
})

#' @rdname Spectra
setReplaceMethod("dataStorageBasePath", "Spectra", function(object, value) {
    dataStorageBasePath(object@backend) <- value
    object
})

################################################################################
##
## Accessing and adding/setting/changing MS data.
##
################################################################################

#' @title Accessing mass spectrometry data
#'
#' @name spectraData
#'
#' @aliases acquisitionNum
#' @aliases centroided
#' @aliases collisionEnergy
#' @aliases dataOrigin
#' @aliases dataStorage
#' @aliases intensity
#' @aliases ionCount
#' @aliases isCentroided
#' @aliases isEmpty
#' @aliases isolationWindowLowerMz
#' @aliases isolationWindowUpperMz
#' @aliases isolationWindowTargetMz
#' @aliases lengths
#' @aliases msLevel
#' @aliases mz
#' @aliases peaksData
#' @aliases peaksVariables
#' @aliases polarity
#' @aliases precursorCharge
#' @aliases precursorIntensity
#' @aliases precursorMz
#' @aliases rtime
#' @aliases scanIndex
#' @aliases smoothed
#' @aliases spectraData
#' @aliases spectraNames
#' @aliases spectraVariables
#' @aliases tic
#' @aliases uniqueMsLevels
#'
#' @description
#'
#' As detailed in the documentation of the [Spectra] class, a `Spectra` object
#' is a container for mass spectrometry (MS) data that includes both the mass
#' peaks data (or *peaks data*, generally *m/z* and intensity values) as well
#' as spectra metadata (so called *spectra variables*). Spectra variables
#' generally define one value per spectrum, while for peaks variables one value
#' per mass peak is defined and hence multiple values per spectrum (depending
#' on the number of mass peaks of a spectrum).
#'
#' Data can be extracted from a `Spectra` object using dedicated accessor
#' functions or also using the `$` operator. Depending on the backend class
#' used by the `Spectra` to represent the data, data can also be added or
#' replaced (again, using dedicated functions or using `$<-`).
#'
#'
#' @param BPPARAM Parallel setup configuration. See [bpparam()] for more
#'     information. See also [processingChunkSize()] for more information
#'     on parallel processing.
#'
#' @param columns For `spectraData()` accessor: optional `character` with
#'     column names (spectra variables) that should be included in the
#'     returned `DataFrame`. By default, all columns are returned.
#'     For `peaksData()` accessor: optional `character` with requested columns
#'     in the individual `matrix` of the returned `list`. Defaults to
#'     `c("mz", "value")` but any values returned by `peaksVariables(object)`
#'     with `object` being the `Spectra` object are supported.
#'
#' @param f For `intensity()`, `mz()` and `peaksData()`: factor defining how
#'     data should be chunk-wise loaded an processed. Defaults to
#'     [processingChunkFactor()].
#'
#' @param i For `asDataFrame()`: A `numeric` indicating which scans to coerce
#'     to a `DataFrame` (default is `seq_along(object)`).
#'
#' @param initial For `tic()`: `logical(1)` whether the initially
#'     reported total ion current should be reported, or whether the
#'     total ion current should be (re)calculated on the actual data
#'     (`initial = FALSE`, same as `ionCount()`).
#'
#' @param j For `[`: not supported.
#'
#' @param name For `$` and `$<-`: the name of the spectra variable to return
#'     or set.
#'
#' @param object A `Spectra` object.
#'
#' @param spectraVars `character()` indicating what spectra variables to add to
#'     the `DataFrame`. Default is `spectraVariables(object)`, i.e. all
#'     available variables.
#'
#' @param use.names For `lengths()`: ignored.
#'
#' @param value A vector with values to replace the respective spectra
#'     variable. Needs to be of the correct data type for the spectra variable.
#'
#' @param x A `Spectra` object.
#'
#' @param ... Additional arguments.
#'
#'
#' @section Spectra variables:
#'
#' A common set of *core spectra variables* are defined for `Spectra`. These
#' have a pre-defined data type and each `Spectra` will return a value for
#' these if requested. If no value for a spectra variable is defined, a missing
#' value (of the correct data type) is returned. The list of core spectra
#' variables and their respective data type is:
#'
#' - *acquisitionNum* `integer(1)`: the index of acquisition of a spectrum
#'   during an MS run.
#' - *centroided* `logical(1)`: whether the spectrum is in profile or centroid
#'   mode.
#' - *collisionEnergy* `numeric(1)`: collision energy used to create an MSn
#'   spectrum.
#' - *dataOrigin* `character(1)`: the *origin* of the spectrum's data, e.g. the
#'   mzML file from which it was read.
#' - *dataStorage* `character(1)`: the (current) storage location of the
#'   spectrum data. This value depends on the backend used to handle and
#'   provide the data. For an *in-memory* backend like the `MsBackendDataFrame`
#'   this will be `"<memory>"`, for an on-disk backend such as the
#'   `MsBackendHdf5Peaks` it will be the name of the HDF5 file where the
#'   spectrum's peak data is stored.
#' - *isolationWindowLowerMz* `numeric(1)`: lower m/z for the isolation
#'   window in which the (MSn) spectrum was measured.
#' - *isolationWindowTargetMz* `numeric(1)`: the target m/z for the isolation
#'   window in which the (MSn) spectrum was measured.
#' - *isolationWindowUpperMz* `numeric(1)`: upper m/z for the isolation window
#'   in which the (MSn) spectrum was measured.
#' - *msLevel* `integer(1)`: the MS level of the spectrum.
#' - *polarity* `integer(1)`: the polarity of the spectrum (`0` and `1`
#'   representing negative and positive polarity, respectively).
#' - *precScanNum* `integer(1)`: the scan (acquisition) number of the precursor
#'   for an MSn spectrum.
#' - *precursorCharge* `integer(1)`: the charge of the precursor of an MSn
#'   spectrum.
#' - *precursorIntensity* `numeric(1)`: the intensity of the precursor of an
#'   MSn spectrum.
#' - *precursorMz* `numeric(1)`: the m/z of the precursor of an MSn spectrum.
#' - *rtime* `numeric(1)`: the retention time of a spectrum.
#' - *scanIndex* `integer(1)`: the index of a spectrum within a (raw) file.
#' - *smoothed* `logical(1)`: whether the spectrum was smoothed.
#'
#' For each of these spectra variable a dedicated accessor function is defined
#' (such as `msLevel()` or `rtime()`) that allows to extract the values of
#' that spectra variable for all spectra in a `Spectra` object. Also,
#' replacement functions are defined, but not all backends might support
#' replacing values for spectra variables. As described above, additional
#' spectra variables can be defined or added. The `spectraVariables()` function
#' can be used to
#'
#' Values for multiple spectra variables, or all spectra vartiables* can be
#' extracted with the `spectraData()` function.
#'
#'
#' @section Peaks variables:
#'
#' `Spectra` also provide mass peak data with the *m/z* and intensity values
#' being the *core* peaks variables:
#'
#' - *intensity* `numeric`: intensity values for the spectrum's peaks.
#' - *mz* `numeric`: the m/z values for the spectrum's peaks.
#'
#' Values for these can be extracted with the `mz()` and `intensity()`
#' functions, or the `peaksData()` function. The former functions return a
#' `NumericList` with the respective values, while the latter returns a `List`
#' with `numeric` two-column matrices. The list of peaks matrices can also
#' be extracted using `as(x, "list")` or `as(x, "SimpleList")` with `x` being
#' a `Spectra` object.
#'
#' Some `Spectra`/backends provide also values for additional peaks variables.
#' The set of available peaks variables can be extracted with the
#' `peaksVariables()` function.
#'
#'
#' @section Functions to access MS data:
#'
#' The set of available functions to extract data from, or set data in, a
#' `Spectra` object are (in alphabetical order) listed below. Note that there
#' are also other functions to extract information from a `Spectra` object
#' documented in [addProcessing()].
#'
#' - `$`, `$<-`: gets (or sets) a spectra variable for all spectra in `object`.
#'   See examples for details. Note that replacing values of a peaks variable
#'   is not supported with a non-empty processing queue, i.e. if any filtering
#'   or data manipulations on the peaks data was performed. In these cases
#'   [applyProcessing()] needs to be called first to apply all cached data
#'   operations.
#'
#' - `[[`, `[[<-`: access or set/add a single spectrum variable (column) in the
#'   backend.
#'
#' - `acquisitionNum()`: returns the acquisition number of each
#'   spectrum. Returns an `integer` of length equal to the number of
#'   spectra (with `NA_integer_` if not available).
#'
#' - `asDataFrame()`: converts the `Spectra` to a `DataFrame` (in long format)
#'   contining all data. Returns a `DataFrame`.
#'
#' - `centroided()`, `centroided<-`: gets or sets the centroiding
#'   information of the spectra. `centroided()` returns a `logical`
#'   vector of length equal to the number of spectra with `TRUE` if a
#'   spectrum is centroided, `FALSE` if it is in profile mode and `NA`
#'   if it is undefined. See also `isCentroided()` for estimating from
#'   the spectrum data whether the spectrum is centroided. `value`
#'   for `centroided<-` is either a single `logical` or a `logical` of
#'   length equal to the number of spectra in `object`.
#'
#' - `collisionEnergy()`, `collisionEnergy<-`: gets or sets the
#'   collision energy for all spectra in `object`. `collisionEnergy()`
#'   returns a `numeric` with length equal to the number of spectra
#'   (`NA_real_` if not present/defined), `collisionEnergy<-` takes a
#'   `numeric` of length equal to the number of spectra in `object`.
#'
#' - `coreSpectraVariables()`: returns the *core* spectra variables along with
#'   their expected data type.
#'
#' - `dataOrigin()`, `dataOrigin<-`: gets or sets the *data origin* for each
#'   spectrum. `dataOrigin()` returns a `character` vector (same length than
#'   `object`) with the origin of the spectra. `dataOrigin<-` expects a
#'   `character` vector (same length than `object`) with the replacement
#'   values for the data origin of each spectrum.
#'
#' - `dataStorage()`: returns a `character` vector (same length than `object`)
#'   with the data storage location of each spectrum.
#'
#' - `intensity()`: gets the intensity values from the spectra. Returns
#'   a [NumericList()] of `numeric` vectors (intensity values for each
#'   spectrum). The length of the list is equal to the number of
#'   `spectra` in `object`.
#'
#' - `ionCount()`: returns a `numeric` with the sum of intensities for
#'   each spectrum. If the spectrum is empty (see `isEmpty()`),
#'   `NA_real_` is returned.
#'
#' - `isCentroided()`: a heuristic approach assessing if the spectra in
#'   `object` are in profile or centroided mode. The function takes
#'   the `qtl`th quantile top peaks, then calculates the difference
#'   between adjacent m/z value and returns `TRUE` if the first
#'   quartile is greater than `k`. (See `Spectra:::.isCentroided()` for
#'   the code.)
#'
#' - `isEmpty()`: checks whether a spectrum in `object` is empty
#'   (i.e. does not contain any peaks). Returns a `logical` vector of
#'   length equal number of spectra.
#'
#' - `isolationWindowLowerMz()`, `isolationWindowLowerMz<-`: gets or sets the
#'   lower m/z boundary of the isolation window.
#'
#' - `isolationWindowTargetMz()`, `isolationWindowTargetMz<-`: gets or sets the
#'   target m/z of the isolation window.
#'
#' - `isolationWindowUpperMz()`, `isolationWindowUpperMz<-`: gets or sets the
#'   upper m/z boundary of the isolation window.
#'
#' - `length()`: gets the number of spectra in the object.
#'
#' - `lengths()`: gets the number of peaks (m/z-intensity values) per
#'   spectrum. Returns an `integer` vector (length equal to the
#'   number of spectra). For empty spectra, `0` is returned.
#'
#' - `msLevel()`: gets the spectra's MS level. Returns an integer vector (names
#'   being spectrum names, length equal to the number of spectra) with the MS
#'   level for each spectrum.
#'
#' - `mz()`: gets the mass-to-charge ratios (m/z) from the
#'   spectra. Returns a [NumericList()] or length equal to the number of
#'   spectra, each element a `numeric` vector with the m/z values of
#'   one spectrum.
#'
#' - `peaksData()`: gets the *peaks* data for all spectra in `object`. Peaks
#'   data consist of the m/z and intensity values as well as possible additional
#'   annotations (variables) of all peaks of each spectrum. The function
#'   returns a [SimpleList()] of two dimensional arrays (either `matrix` or
#'   `data.frame`), with each array providing the values for the requested
#'   *peak variables* (by default `"mz"` and `"intensity"`). Optional parameter
#'   `columns` is passed to the backend's `peaksData()` function to allow
#'   the selection of specific (or additional) peaks variables (columns) that
#'   should be extracted (if available). Importantly,
#'   it is **not** guaranteed that each backend supports this parameter (while
#'   each backend must support extraction of `"mz"` and `"intensity"` columns).
#'   Parameter `columns` defaults to `c("mz", "intensity")` but any value
#'   returned by `peaksVariables(object)` is supported.
#'   Note also that it is possible to extract the peak data with
#'   `as(x, "list")` and `as(x, "SimpleList")` as a `list` and `SimpleList`,
#'   respectively. Note however that, in contrast to `peaksData()`, `as()`
#'   does not support the parameter `columns`.
#'
#' - `peaksVariables()`: lists the available variables for mass peaks provided
#'   by the backend. Default peak variables are `"mz"` and `"intensity"` (which
#'   all backends need to support and provide), but some backends might provide
#'   additional variables.
#'   These variables correspond to the column names of the peak data array
#'   returned by `peaksData()`.
#'
#' - `polarity()`, `polarity<-`: gets or sets the polarity for each
#'   spectrum. `polarity()` returns an `integer` vector (length equal
#'   to the number of spectra), with `0` and `1` representing negative
#'   and positive polarities, respectively. `polarity<-` expects an
#'   `integer` vector of length 1 or equal to the number of spectra.
#'
#' - `precursorCharge()`, `precursorIntensity()`, `precursorMz()`,
#'   `precScanNum()`, `precAcquisitionNum()`: gets the charge (`integer`),
#'   intensity (`numeric`), m/z (`numeric`), scan index (`integer`)
#'   and acquisition number (`interger`) of the precursor for MS level >
#'   2 spectra from the object. Returns a vector of length equal to
#'   the number of spectra in `object`. `NA` are reported for MS1
#'   spectra of if no precursor information is available.
#'
#' - `rtime()`, `rtime<-`: gets or sets the retention times (in seconds)
#'   for each spectrum.  `rtime()` returns a `numeric` vector (length
#'   equal to the number of spectra) with the retention time for each
#'   spectrum.  `rtime<-` expects a numeric vector with length equal
#'   to the number of spectra.
#'
#' - `scanIndex()`: returns an `integer` vector with the *scan index*
#'   for each spectrum. This represents the relative index of the
#'   spectrum within each file. Note that this can be different to the
#'   `acquisitionNum` of the spectrum which represents the index of the
#'   spectrum during acquisition/measurement (as reported in the mzML file).
#'
#' - `smoothed()`,`smoothed<-`: gets or sets whether a spectrum is
#'   *smoothed*. `smoothed()` returns a `logical` vector of length equal
#'   to the number of spectra. `smoothed<-` takes a `logical` vector
#'   of length 1 or equal to the number of spectra in `object`.
#'
#' - `spectraData()`: gets general spectrum metadata (annotation, also called
#'   header). `spectraData()` returns a `DataFrame`. Note that this
#'   method does by default **not** return m/z or intensity values.
#'
#' - `spectraData<-`: **replaces** the full spectra data of the `Spectra`
#'   object with the one provided with `value`. The `spectraData<-` function
#'   expects a `DataFrame` to be passed as value with the same number of rows
#'   as there a spectra in `object`. Note that replacing values of
#'   peaks variables is not supported with a non-empty processing queue, i.e.
#'   if any filtering or data manipulations on the peaks data was performed.
#'   In these cases [applyProcessing()] needs to be called first to apply all
#'   cached data operations and empty the processing queue.
#'
#' - `spectraNames()`, `spectraNames<-`: gets or sets the spectra names.
#'
#' - `spectraVariables()`: returns a `character` vector with the
#'   available spectra variables (columns, fields or attributes of each
#'   spectrum) available in `object`. Note that `spectraVariables()` does not
#'   list the *peak variables* (`"mz"`, `"intensity"` and eventual additional
#'   annotations for each MS peak). Peak variables are returned by
#'   `peaksVariables()`.
#'
#' - `tic()`: gets the total ion current/count (sum of signal of a
#'   spectrum) for all spectra in `object`. By default, the value
#'   reported in the original raw data file is returned. For an empty
#'   spectrum, `0` is returned.
#'
#' - `uniqueMsLevels()`: get the unique MS levels available in `object`. This
#'   function is supposed to be more efficient than `unique(msLevel(object))`.
#'
#' @md
#'
#' @seealso
#'
#' - [addProcessing()] for functions to analyze `Spectra`.
#'
#' - [Spectra] for a general description of the `Spectra` object.
#'
#' @author Sebastian Gibb, Johannes Rainer, Laurent Gatto, Philippine Louail
#'
#' @examples
#'
#' ## Create a Spectra from mzML files and use the `MsBackendMzR` on-disk
#' ## backend.
#' sciex_file <- dir(system.file("sciex", package = "msdata"),
#'     full.names = TRUE)
#' sciex <- Spectra(sciex_file, backend = MsBackendMzR())
#' sciex
#'
#' ## Get the number of spectra in the data set
#' length(sciex)
#'
#' ## Get the number of mass peaks per spectrum - limit to the first 6
#' lengths(sciex) |> head()
#'
#' ## Get the MS level for each spectrum - limit to the first 6 spectra
#' msLevel(sciex) |> head()
#'
#' ## Alternatively, we could also use $ to access a specific spectra variable.
#' ## This could also be used to add additional spectra variables to the
#' ## object (see further below).
#' sciex$msLevel |> head()
#'
#' ## Get the intensity and m/z values.
#' intensity(sciex)
#' mz(sciex)
#'
#' ## Convert a subset of the Spectra object to a long DataFrame.
#' asDataFrame(sciex, i = 1:3, spectraVars = c("rtime", "msLevel"))
#'
#' ## Create a Spectra providing a `DataFrame` containing the spectrum data.
#'
#' spd <- DataFrame(msLevel = c(1L, 2L), rtime = c(1.1, 1.2))
#' spd$mz <- list(c(100, 103.2, 104.3, 106.5), c(45.6, 120.4, 190.2))
#' spd$intensity <- list(c(200, 400, 34.2, 17), c(12.3, 15.2, 6.8))
#'
#' s <- Spectra(spd)
#' s
#'
#' ## List all available spectra variables (i.e. spectrum data and metadata).
#' spectraVariables(s)
#'
#' ## For all *core* spectrum variables accessor functions are available. These
#' ## return NA if the variable was not set.
#' centroided(s)
#' dataStorage(s)
#' rtime(s)
#' precursorMz(s)
#'
#' ## The core spectra variables are:
#' coreSpectraVariables()
#'
#' ## Add an additional metadata column.
#' s$spectrum_id <- c("sp_1", "sp_2")
#'
#' ## List spectra variables, "spectrum_id" is now also listed
#' spectraVariables(s)
#'
#' ## Get the values for the new spectra variable
#' s$spectrum_id
#'
#' ## Extract specific spectra variables.
#' spectraData(s, columns = c("spectrum_id", "msLevel"))
#'
#'
#' ##  --------  PEAKS VARIABLES AND DATA  --------
#'
#' ## Get the peak data (m/z and intensity values).
#' pks <- peaksData(s)
#' pks
#' pks[[1]]
#' pks[[2]]
#'
#' ## Note that we could get the same resulb by coercing the `Spectra` to
#' ## a `list` or `SimpleList`:
#' as(s, "list")
#' as(s, "SimpleList")
#'
#' ## Or use `mz()` and `intensity()` to extract the m/z and intensity values
#' ## separately
#' mz(s)
#' intensity(s)
#'
#' ## Some `MsBackend` classes provide support for arbitrary peaks variables
#' ## (in addition to the mandatory `"mz"` and `"intensity"` values. Below
#' ## we create a simple data frame with an additional peak variable `"pk_ann"`
#' ## and create a `Spectra` with a `MsBackendMemory` for that data.
#' ## Importantly the number of values (per spectrum) need to be the same
#' ## for all peak variables.
#'
#' tmp <- data.frame(msLevel = c(2L, 2L), rtime = c(123.2, 123.5))
#' tmp$mz <- list(c(103.1, 110.4, 303.1), c(343.2, 453.1))
#' tmp$intensity <- list(c(130.1, 543.1, 40), c(0.9, 0.45))
#' tmp$pk_ann <- list(c(NA_character_, "A", "P"), c("B", "P"))
#'
#' ## Create the Spectra. With parameter `peaksVariables` we can define
#' ## the columns in `tmp` that contain peaks variables.
#' sps <- Spectra(tmp, source = MsBackendMemory(),
#'     peaksVariables = c("mz", "intensity", "pk_ann"))
#' peaksVariables(sps)
#'
#' ## Extract just the m/z and intensity values
#' peaksData(sps)[[1L]]
#'
#' ## Extract the full peaks data
#' peaksData(sps, columns = peaksVariables(sps))[[1L]]
#'
#' ## Access just the pk_ann variable
#' sps$pk_ann
#'
#'
NULL

#' @importFrom methods setAs
setAs("Spectra", "list", function(from, to) {
    .peaksapply(from)
})

setAs("Spectra", "SimpleList", function(from, to) {
    peaksData(from)
})

#' @export
#'
#' @rdname spectraData
asDataFrame <- function(object, i = seq_along(object),
                        spectraVars = spectraVariables(object)) {
    stopifnot(inherits(object, "Spectra"))
    object <- object[i]
    n <- sapply(peaksData(object), nrow)
    v <- spectraData(object)[rep(seq_along(object), n), spectraVars]
    p <- do.call(rbind, as.list(peaksData(object)))
    cbind(p, v)
}

#' @rdname spectraData
#'
#' @export
setMethod("acquisitionNum", "Spectra", function(object)
    acquisitionNum(object@backend))

#' @rdname spectraData
setMethod("centroided", "Spectra", function(object) {
    centroided(object@backend)
})

#' @rdname spectraData
setReplaceMethod("centroided", "Spectra", function(object, value) {
    centroided(object@backend) <- value
    object
})

#' @rdname spectraData
setMethod("collisionEnergy", "Spectra", function(object) {
    collisionEnergy(object@backend)
})

#' @rdname spectraData
setReplaceMethod("collisionEnergy", "Spectra", function(object, value) {
    collisionEnergy(object@backend) <- value
    object
})

#' @rdname spectraData
#'
#' @export
coreSpectraVariables <- function() .SPECTRA_DATA_COLUMNS

#' @rdname spectraData
setMethod("dataOrigin", "Spectra", function(object) dataOrigin(object@backend))

#' @rdname spectraData
setReplaceMethod("dataOrigin", "Spectra", function(object, value) {
    dataOrigin(object@backend) <- value
    object
})

#' @rdname spectraData
setMethod("dataStorage", "Spectra",
          function(object) dataStorage(object@backend))

#' @rdname spectraData
setMethod("intensity", "Spectra", function(object,
                                           f = processingChunkFactor(object),
                                           ...) {
    if (length(object@processingQueue) || length(f))
        NumericList(.peaksapply(object, FUN = function(z, ...) z[, 2],
                                f = f, ...), compress = FALSE)
    else intensity(object@backend)
})

#' @rdname spectraData
setMethod("ionCount", "Spectra", function(object) {
    if (length(object))
        unlist(.peaksapply(
            object, FUN = function(pks, ...) sum(pks[, 2], na.rm = TRUE)),
            use.names = FALSE)
    else numeric()
})

#' @rdname spectraData
setMethod("isCentroided", "Spectra", function(object, ...) {
    if (length(object))
        unlist(.peaksapply(object, FUN = .peaks_is_centroided),
               use.names = FALSE)
    else logical()
})

#' @rdname spectraData
setMethod("isEmpty", "Spectra", function(x) {
    if (length(x))
        unlist(.peaksapply(x, FUN = function(pks, ...) nrow(pks) == 0),
               use.names = FALSE)
    else logical()
})

#' @rdname spectraData
setMethod("isolationWindowLowerMz", "Spectra", function(object) {
    isolationWindowLowerMz(object@backend)
})

#' @rdname spectraData
setReplaceMethod("isolationWindowLowerMz", "Spectra", function(object, value) {
    isolationWindowLowerMz(object@backend) <- value
    object
})

#' @rdname spectraData
setMethod("isolationWindowTargetMz", "Spectra", function(object) {
    isolationWindowTargetMz(object@backend)
})

#' @rdname spectraData
setReplaceMethod("isolationWindowTargetMz", "Spectra", function(object, value) {
    isolationWindowTargetMz(object@backend) <- value
    object
})

#' @rdname spectraData
setMethod("isolationWindowUpperMz", "Spectra", function(object) {
    isolationWindowUpperMz(object@backend)
})

#' @rdname spectraData
setReplaceMethod("isolationWindowUpperMz", "Spectra", function(object, value) {
    isolationWindowUpperMz(object@backend) <- value
    object
})

#' @rdname spectraData
#'
#' @exportMethod length
setMethod("length", "Spectra", function(x) length(x@backend))

#' @rdname spectraData
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

#' @rdname spectraData
setMethod("msLevel", "Spectra", function(object) msLevel(object@backend))

#' @rdname spectraData
setMethod("mz", "Spectra", function(object, f = processingChunkFactor(object),
                                    ...) {
    if (length(object@processingQueue) || length(f))
        NumericList(.peaksapply(object, FUN = function(z, ...) z[, 1],
                                f = f, ...), compress = FALSE)
    else mz(object@backend)
})

#' @rdname spectraData
#'
#' @export
setMethod(
    "peaksData", "Spectra",
    function(object, columns = c("mz", "intensity"),
             f = processingChunkFactor(object), ..., BPPARAM = bpparam()) {
        if (length(object@processingQueue) || length(f))
            SimpleList(.peaksapply(object, columns = columns, f = f))
        else SimpleList(peaksData(object@backend, columns = columns))
    })

#' @rdname spectraData
setMethod("peaksVariables", "Spectra", function(object)
    peaksVariables(object@backend))

#' @rdname spectraData
setMethod("polarity", "Spectra", function(object) {
    polarity(object@backend)
})

#' @rdname spectraData
setReplaceMethod("polarity", "Spectra", function(object, value) {
    polarity(object@backend) <- value
    object
})

#' @rdname spectraData
setMethod("precScanNum", "Spectra", function(object) {
    precScanNum(object@backend)
})

#' @rdname spectraData
setMethod("precursorCharge", "Spectra", function(object) {
    precursorCharge(object@backend)
})

#' @rdname spectraData
setMethod("precursorIntensity", "Spectra", function(object) {
    precursorIntensity(object@backend)
})

#' @rdname spectraData
setMethod("precursorMz", "Spectra", function(object) {
    precursorMz(object@backend)
})

#' @rdname spectraData
setReplaceMethod("precursorMz", "Spectra", function(object, ..., value) {
    precursorMz(object@backend) <- value
    object
})

#' @rdname spectraData
setMethod("rtime", "Spectra", function(object) {
    rtime(object@backend)
})

#' @rdname spectraData
setReplaceMethod("rtime", "Spectra", function(object, value) {
    rtime(object@backend) <- value
    object
})

#' @rdname spectraData
setMethod("scanIndex", "Spectra", function(object) {
    scanIndex(object@backend)
})

#' @rdname spectraData
setMethod("smoothed", "Spectra", function(object) {
    smoothed(object@backend)
})

#' @rdname spectraData
setReplaceMethod("smoothed", "Spectra", function(object, value) {
    smoothed(object@backend) <- value
    object
})

#' @rdname spectraData
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

#' @rdname spectraData
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

#' @rdname spectraData
setMethod("spectraNames", "Spectra", function(object) {
    spectraNames(object@backend)
})

#' @rdname spectraData
setReplaceMethod("spectraNames", "Spectra", function(object, value) {
    spectraNames(object@backend) <- value
    object
})

#' @rdname spectraData
setMethod("spectraVariables", "Spectra", function(object) {
    setdiff(spectraVariables(object@backend), peaksVariables(object@backend))
})

#' @rdname spectraData
setMethod("tic", "Spectra", function(object, initial = TRUE) {
    if (!length(object))
        return(numeric())
    if (initial)
        tic(object@backend, initial = initial)
    else ionCount(object)
})

#' @rdname spectraData
setMethod("uniqueMsLevels", "Spectra", function(object, ...) {
    uniqueMsLevels(object@backend, ...)
})

#' @rdname spectraData
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

#' @rdname spectraData
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

#' @rdname spectraData
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

#' @rdname spectraData
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


################################################################################
##
## Merging, splitting and aggregating Spectra: length of Spectra is changed
##
################################################################################

#' @title Merging, aggregating and splitting Spectra
#'
#' @name combineSpectra
#'
#' @aliases combineSpectra
#' @aliases split
#' @aliases joinSpectraData
#' @aliases cbind2
#'
#' @description
#'
#' Various functions are availabe to combine, aggregate or split data from one
#' of more `Spectra` objects. These are:
#'
#' - `c()` and `concatenateSpectra()`: combines several `Spectra` objects into
#'   a single object. The resulting `Spectra` contains all data from all
#'   individual `Spectra`, i.e. the union of all their spectra variables.
#'   Concatenation will fail if the processing queue of any of the `Spectra`
#'   objects is not empty or if different backends are used for the `Spectra`
#'   objects. In such cases it is suggested to first change the backends of
#'   all `Spectra` to the same type of backend (using the [setBackend()]
#'   function and to eventually (if needed) apply the processing queue using
#'   the [applyProcessing()] function.
#'
#' - `cbind2()`: Appends multiple spectra variables from a `data.frame`,
#'   `DataFrame` or `matrix` to the `Spectra` object at once. The order of
#'   the values (rows) in `y` has to match the order of spectra in `x`. The
#'   function does not allow to replace existing spectra variables. `cbind2()`
#'   returns a `Spectra` object with the appended spectra variables. For a more
#'   controlled way of adding spectra variables, see the `joinSpectraData()`
#'   function.
#'
#' - `combineSpectra()`: combines sets of spectra (defined with parameter `f`)
#'   into a single spectrum per set aggregating their MS data (i.e. their
#'   *peaks data* matrices with the *m/z* and intensity values of their
#'   mass peaks). The spectra variable values of the first spectrum per set
#'   are reported for the combined spectrum. The peak matrices of the spectra
#'   per set are combined using the function specified with parameter `FUN`
#'   which uses by default the [combinePeaksData()] function. See the
#'   documentation of [combinePeaksData()] for details on the aggregation of
#'   the peak data and the package vignette for examples.
#'   The sets of spectra can be specified with parameter `f` which is expected
#'   to be a `factor` or `vector` of length equal to the length of the
#'   `Spectra` specifying to which set a spectrum belongs to. The function
#'   returns a `Spectra` of length equal to the unique levels of `f`. The
#'   optional parameter `p` allows to define how the `Spectra` should be
#'   split for potential parallel processing. The default is
#'   `p = x$dataStorage` and hence a per storage file parallel processing is
#'   applied for `Spectra` with on disk data representations (such as the
#'   [MsBackendMzR()]). This also prevents that spectra from different data
#'   files/samples are combined (eventually use e.g. `p = x$dataOrigin` or any
#'   other spectra variables defining the originating samples for a spectrum).
#'   Before combining the peaks data, all eventual present processing steps are
#'   applied (by calling [applyProcessing()] on the `Spectra`). This function
#'   will replace the original *m/z* and intensity values of a `Spectra` hence
#'   it can not be called on a `Spectra` with a *read-only* backend. In such
#'   cases, the backend should be changed to a *writeable* backend before
#'   using the [setBackend()] function (to e.g. a [MsBackendMemory()] backend).
#'
#' - `joinSpectraData()`: Individual spectra variables can be directly
#'    added with the `$<-` or `[[<-` syntax. The `joinSpectraData()`
#'    function allows to merge a `DataFrame` to the existing spectra
#'    data of a `Spectra`. This function diverges from the [merge()] method in
#'    two main ways:
#'    - The `by.x` and `by.y` column names must be of length 1.
#'    - If variable names are shared in `x` and `y`, the spectra
#'      variables of `x` are not modified. It's only the `y`
#'      variables that are appended with the suffix defined in
#'      `suffix.y`. This is to avoid modifying any core spectra
#'      variables that would lead to an invalid object.
#'    - Duplicated Spectra keys (i.e. `x[[by.x]]`) are not
#'      allowed. Duplicated keys in the `DataFrame` (i.e `y[[by.y]]`)
#'      throw a warning and only the last occurrence is kept. These
#'      should be explored and ideally be removed using for
#'      `QFeatures::reduceDataFrame()`, `PMS::reducePSMs()` or similar
#'      functions.
#'    For a more general function that allows to append `data.frame`,
#'    `DataFrame` and `matrix` see `cbind2()`.
#'
#' - `split()`: splits the `Spectra` object based on parameter `f` into a `list`
#'   of `Spectra` objects.
#'
#' @param BPPARAM Parallel setup configuration. See [bpparam()] for more
#'     information. This is passed directly to the [backendInitialize()] method
#'     of the [MsBackend-class].
#'
#' @param by.x A `character(1)` specifying the spectra variable used
#'     for merging. Default is `"spectrumId"`.
#'
#' @param by.y A `character(1)` specifying the column used for
#'     merging. Set to `by.x` if missing.
#'
#' @param drop For `split()`: not considered.
#'
#' @param f For `split()`: factor defining how to split `x`. See [base::split()]
#'     for details.
#'     For `combineSpectra()`: `factor` defining the grouping of the spectra
#'     that should be combined. Defaults to `x$dataStorage`.
#'
#' @param FUN For `combineSpectra()`: function to combine the (peak matrices)
#'     of the spectra. Defaults to [combinePeaksData()].
#'
#' @param p For `combineSpectra()`: `factor` defining how to split the input
#'     `Spectra` for parallel processing. Defaults to `x$dataStorage`, i.e.,
#'     depending on the used backend, per-file parallel processing will be
#'     performed.
#'
#' @param suffix.y A `character(1)` specifying the suffix to be used
#'     for making the names of columns in the merged spectra variables
#'     unique. This suffix will be used to amend `names(y)`, while
#'     `spectraVariables(x)` will remain unchanged.
#'
#' @param x A `Spectra` object.
#'
#' @param y For `joinSpectraData()`: `DataFrame` with the spectra variables
#'     to join/add. For `cbind2()`: a `data.frame`, `DataFrame` or
#'     `matrix`. The number of rows and their order has to match the
#'     number of spectra in `x`, respectively their order.
#'
#' @param ... Additional arguments.
#'
#' @seealso
#'
#' - [combinePeaks()] for functions to aggregate mass peaks data.
#'
#' - [Spectra] for a general description of the `Spectra` object.
#'
#' @importFrom MsCoreUtils vapply1c
#'
#' @author Sebastian Gibb, Johannes Rainer, Laurent Gatto
#'
#' @examples
#'
#' ## Create a Spectra providing a `DataFrame` containing a MS data.
#'
#' spd <- DataFrame(msLevel = c(1L, 2L), rtime = c(1.1, 1.2))
#' spd$mz <- list(c(100, 103.2, 104.3, 106.5), c(45.6, 120.4, 190.2))
#' spd$intensity <- list(c(200, 400, 34.2, 17), c(12.3, 15.2, 6.8))
#'
#' s <- Spectra(spd)
#' s
#'
#' ## Create a second Spectra from mzML files and use the `MsBackendMzR`
#' ## on-disk backend.
#' sciex_file <- dir(system.file("sciex", package = "msdata"),
#'     full.names = TRUE)
#' sciex <- Spectra(sciex_file, backend = MsBackendMzR())
#' sciex
#'
#' ## Subset to the first 100 spectra to reduce running time of the examples
#' sciex <- sciex[1:100]
#'
#'
#' ##  --------  COMBINE SPECTRA  --------
#'
#' ## Combining the `Spectra` object `s` with the MS data from `sciex`.
#' ## Calling directly `c(s, sciex)` would result in an error because
#' ## both backends use a different backend. We thus have to first change
#' ## the backends to the same backend. We change the backend of the `sciex`
#' ## `Spectra` to a `MsBackendMemory`, the backend used by `s`.
#'
#' sciex <- setBackend(sciex, MsBackendMemory())
#'
#' ## Combine the two `Spectra`
#' all <- c(s, sciex)
#' all
#'
#' ## The new `Spectra` objects contains the union of spectra variables from
#' ## both:
#' spectraVariables(all)
#'
#' ## The spectra variables that were not present in `s`:
#' setdiff(spectraVariables(all), spectraVariables(s))
#'
#' ## The values for these were filled with missing values for spectra from
#' ## `s`:
#' all$peaksCount |> head()
#'
#'
#' ##  --------  AGGREGATE SPECTRA  --------
#'
#' ## Sets of spectra can be combined into a single, representative spectrum
#' ## per set using `combineSpectra()`. This aggregates the peaks data (i.e.
#' ## the spectra's m/z and intensity values) while using the values for all
#' ## spectra variables from the first spectrum per set. Below we define the
#' ## sets as all spectra measured in the *same second*, i.e. rounding their
#' ## retention time to the next closer integer value.
#' f <- round(rtime(sciex))
#' head(f)
#'
#' cmp <- combineSpectra(sciex, f = f)
#'
#' ## The length of `cmp` is now equal to the length of unique levels in `f`:
#' length(cmp)
#'
#' ## The spectra variable value from the first spectrum per set is used in
#' ## the representative/combined spectrum:
#' cmp$rtime
#'
#' ## The peaks data was aggregated: the number of mass peaks of the first six
#' ## spectra from the original `Spectra`:
#' lengths(sciex) |> head()
#'
#' ## and for the first aggreagated spectra:
#' lengths(cmp) |> head()
#'
#' ## The default peaks data aggregation method joins all mass peaks. See
#' ## documentation of the `combinePeaksData()` function for more options.
#'
#'
#' ##  --------  SPLITTING DATA  --------
#'
#' ## A `Spectra` can be split into a `list` of `Spectra` objects using the
#' ## `split()` function defining the sets into which the `Spectra` should
#' ## be splitted into with parameter `f`.
#' sciex_split <- split(sciex, f)
#'
#' length(sciex_split)
#' sciex_split |> head()
#'
#'
#' ##  --------  ADDING SPECTRA DATA  --------
#'
#' ## Adding new spectra variables
#' sciex1 <- filterDataOrigin(sciex, dataOrigin(sciex)[1])
#' spv <- DataFrame(spectrumId = sciex1$spectrumId[3:12], ## used for merging
#'                  var1 = rnorm(10),
#'                  var2 = sample(letters, 10))
#' spv
#'
#' sciex2 <- joinSpectraData(sciex1, spv, by.y = "spectrumId")
#'
#' spectraVariables(sciex2)
#' spectraData(sciex2)[1:13, c("spectrumId", "var1", "var2")]
#'
#' ## Append new spectra variables with cbind2()
#' df <- data.frame(cola = seq_len(length(sciex1)), colb = "b")
#' data_append <- cbind2(sciex1, df)
NULL

#' @rdname combineSpectra
#'
#' @exportMethod c
setMethod("c", "Spectra", function(x, ...) {
    .concatenate_spectra(unname(list(unname(x), ...)))
})

#' @rdname combineSpectra
#'
#' @export
setMethod("cbind2", signature(x = "Spectra",
                              y = "dataframeOrDataFrameOrmatrix"),
          function(x, y, ...) {
              x@backend <- cbind2(x@backend, y, ...)
          })

#' @rdname combineSpectra
setMethod("split", "Spectra", function(x, f, drop = FALSE, ...) {
    bcknds <- split(x@backend, f, ...)
    lapply(bcknds, function(b) {
        slot(x, "backend", check = FALSE) <- b
        x
    })
})


################################################################################
##
## Aggregating peaks data
##
################################################################################

#' @title Aggregating and combining mass peaks data
#'
#' @name combinePeaks
#'
#' @description
#'
#' In addition to aggregating content of spectra variables (describe in
#' [combineSpectra()]) it is also possible to aggregate and combine mass peaks
#' data from individual spectra within a `Spectra`. These `combinePeaks()`
#' function combines mass peaks **within each spectrum** with a difference in
#' their m/z values that is smaller than the maximal acceptable difference
#' defined by `ppm` and `tolerance`. Parameters `intensityFun` and `mzFun`
#' allow to define functions to aggregate the intensity and m/z values for
#' each such group of peaks. With `weighted = TRUE` (the default), the m/z
#' value of the combined peak is calculated using an intensity-weighted mean
#' and parameter `mzFun` is ignored. The [MsCoreUtils::group()] function is
#' used for the grouping of mass peaks. Parameter `msLevel.` allows to define
#' selected MS levels for which peaks should be combined. This function
#' returns a `Spectra` with the same number of spectra than the input object,
#' but with possibly combined peaks within each spectrum.
#' Additional peak variables (other than `"mz"` and `"intensity"`) are
#' dropped (i.e. their values are replaced with `NA`) for combined peaks
#' unless they are constant across the combined peaks. See also
#' [reduceSpectra()] for a function to select a single *representative*
#' mass peak for each peak group.
#'
#' @param intensityFun Function to aggregate intensities for all peaks in
#'     each peak group into a single intensity value.
#'
#' @param msLevel. `integer` defining the MS level(s) of the spectra to which
#'     the function should be applied (defaults to all MS levels of `object`.
#'
#' @param mzFun Function to aggregate m/z values for all mass peaks within
#'     each peak group into a single m/z value. This parameter is ignored if
#'     `weighted = TRUE` (the default).
#'
#' @param object A `Spectra` object.
#'
#' @param ppm `numeric(1)` defining a relative, m/z-dependent, maximal
#'     accepted difference between m/z values for peaks to be grouped. Default
#'     is `ppm = 20`.
#'
#' @param tolerance `numeric(1)` allowing to define a constant maximal
#'     accepted difference between m/z values for peaks to be grouped. Default
#'     is `tolerance = 0`.
#'
#' @param weighted `logical(1)` whether m/z values of peaks within each peak
#'     group should be aggregated into a single m/z value using an
#'     intensity-weighted mean. Defaults to `weighted = TRUE`.
#'
#' @param ... ignored.
#'
#' @md
#'
#' @seealso
#'
#' - [combineSpectra()] for functions to combine or aggregate `Spectra`'s
#'   spectra data.
#'
#' - [combinePeaksData()] for the function to combine the mass peaks data.
#'
#' - [reduceSpectra()] and similar functions to filter mass peaks data.
#'
#' - [Spectra] for a general description of the `Spectra` object.
#'
#' @author Sebastian Gibb, Johannes Rainer, Laurent Gatto
#'
#' @examples
#'
#' ## Create a Spectra from mzML files and use the `MsBackendMzR` on-disk
#' ## backend.
#' sciex_file <- dir(system.file("sciex", package = "msdata"),
#'     full.names = TRUE)
#' sciex <- Spectra(sciex_file, backend = MsBackendMzR())
#'
#' ## Combine mass peaks per spectrum with a difference in their m/z value
#' ## that is smaller than 20 ppm. The intensity values of such peaks are
#' ## combined by summing their values, while for the m/z values the median
#' ## is reported
#' sciex_comb <- combinePeaks(sciex, ppm = 20,
#'     intensityFun = sum, mzFun = median)
#'
#' ## Comparing the number of mass peaks before and after aggregation
#' lengths(sciex) |> head()
#' lengths(sciex_comb) |> head()
#'
#' ## Plotting the first spectrum before and after aggregation
#' par(mfrow = c(1, 2))
#' plotSpectra(sciex[2L])
#' plotSpectra(sciex_comb[2L])
#'
#' ## Using `reduceSpectra()` to keep for each group of mass peaks with a
#' ## difference in their m/z values < 20ppm the one with the highest intensity.
#' sciex_red <- reduceSpectra(sciex, ppm = 20)
#'
#' ## Comparing the number of mass peaks before and after the operation
#' lengths(sciex) |> head()
#' lengths(sciex_red) |> head()
NULL

#' @rdname hidden_aliases
setMethod("combinePeaks", "list", function(object, ...) {
    .Deprecated("combinePeaksData", old = "combinePeaks",
                msg = paste0("'combinePeaks' for lists of peak matrices is ",
                             "deprecated; please use 'combinePeaksData' ",
                             "instead."))
    combinePeaksData(object, ...)
})

#' @rdname combinePeaks
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


################################################################################
##
## Filtering, subsetting Spectra: subsetting Spectra and its data content.
##
################################################################################

#' @title Filter and subset Spectra objects
#'
#' @name filterMsLevel
#'
#' @aliases [,Spectra-method
#' @aliases filterAcquisitionNum
#' @aliases filterDataOrigin
#' @aliases filterDataStorage
#' @aliases filterEmptySpectra
#' @aliases filterIsolationWindow
#' @aliases filterMsLevel
#' @aliases filterPolarity
#' @aliases filterPrecursorCharge
#' @aliases filterPrecursorIsotopes
#' @aliases filterPrecursorMzRange
#' @aliases filterPrecursorMzValues
#' @aliases filterPrecursorScan
#' @aliases filterRanges
#' @aliases filterRt
#' @aliases filterValues
#' @aliases dropNaSpectraVariables
#' @aliases selectSpectraVariables
#' @aliases filterIntensity
#' @aliases filterMzRange
#' @aliases filterMzValues
#' @aliases reduceSpectra
#'
#' @description
#'
#' A variety of functions to filter or subset `Spectra` objects are available.
#' These can be generally separated into two main classes: I) *classical*
#' subset operations that immediately reduce the number of spectra in the
#' object and II) filters that reduce the **content** of the object without
#' changing its length (i.e. the number of spectra). The latter can be further
#' subdivided into functions that affect the content of the `spectraData` (i.e.
#' the general spectrum metadata) and those that reduce the content of the
#' object's `peaksData` (i.e. the m/z and intensity values of a spectrum's
#' mass peaks).
#'
#' A description of functions from these 3 different categories are given below
#' in sections *Subset `Spectra`*, *Filter content of `spectraData()`* and
#' *Filter content of `peaksData()`*, respectively.
#'
#'
#' @section Subset `Spectra`:
#'
#' These functions affect the number of spectra in a `Spectra` object creating
#' a subset of the original object without affecting its content.
#'
#' - `[`: subsets the spectra keeping only selected elements (`i`). The method
#'   **always** returns a `Spectra` object.
#'
#' - `filterAcquisitionNum()`: filters the object keeping only spectra matching
#'   the provided acquisition numbers (argument `n`). If `dataOrigin` or
#'   `dataStorage` is also provided, `object` is subsetted to the spectra with
#'   an acquisition number equal to `n` **in spectra with matching dataOrigin
#'   or dataStorage values** retaining all other spectra.
#'   Returns the filtered `Spectra`.
#'
#' - `filterDataOrigin()`: filters the object retaining spectra matching the
#'   provided `dataOrigin`. Parameter `dataOrigin` has to be of type
#'   `character` and needs to match exactly the data origin value of the
#'   spectra to subset.
#'   Returns the filtered `Spectra` object (with spectra ordered according to
#'   the provided `dataOrigin` parameter).
#'
#' - `filterDataStorage()`: filters the object retaining spectra stored in the
#'   specified `dataStorage`. Parameter `dataStorage` has to be of type
#'   `character` and needs to match exactly the data storage value of the
#'   spectra to subset.
#'   Returns the filtered `Spectra` object (with spectra ordered according to
#'   the provided `dataStorage` parameter).
#'
#' - `filterEmptySpectra()`: removes empty spectra (i.e. spectra without peaks).
#'   Returns the filtered `Spectra` object (with spectra in their
#'   original order).
#'
#' - `filterIsolationWindow()`: retains spectra that contain `mz` in their
#'   isolation window m/z range (i.e. with an `isolationWindowLowerMz` <= `mz`
#'   and `isolationWindowUpperMz` >= `mz`. Returns the filtered `Spectra`
#'   object (with spectra in their original order).
#'
#' - `filterMsLevel()`: filters object by MS level keeping only spectra matching
#'   the MS level specified with argument `msLevel`. Returns the filtered
#'   `Spectra` (with spectra in their original order).
#'
#' - `filterPolarity()`: filters the object keeping only spectra matching the
#'   provided polarity. Returns the filtered `Spectra` (with spectra in their
#'   original order).
#'
#' - `filterPrecursorCharge()`: retains spectra with the defined precursor
#'   charge(s).
#'
#' - `filterPrecursorIsotopes()`: groups MS2 spectra based on their precursor
#'   m/z and precursor intensity into predicted isotope groups and keep for each
#'   only the spectrum representing the monoisotopic precursor. MS1 spectra
#'   are returned as is. See documentation for `deisotopeSpectra()` below for
#'   details on isotope prediction and parameter description.
#'
#' - `filterPrecursorMaxIntensity()`: filters the `Spectra` keeping for groups
#'   of (MS2) spectra with similar precursor m/z values (given parameters
#'   `ppm` and `tolerance`) the one with the highest precursor intensity. The
#'   function filters only MS2 spectra and returns all MS1 spectra. If
#'   precursor intensities are `NA` for all spectra within a spectra group, the
#'   first spectrum of that groups is returned.
#'   Note: some manufacturers don't provide precursor intensities. These can
#'   however also be estimated with [estimatePrecursorIntensity()].
#'
#' - `filterPrecursorMzRange()` (previously `filterPrecursorMz()` which is now
#'   deprecated): retains spectra with a precursor m/z within the
#'   provided m/z range. See examples for details on selecting spectra with
#'   a precursor m/z for a target m/z accepting a small difference in *ppm*.
#'
#' - `filterPrecursorMzValues()`: retains spectra with precursor m/z matching
#'   any of the provided m/z values (given `ppm` and `tolerance`). Spectra with
#'   missing precursor m/z value (e.g. MS1 spectra) are dropped.
#'
#' - `filterPrecursorScan()`: retains parent (e.g. MS1) and children scans (e.g.
#'   MS2) of acquisition number `acquisitionNum`. Returns the filtered
#'   `Spectra` (with spectra in their original order). Parameter `f` allows to
#'   define which spectra belong to the same sample or original data file (
#'   defaults to `f = dataOrigin(object)`).
#'
#' - `filterRanges()`: allows filtering of the `Spectra` object based on user
#'   defined *numeric* ranges (parameter `ranges`) for one or more available
#'   spectra variables in object (spectra variable names can be specified with
#'   parameter `spectraVariables`). Spectra for which the value of a spectra
#'   variable is within it's defined range are retained. If multiple
#'   ranges/spectra variables are defined, the `match` parameter can be used
#'   to specify whether all conditions (`match = "all"`; the default) or if
#'   any of the conditions must match (`match = "any"`; all spectra for which
#'   values are within any of the provided ranges are retained).
#'
#' - `filterRt()`: retains spectra of MS level `msLevel` with retention
#'   times (in seconds) within (`>=`) `rt[1]` and (`<=`)
#'   `rt[2]`. Returns the filtered `Spectra` (with spectra in their
#'   original order).
#'
#' - `filterValues()`: allows filtering of the `Spectra` object based on
#'   similarities of *numeric* values of one or more `spectraVariables(object)`
#'   (parameter `spectraVariables`) to provided values (parameter `values`)
#'   given acceptable differences (parameters tolerance and ppm). If multiple
#'   values/spectra variables are defined, the `match` parameter can be used
#'   to specify whether all conditions (`match = "all"`; the default) or if
#'   any of the conditions must match (`match = "any"`; all spectra for which
#'   values are within any of the provided ranges are retained).
#'
#'
#' @section Filter content of `spectraData()`:
#'
#' The functions described in this section filter the content from a
#' `Spectra`'s spectra data, i.e. affect values of, or complete, spectra
#' variables. None of these functions reduces the object's number of spectra.
#'
#' - `dropNaSpectraVariables()`: removes spectra variables (i.e. columns in the
#'   object's `spectraData` that contain only missing values (`NA`). Note that
#'   while columns with only `NA`s are removed, a `spectraData()` call after
#'   `dropNaSpectraVariables()` might still show columns containing `NA` values
#'   for *core* spectra variables. The total number of spectra is not changed
#'   by this function.
#'
#' - `selectSpectraVariables()`: reduces the information within the object to
#'   the selected spectra variables: all data for variables not specified will
#'   be dropped. For mandatory columns (i.e., those listed by
#'   [coreSpectraVariables()], such as *msLevel*, *rtime* ...) only
#'   the values will be dropped but not the variable itself. Additional (or
#'   user defined) spectra variables will be completely removed.
#'   Returns the filtered `Spectra`.
#'
#'
#' @section Filter content of `peaksData()`:
#'
#' The functions described in this section filter the content of the
#' `Spectra`'s peaks data, i.e. either the number or the values (*m/z* or
#' intensity values) of the mass peaks. Also, the actual operation is only
#' executed once peaks data is accessed (through `peaksData()`,
#' `mz()` or `intensity()`) or `applyProcessing()` is called.
#' These operations don't affect the number of spectra in the `Spectra` object.
#'
#' - `deisotopeSpectra()`: *deisotopes* each spectrum keeping only the
#'   monoisotopic peak for groups of isotopologues. Isotopologues are
#'   estimated using the [isotopologues()] function from the
#'   *MetaboCoreUtils* package. Note that
#'   the default parameters for isotope prediction/detection have been
#'   determined using data from the Human Metabolome Database (HMDB) and
#'   isotopes for elements other than CHNOPS might not be detected. See
#'   parameter `substDefinition` in the documentation of [isotopologues()] for
#'   more information. The approach and code to define the parameters for
#'   isotope prediction is described
#'   [here](https://github.com/EuracBiomedicalResearch/isotopologues).
#'
#' - `filterFourierTransformArtefacts()`: removes (Orbitrap) fast fourier
#'   artefact peaks from spectra (see examples below). The function iterates
#'   through all intensity ordered peaks in a spectrum and removes all peaks
#'   with an m/z within +/- `halfWindowSize` of the current peak if their
#'   intensity is lower than `threshold` times the current peak's intensity.
#'   Additional parameters `keepIsotopes`, `maxCharge` and `isotopeTolerance`
#'   allow to avoid removing of potential `[13]C` isotope peaks (`maxCharge`
#'   being the maximum charge that should be considered and `isotopeTolerance`
#'   the absolute acceptable tolerance for matching their m/z).
#'   See [filterFourierTransformArtefacts()] for details and background and
#'   `deisitopeSpectra()` for an alternative.
#'
#' - `filterIntensity()`: filters mass peaks in each spectrum keeping only
#'   those with intensities that are within the provided range or match the
#'   criteria of the provided function. For the former, parameter `intensity`
#'   has to be a `numeric` defining the intensity range, for the latter a
#'   `function` that takes the intensity values of the spectrum and returns
#'   a `logical` whether the peak should be retained or not (see examples
#'   below for details) - additional parameters to the function can be passed
#'   with `...`.
#'   To remove only peaks with intensities below a certain threshold, say
#'   100, use `intensity = c(100, Inf)`. Note: also a single value can be
#'   passed with the `intensity` parameter in which case an upper limit of
#'   `Inf` is used.
#'   Note that this function removes also peaks with missing intensities
#'   (i.e. an intensity of `NA`). Parameter `msLevel.` allows to restrict the
#'   filtering to spectra of the specified MS level(s).
#'
#' - `filterMzRange()`: filters mass peaks in the object keeping or removing
#'   those in each spectrum that are within the provided m/z range. Whether
#'   peaks are retained or removed can be configured with parameter `keep`
#'   (default `keep = TRUE`).
#'
#' - `filterMzValues()`: filters mass peaks in the object keeping all
#'   peaks in each spectrum that match the provided m/z value(s) (for
#'   `keep = TRUE`, the default) or removing all of them (for `keep = FALSE`).
#'   The m/z matching considers also the absolute `tolerance` and m/z-relative
#'   `ppm` values. `tolerance` and `ppm` have to be of length 1.
#'
#' - `filterPeaksRanges()`: filters mass peaks of a `Spectra` object using any
#'   set of range-based filters on numeric spectra or peaks variables. See
#'   [filterPeaksRanges()] for more information.
#'
#' - `filterPrecursorPeaks()`: removes peaks from each spectrum in `object` with
#'   an m/z equal or larger than the m/z of the precursor, depending on the
#'   value of parameter `mz`: for `mz = ==" (the default) peaks with matching
#'   m/z (considering an absolute and relative acceptable difference depending
#'   on `tolerance` and `ppm`, respectively) are removed. For `mz = ">="` all
#'   peaks with an m/z larger or equal to the precursor m/z (minus `tolerance`
#'   and the `ppm` of the precursor m/z) are removed. Parameter `msLevel.`
#'   allows to restrict the filter to certain MS levels (by default the filter
#'   is applied to all MS levels). Note that no peaks are removed if the
#'   precursor m/z is `NA` (e.g. typically for MS1 spectra).
#'
#' - `reduceSpectra()`: keeps for groups of peaks with similar m/z values in
#'   (given `ppm` and `tolerance`) in each spectrum only the mass peak with the
#'   highest intensity removing all other peaks hence *reducing* each
#'   spectrum to the highest intensity peaks per *peak group*.
#'   Peak groups are defined using the [group()] function from the
#'   *MsCoreUtils* package. See also the [combinePeaks()] function for an
#'   alternative function to combine peaks within each spectrum.
#'
#' @param acquisitionNum for `filterPrecursorScan()`: `integer` with the
#'     acquisition number of the spectra to which the object should be
#'     subsetted.
#'
#' @param charge For `deisotopeSpectra()`: expected charge of the ionized
#'     compounds. See [isotopologues()] for details.
#'
#' @param dataOrigin For `filterDataOrigin()`: `character` to define which
#'     spectra to keep.
#'     For `filterAcquisitionNum()`: optionally specify if filtering should
#'     occurr only for spectra of selected `dataOrigin`.
#'
#' @param dataStorage For `filterDataStorage()`: `character` to define which
#'     spectra to keep.
#'     For `filterAcquisitionNum()`: optionally specify if filtering should
#'     occur only for spectra of selected `dataStorage`.
#'
#' @param drop For `[`: not considered.
#'
#' @param f For `filterPrecursorScan()`: defining which spectra
#'     belong to the same original data file (sample): Defaults to
#'     `f = dataOrigin(x)`.
#'
#' @param halfWindowSize For `filterFourierTransformArtefacts()`: `numeric(1)`
#'     defining the m/z window left and right of a peak where to remove
#'     fourier transform artefacts.
#'
#' @param i For `[`: `integer`, `logical` or `character` to subset the
#'     object.
#'
#' @param intensity For `filterIntensity()`: `numeric` of length 1 or 2
#'     defining either the lower or the lower and upper intensity limit for the
#'     filtering, or a `function` that takes the intensities as input and
#'     returns a `logical` (same length then peaks in the spectrum) whether the
#'     peak should be retained or not. Defaults to `intensity = c(0, Inf)` thus
#'     only peaks with `NA` intensity are removed.
#'
#' @param isotopeTolerance For `filterFourierTransformArtefacts()`: the m/z
#'     `tolerance` to be used to define whether peaks might be isotopes of
#'     the current tested peak.
#'
#' @param j For `[`: not supported.
#'
#' @param keep For `filterMzValues()` and `filterMzRange()`: `logical(1)`
#'     whether  the matching peaks should be retained (`keep = TRUE`, the
#'     default) or dropped (`keep = FALSE`).
#'
#' @param keepIsotopes For `filterFourierTransformArtefacts()`: whether isotope
#'     peaks should not be removed as fourier artefacts.
#'
#' @param match For `filterRanges()` and `filterValues()`: `character(1) `
#'     defining whether the condition has to match for all provided
#'     `ranges`/`values` (`match = "all"`; the default), or for any of them
#'     (`match = "any"`) for spectra to be retained.
#'
#' @param maxCharge For `filterFourierTransformArtefacts()`: the maximum charge
#'     to be considered for isotopes.
#'
#' @param msLevel. `integer` defining the MS level(s) of the spectra to which
#'     the function should be applied (defaults to all MS levels of `object`.
#'     For `filterMsLevel()`: the MS level to which `object` should be
#'     subsetted.
#'
#' @param mz For `filterIsolationWindow()`: `numeric(1)` with the m/z value to
#'     filter the object. For `filterPrecursorMz()` and `filterMzRange()`:
#'     `numeric(2)` defining the lower and upper m/z boundary.
#'     For `filterMzValues()` and `filterPrecursorMzValues()`: `numeric` with
#'     the m/z values to match peaks or precursor m/z against.
#'     For `filterPrecursorPeaks()`: `character(1)` defining whether mass peaks
#'     with an m/z matching the spectrum's precursor m/z (`mz = "=="`,
#'     the default) or mass peaks with a m/z that is equal or larger
#'     (`mz = ">="`) should be removed.
#'
#' @param n for `filterAcquisitionNum()`: `integer` with the acquisition
#'     numbers to filter for.
#'
#' @param object `Spectra` object.
#'
#' @param polarity for `filterPolarity()`: `integer` specifying the polarity to
#'     to subset `object`.
#'
#' @param ppm For `filterMzValues()` and `reduceSpectra()`: `numeric(1)`
#'     defining a relative, m/z-dependent, maximal accepted difference between
#'     m/z values for peaks to be matched (or grouped).
#'     For `filterPrecursorMaxIntensity()`: `numeric(1)` defining the relative
#'     maximal accepted difference of precursor m/z values of spectra for
#'     grouping them into *precursor groups*. For `filterPrecursorIsotopes()`:
#'     passed directly to the [isotopologues()] function.
#'     For `filterValues()`: `numeric` of any length allowing to define
#'     a maximal accepted difference between user input `values` and the
#'     `spectraVariables` values. If it is not equal to the length of the
#'     value provided with parameter `spectraVariables`, `ppm[1]` will be
#'     recycled.
#'
#' @param ranges for `filterRanges()`: A `numeric` vector of paired values
#'     (upper and lower boundary) that define the ranges to filter the `object`.
#'     These paired values need to be in the same order as the
#'     `spectraVariables` parameter (see below).
#'
#' @param rt for `filterRt()`: `numeric(2)` defining the retention time range to
#'     be used to subset/filter `object`.
#'
#' @param spectraVariables For `selectSpectraVariables()`: `character` with the
#'     names of the spectra variables to which the backend should be
#'     subsetted. For `filterRanges()` and `filterValues()`: `character`
#'     vector specifying the column(s) from `spectraData(object)` on which
#'     to filter the data and that correspond to the the names of the
#'     spectra variables that should be used for the filtering.
#'
#' @param substDefinition For `deisotopeSpectra()` and
#'     `filterPrecursorIsotopes()`: `matrix` or `data.frame` with definitions
#'     of isotopic substitutions. Uses by default isotopic substitutions
#'     defined from all compounds in the Human Metabolome Database (HMDB). See
#'     [isotopologues()] or [isotopicSubstitutionMatrix()] in the
#'     *MetaboCoreUtils* for details.
#'
#' @param threshold For `filterFourierTransformArtefacts()`: the relative
#'     intensity (to a peak) below which peaks are considered fourier
#'     artefacts. Defaults to `threshold = 0.2` hence removing peaks that
#'     have an intensity below 0.2 times the intensity of the tested peak
#'     (within the selected `halfWindowSize`).
#'
#' @param tolerance For `filterMzValues()` and `reduceSpectra()`:
#'     `numeric(1)` allowing to define a constant maximal accepted difference
#'     between m/z values for peaks to be matched (or grouped). For
#'     `containsMz()` it can also be of length equal `mz` to specify a different
#'     tolerance for each m/z value.
#'     For `filterPrecursorMaxIntensity()`: `numeric(1)` defining the
#'     (constant) maximal accepted difference of precursor m/z values of
#'     spectra for grouping them into *precursor groups*. For
#'     `filterPrecursorIsotopes()`: passed directly to the [isotopologues()]
#'     function. For `filterValues()`: `numeric` of any length allowing to
#'     define a maximal accepted difference between user input `values` and the
#'     `spectraVariables` values. If it is not equal to the length of the
#'     value provided with parameter `spectraVariables`, `tolerance[1]` will be
#'     recycled. Default is `tolerance = 0`.
#'
#' @param values for `filterValues()`: A `numeric` vector that define the
#'     values to filter the Spectra data. These values need to be in the same
#'     order as the `spectraVariables` parameter.
#'
#' @param x `Spectra` object.
#'
#' @param z For `filterPrecursorCharge()`: `integer()` with the precursor
#'     charges to be used as filter.
#'
#' @param ... Additional arguments.
#'
#' @seealso
#'
#' - [combineSpectra()] for functions to combine or aggregate `Spectra`.
#'
#' - [combinePeaks()] for functions to combine or aggregate a `Spectra`'s
#'   `peaksData()`
#'
#' @md
#'
#' @author Sebastian Gibb, Johannes Rainer, Laurent Gatto, Philippine Louail, Nir Shahaf
#'
#' @examples
#'
#' ## Load a `Spectra` object with LC-MS/MS data.
#' fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML",
#'     package = "msdata")
#' sps_dda <- Spectra(fl)
#' sps_dda
#'
#'
#' ##  --------  SUBSET SPECTRA  --------
#'
#' ## Subset to the first 3 spectra
#' tmp <- sps_dda[1:3]
#' tmp
#' length(tmp)
#'
#' ## Subset to all MS2 spectra; this could be done with [, or, more
#' ## efficiently, with the `filterMsLevel` function:
#' sps_dda[msLevel(sps_dda) == 2L]
#' filterMsLevel(sps_dda, 2L)
#'
#' ## Filter the object keeping only MS2 spectra with an precursor m/z value
#' ## between a specified range:
#' filterPrecursorMzRange(sps_dda, c(80, 90))
#'
#' ## Filter the object to MS2 spectra with an precursor m/z matching a
#' ## pre-defined value (given ppm and tolerance)
#' filterPrecursorMzValues(sps_dda, 85, ppm = 5, tolerance = 0.1)
#'
#' ## The `filterRanges()` function allows to filter a `Spectra` based on
#' ## numerical ranges of any of its (numerical) spectra variables.
#' ## First, determine the variable(s) on which to base the filtering:
#' sv <- c("rtime", "precursorMz", "peaksCount")
#' ## Note that ANY variables can be chosen here, and as many as wanted.
#'
#' ## Define the ranges (pairs of values with lower and upper boundary) to be
#' ## used for the individual spectra variables. The first two values will be
#' ## used for the first spectra variable (e.g., `"rtime"` here), the next two
#' ## for the second (e.g. `"precursorMz"` here) and so on:
#' ranges <- c(30, 350, 200, 500, 350, 600)
#'
#' ## Input the parameters within the filterRanges function:
#' filt_spectra <- filterRanges(sps_dda, spectraVariables = sv,
#'                 ranges = ranges)
#' filt_spectra
#'
#' ## `filterRanges()` can also be used to filter a `Spectra` object with
#' ## multiple ranges for the same `spectraVariable` (e.g, here `"rtime"`)
#' sv <- c("rtime", "rtime")
#' ranges <- c(30, 100, 200, 300)
#' filt_spectra <- filterRanges(sps_dda, spectraVariables = sv,
#'                 ranges = ranges, match = "any")
#' filt_spectra
#'
#' ## While `filterRanges()` filtered on numeric ranges, `filterValues()`
#' ## allows to filter an object matching spectra variable values to user
#' ## provided values (allowing to configure allowed differences using the
#' ## `ppm` and `tolerance` parameters).
#' ## First determine the variable(s) on which to base the filtering:
#' sv <- c("rtime", "precursorMz")
#' ## Note that ANY variables can be chosen here, and as many as wanted.
#'
#' ## Define the values that will be used to filter the spectra based on their
#' ## similarities to their respective `spectraVariables`.
#' ## The first values in the parameters values, tolerance and ppm will be
#' ## used for the first spectra variable (e.g. `"rtime"` here), the next for
#' ## the second (e.g. `"precursorMz"` here) and so on:
#' values <- c(350, 80)
#' tolerance <- c(100, 0.1)
#' ppm <- c(0, 50)
#'
#' ## Input the parameters within the `filterValues()` function:
#' filt_spectra <- filterValues(sps_dda, spectraVariables = sv,
#'                 values = values, tolerance = tolerance, ppm = ppm)
#' filt_spectra
#'
#'
#' ##  --------  FILTER SPECTRA DATA  --------
#'
#' ## Remove spectra variables without content (i.e. with only missing values)
#' sps_noNA <- dropNaSpectraVariables(sps_dda)
#'
#' ## This reduced the size of the object slightly
#' print(object.size(sps_dda), unit = "MB")
#' print(object.size(sps_noNA), unit = "MB")
#'
#' ## With the `selectSpectraVariables()` function it is in addition possible
#' ## to subset the data of a `Spectra` to the selected columns/variables,
#' ## keeping only their data:
#' tmp <- selectSpectraVariables(sps_dda, c("msLevel", "mz", "intensity",
#'     "scanIndex"))
#' print(object.size(tmp), units = "MB")
#'
#' ## Except the selected variables, all data is now removed. Accessing
#' ## core spectra variables still works, but returns only NA
#' rtime(tmp) |> head()
#'
#'
#' ##  --------  FILTER PEAKS DATA  --------
#'
#' ## `filterMzValues()` filters the mass peaks data of a `Spectra` retaining
#' ## only those mass peaks with an m/z value matching the provided value(s).
#' sps_sub <- filterMzValues(sps_dda, mz = c(103, 104), tolerance = 0.3)
#'
#' ## The filtered `Spectra` has the same length
#' length(sps_dda)
#' length(sps_sub)
#'
#' ## But the number of mass peaks changed
#' lengths(sps_dda) |> head()
#' lengths(sps_sub) |> head()
#'
#' ## This function can also be used to remove specific peaks from a spectrum
#' ## by setting `keep = FALSE`.
#' sps_sub <- filterMzValues(sps_dda, mz = c(103, 104),
#'     tolerance = 0.3, keep = FALSE)
#' lengths(sps_sub) |> head()
#'
#' ## With the `filterMzRange()` function it is possible to keep (or remove)
#' ## mass peaks with m/z values within a specified numeric range.
#' sps_sub <- filterMzRange(sps_dda, mz = c(100, 150))
#' lengths(sps_sub) |> head()
#'
#' ## See also the `filterPeaksRanges()` function for a more flexible framework
#' ## to filter mass peaks
#'
#'
#' ## Removing fourier transform artefacts seen in Orbitra data.
#'
#' ## Loading an Orbitrap spectrum with artefacts.
#' data(fft_spectrum)
#' plotSpectra(fft_spectrum, xlim = c(264.5, 265.5))
#' plotSpectra(fft_spectrum, xlim = c(264.5, 265.5), ylim = c(0, 5e6))
#'
#' fft_spectrum <- filterFourierTransformArtefacts(fft_spectrum)
#' fft_spectrum
#' plotSpectra(fft_spectrum, xlim = c(264.5, 265.5), ylim = c(0, 5e6))
#'
#' ## Using a few examples peaks in your data you can optimize the parameters
#' fft_spectrum_filtered <- filterFourierTransformArtefacts(fft_spectrum,
#'                                                halfWindowSize = 0.2,
#'                                                threshold = 0.005,
#'                                                keepIsotopes = TRUE,
#'                                                maxCharge = 5,
#'                                                isotopeTolerance = 0.005
#'                                                )
#'
#' fft_spectrum_filtered
#' length(mz(fft_spectrum_filtered)[[1]])
#' plotSpectra(fft_spectrum_filtered, xlim = c(264.5, 265.5), ylim = c(0, 5e6))
#'
#'
#' ## *Reducing* a `Spectra` keeping for groups of mass peaks (characterized
#' ## by similarity of their m/z values) only one representative peak. This
#' ## function helps cleaning fragment spectra.
#' ## Filter the data set to MS2 spectra
#' ms2 <- filterMsLevel(sps_dda, 2L)
#'
#' ## For groups of fragment peaks with a difference in m/z < 0.1, keep only
#' ## the largest one.
#' ms2_red <- reduceSpectra(ms2, ppm = 0, tolerance = 0.1)
#' lengths(ms2) |> tail()
#' lengths(ms2_red) |> tail()
NULL

#' @rdname filterMsLevel
setMethod("dropNaSpectraVariables", "Spectra", function(object) {
    object@backend <- dropNaSpectraVariables(object@backend)
    object
})

#' @rdname filterMsLevel
setMethod(
    "selectSpectraVariables", "Spectra",
    function(object, spectraVariables = union(spectraVariables(object),
                                              peaksVariables(object))) {
        spectraVariables <- union(spectraVariables, "dataStorage")
        object@backend <- selectSpectraVariables(
            object@backend, spectraVariables = spectraVariables)
        object
    })

#' @rdname filterMsLevel
#'
#' @export
setMethod("[", "Spectra", function(x, i, j, ..., drop = FALSE) {
    if (!missing(j))
        stop("Subsetting 'Spectra' by columns is not (yet) supported")
    if (missing(i))
        return(x)
    slot(x, "backend", check = FALSE) <- extractByIndex(
        x@backend, i2index(i, length(x)))
    x
})

#' @rdname filterMsLevel
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

#' @rdname filterMsLevel
setMethod("filterEmptySpectra", "Spectra", function(object) {
    object@backend <- extractByIndex(object@backend,
                                     which(as.logical(lengths(object))))
    object@processing <- .logging(object@processing,
                                  "Filter: removed empty spectra.")
    object
})

#' @rdname filterMsLevel
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

#' @rdname filterMsLevel
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

#' @rdname filterMsLevel
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

#' @rdname filterMsLevel
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


#' @rdname filterMsLevel
setMethod("filterIsolationWindow", "Spectra", function(object, mz = numeric()) {
    object@backend <- filterIsolationWindow(object@backend, mz = mz)
    object@processing <- .logging(object@processing,
                                  "Filter: select spectra containing m/z ",
                                  mz, " in their isolation window")
    object
})

#' @rdname filterMsLevel
setMethod("filterMsLevel", "Spectra", function(object, msLevel. = integer()) {
    object@backend <- filterMsLevel(object@backend, msLevel = msLevel.)
    object@processing <- .logging(object@processing,
                                  "Filter: select MS level(s) ",
                                  paste0(unique(msLevel.), collapse = " "))
    object
})

#' @rdname filterMsLevel
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

#' @rdname filterMsLevel
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

#' @rdname filterMsLevel
setMethod("filterPolarity", "Spectra", function(object, polarity = integer()) {
    object@backend <- filterPolarity(object@backend, polarity = polarity)
    object@processing <- .logging(object@processing,
                                  "Filter: select spectra with polarity ",
                                  paste0(polarity, collapse = " "))
    object
})

#' @rdname filterMsLevel
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

#' @rdname filterMsLevel
setMethod("filterPrecursorMzRange", "Spectra",
          function(object, mz = numeric()) {
              object@backend <- filterPrecursorMzRange(object@backend, mz)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select spectra with a precursor m/z within [",
                  paste0(mz, collapse = ", "), "]")
              object
          })

#' @rdname filterMsLevel
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

#' @rdname filterMsLevel
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

#' @rdname filterMsLevel
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

#' @rdname filterMsLevel
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

#' @rdname filterMsLevel
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

#' @rdname filterMsLevel
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


################################################################################
##
## Data manipulation and analysis operations (lazy processing)
##
################################################################################

#' @title Data manipulation and analysis methods
#'
#' @name addProcessing
#'
#' @aliases addProcessing
#' @aliases applyProcessing
#' @aliases bin
#' @aliases containsMz
#' @aliases containsNeutralLoss
#' @aliases entropy
#' @aliases pickPeaks
#' @aliases replaceIntensitiesBelow
#' @aliases reset
#' @aliases smooth
#' @aliases spectrapply
#'
#' @description
#'
#' Various data analysis functions are available for `Spectra` objects. These
#' can be categorized into functions that either return a `Spectra` object
#' (with the manipulated data) and functions that directly return the
#' result from the calculation. For the former category, the data manipulations
#' are cached in the result object's *processing queue* and only exectuted
#' on-the-fly when the respective data gets extracted from the `Spectra` (see
#' section *The processing queue* for more information).
#'
#' For the second category, the calculations are directly executed and the
#' result, usually one value per spectrum, returned. Generally, to reduce
#' memory demand, a chunk-wise processing of the data is performed.
#'
#'
#' @section Data analysis methods returning a `Spectra`:
#'
#' The methods listed here return a `Spectra` object as a result.
#'
#' - `addProcessing()`: adds an arbitrary function that should be applied to the
#'   peaks matrix of every spectrum in `object`. The function (can be passed
#'   with parameter `FUN`) is expected to take a peaks matrix as input and to
#'   return a peaks matrix. A peaks matrix is a numeric matrix with two columns,
#'   the first containing the m/z values of the peaks and the second the
#'   corresponding intensities. The function has to have `...` in its
#'   definition. Additional arguments can be passed with `...`. With parameter
#'   `spectraVariables` it is possible to define additional spectra variables
#'   from `object` that should be passed to the function `FUN`. These will be
#'   passed by their name (e.g. specifying `spectraVariables = "precursorMz"`
#'   will pass the spectra's precursor m/z as a parameter named `precursorMz`
#'   to the function. The only exception is the spectra's MS level, these will
#'   be passed to the function as a parameter called `spectrumMsLevel` (i.e.
#'   with `spectraVariables = "msLevel"` the MS levels of each spectrum will be
#'   submitted to the function as a parameter called `spectrumMsLevel`).
#'   Examples are provided in the package vignette.
#'
#' - `bin()`: aggregates individual spectra into discrete (m/z) bins. Binning is
#'   performed only on spectra of the specified MS level(s) (parameter
#'   `msLevel`, by default all MS levels of `x`). The bins can be defined with
#'   parameter `breaks` which by default are equally sized bins, with size
#'   being defined by parameter `binSize`, from the minimal to the maximal m/z
#'   of all spectra (of MS level `msLevel`) within `x`. The same bins are used
#'   for all spectra in `x`. All intensity values for peaks falling into the
#'   same bin are aggregated using the function provided with parameter `FUN`
#'   (defaults to `FUN = sum`, i.e. all intensities are summed up). Note that
#'   the binning operation is applied to the peak data on-the-fly upon data
#'   access and it is possible to *revert* the operation with the `reset()`
#'   function (see description of `reset()` below).
#'
#' - `countIdentifications`: counts the number of identifications each scan has
#'   led to. See [countIdentifications()] for more details.
#'
#' - `pickPeaks()`: picks peaks on individual spectra using a moving
#'   window-based approach (window size = `2 * halfWindowSize`). For noisy
#'   spectra there are currently two different noise estimators available,
#'   the *M*edian *A*bsolute *D*eviation (`method = "MAD"`) and
#'   Friedman's Super Smoother (`method = "SuperSmoother"`),
#'   as implemented in the [`MsCoreUtils::noise()`].
#'   The method supports also to optionally *refine* the m/z value of
#'   the identified centroids by considering data points that belong (most
#'   likely) to the same mass peak. Therefore the m/z value is calculated as an
#'   intensity weighted average of the m/z values within the peak region.
#'   The peak region is defined as the m/z values (and their respective
#'   intensities) of the `2 * k` closest signals to the centroid or the closest
#'   valleys (`descending = TRUE`) in the `2 * k` region. For the latter the `k`
#'   has to be chosen general larger. See [`MsCoreUtils::refineCentroids()`] for
#'   details.
#'   If the ratio of the signal to the highest intensity of the peak is below
#'   `threshold` it will be ignored for the weighted average.
#'
#' - `replaceIntensitiesBelow()`: replaces intensities below a specified
#'   threshold with the provided `value`. Parameter `threshold` can be either
#'   a single numeric value or a function which is applied to all non-`NA`
#'   intensities of each spectrum to determine a threshold value for each
#'   spectrum. The default is `threshold = min` which replaces all values
#'   which are <= the minimum intensity in a spectrum with `value` (the
#'   default for `value` is `0`). Note that the function specified with
#'   `threshold` is expected to have a parameter `na.rm` since `na.rm = TRUE`
#'   will be passed to the function. If the spectrum is in profile mode,
#'   ranges of successive non-0 peaks <= `threshold` are set to 0.
#'   Parameter `msLevel.` allows to apply this to only spectra of certain MS
#'   level(s).
#'
#' - `scalePeaks()`: scales intensities of peaks within each spectrum depending
#'   on parameter `by`. With `by = sum` (the default) peak intensities are
#'   divided by the sum of peak intensities within each spectrum. The sum of
#'   intensities is thus 1 for each spectrum after scaling. Parameter
#'   `msLevel.` allows to apply the scaling of spectra of a certain MS level.
#'   By default (`msLevel. = uniqueMsLevels(x)`) intensities for all
#'   spectra will be scaled.
#'
#' - `smooth()`: smooths individual spectra using a moving window-based approach
#'   (window size = `2 * halfWindowSize`). Currently, the
#'   Moving-Average- (`method = "MovingAverage"`),
#'   Weighted-Moving-Average- (`method = "WeightedMovingAverage")`,
#'   weights depending on the distance of the center and calculated
#'   `1/2^(-halfWindowSize:halfWindowSize)`) and
#'   Savitzky-Golay-Smoothing (`method = "SavitzkyGolay"`) are supported.
#'   For details how to choose the correct `halfWindowSize` please see
#'   [`MsCoreUtils::smooth()`].
#'
#'
#' @section Data analysis methods returning the result from the calculation:
#'
#' The functions listed in this section return immediately the result from the
#' calculation. To reduce memory demand (and allow parallel processing) the
#' calculations a chunk-wise processing is generally performed.
#'
#' - `chunkapply()`: apply an arbitrary function to chunks of spectra. See
#'   [chunkapply()] for details and examples.
#'
#' - `containsMz()`: checks for each of the spectra whether they contain mass
#'   peaks with an m/z equal to `mz` (given acceptable difference as defined by
#'   parameters `tolerance` and `ppm` - see [common()] for details). Parameter
#'   `which` allows to define whether any (`which = "any"`, the default) or
#'   all (`which = "all"`) of the `mz` have to match. The function returns
#'   `NA` if `mz` is of length 0 or is `NA`.
#'
#' - `containsNeutralLoss()`: checks for each spectrum in `object` if it has a
#'   peak with an m/z value equal to its precursor m/z - `neutralLoss` (given
#'   acceptable difference as defined by parameters `tolerance` and `ppm`).
#'   Returns `NA` for MS1 spectra (or spectra without a precursor m/z).
#'
#' - `entropy()`: calculates the entropy of each spectra based on the metrics
#'   suggested by Li et al. (https://doi.org/10.1038/s41592-021-01331-z).
#'   See also [nentropy()] in the *MsCoreUtils* package for details.
#'
#' - `estimatePrecursorIntensity()`: defines the precursor intensities for MS2
#'   spectra using the intensity of the matching MS1 peak from the
#'   closest MS1 spectrum (i.e. the last MS1 spectrum measured before the
#'   respective MS2 spectrum). With `method = "interpolation"` it is also
#'   possible to calculate the precursor intensity based on an interpolation of
#'   intensity values (and retention times) of the matching MS1 peaks from the
#'   previous and next MS1 spectrum. See [estimatePrecursorIntensity()] for
#'   examples and more details.
#'
#' - `estimatePrecursorMz()`: **for DDA data**: allows to estimate a fragment
#'   spectra's precursor m/z based on the reported precursor m/z and the data
#'   from the previous MS1 spectrum. See [estimatePrecursorMz()] for details.
#'
#' - `neutralLoss()`: calculates neutral loss spectra for fragment spectra. See
#'   [neutralLoss()] for detailed documentation.
#'
#' - `spectrapply()`: applies a given function to each individual spectrum or
#'   sets of a `Spectra` object. By default, the `Spectra` is split into
#'   individual spectra (i.e. `Spectra` of length 1) and the function `FUN`
#'   is applied to each of them. An alternative splitting can be defined with
#'   parameter `f`. Parameters for `FUN` can be passed using `...`.
#'   The returned result and its order depend on the function `FUN` and how
#'   `object` is split (hence on `f`, if provided). Parallel processing is
#'   supported and can be configured with parameter `BPPARAM`, is however only
#'   suggested for computational intense `FUN`.
#'   As an alternative to the (eventual parallel) processing of the full
#'   `Spectra`, `spectrapply()` supports also a chunk-wise processing. For this,
#'   parameter `chunkSize` needs to be specified. `object` is then split into
#'   chunks of size `chunkSize` which are then (stepwise) processed by `FUN`.
#'   This guarantees a lower memory demand (especially for on-disk backends)
#'   since only the data for one chunk needs to be loaded into memory in each
#'   iteration. Note that by specifying `chunkSize`, parameters `f` and
#'   `BPPARAM` will be ignored.
#'   See also `chunkapply()` above or examples below for details on chunk-wise
#'   processing.
#'
#'
#' @section The processing queue:
#'
#' Operations that modify mass peak data, i.e. the m/z and intensity values of
#' a `Spectra` are generally not applied immediately to the data but are
#' *cached* within the object's *processing queue*. These operations are then
#' applied to the data only upon request, for example when m/z and/or
#' intensity values are extracted. This lazy execution guarantees that the
#' same functionality can be applied to any `Spectra` object, regardless of
#' the type of backend that is used. Thus, data manipulation operations can
#' also be applied to data that is *read only*. As a side effect, this enables
#' also to *undo* operations using the `reset()` function.
#'
#' Functions related to the processing queue are:
#'
#' - `applyProcessing()`: for `Spectra` objects that use a **writeable** backend
#'   only: apply all steps from the lazy processing queue to the peak data and
#'   write it back to the data storage. Parameter `f` allows to specify how
#'   `object` should be split for parallel processing. This should either be
#'   equal to the `dataStorage`, or `f = rep(1, length(object))` to disable
#'   parallel processing alltogether. Other partitionings might result in
#'   errors (especially if a `MsBackendHdf5Peaks` backend is used).
#'
#' - `processingLog()`: returns a `character` vector with the processing log
#'   messages.
#'
#' - `reset()`: restores the data to its original state (as much as possible):
#'   removes any processing steps from the lazy processing queue and calls
#'   `reset()` on the backend which, depending on the backend, can also undo
#'   e.g. data filtering operations. Note that a `reset*(` call after
#'   `applyProcessing()` will not have any effect. See examples below for more
#'   information.
#'
#' @param binSize For `bin()`: `numeric(1)` defining the size for the m/z bins.
#'     Defaults to `binSize = 1`.
#'
#' @param BPPARAM Parallel setup configuration. See [bpparam()] for more
#'     information. This is passed directly to the [backendInitialize()] method
#'     of the [MsBackend-class]. See also [processingChunkSize()] for
#'     additional information on parallel processing.
#'
#' @param breaks For `bin()`: `numeric` defining the m/z breakpoints between
#'     bins.
#'
#' @param by For `scalePeaks()`: function to calculate a single `numeric` from
#'     intensity values of a spectrum by which all intensities (of
#'     that spectrum) should be divided by. The default `by = sum` will
#'     divide intensities of each spectrum by the sum of intensities of that
#'     spectrum.
#'
#' @param chunkSize For `spectrapply()`: size of the chunks into which the
#'     `Spectra` should be split. This parameter overrides parameters
#'     `f` and `BPPARAM`.
#'
#' @param descending For `pickPeaks()`: `logical`, if `TRUE` just values
#'     betwee the nearest valleys around the peak centroids are used.
#
#' @param f For `spectrapply()` and `applyProcessing()`: `factor` defining
#'     how `object` should be splitted for eventual parallel processing.
#'     Defaults to `factor()` for `spectrapply()` hence the object is not
#'     splitted while it defaults to `f = processingChunkSize(object)` for
#'     `applyProcessing()` splitting thus the object by default into chunks
#'     depending on [processingChunkSize()].
#'
#' @param FUN For `addProcessing()`: function to be applied to the peak matrix
#'     of each spectrum in `object`.
#'     For `bin()`: function to aggregate intensity values of peaks falling
#'     into the same bin. Defaults to `FUN = sum` thus summing up intensities.
#'     For `spectrapply()` and `chunkapply()`: function to be applied to
#'     each individual or each chunk of `Spectra`.
#'
#' @param halfWindowSize For `pickPeaks()`: `integer(1)`, used in the
#'     identification of the mass peaks: a local maximum has to be the
#'     maximum in the window from `(i - halfWindowSize):(i + halfWindowSize)`.
#'     For `smooth()`: `integer(1)`, used in the smoothing algorithm, the
#'     window reaches from `(i - halfWindowSize):(i + halfWindowSize)`.
#'
#' @param k For `pickPeaks()`: `integer(1)`, number of values left and right of
#'     the peak that should be considered in the weighted mean calculation.
#'
#' @param method For `pickPeaks()`: `character(1)`, the noise estimators that
#'     should be used, currently the the *M*edian *A*bsolute *D*eviation
#'     (`method = "MAD"`) and Friedman's Super Smoother
#'     (`method = "SuperSmoother"`) are supported.
#'     For `smooth()`: `character(1)`, the smoothing function that should be
#'     used, currently, the Moving-Average- (`method = "MovingAverage"`),
#'     Weighted-Moving-Average- (`method = "WeightedMovingAverage")`,
#'     Savitzky-Golay-Smoothing (`method = "SavitzkyGolay"`) are supported.
#'
#' @param msLevel. `integer` defining the MS level(s) of the spectra to which
#'     the function should be applied (defaults to all MS levels of `object`.
#'
#' @param mz For `containsMz()`: `numeric` with the m/z value(s) of the mass
#'     peaks to check.
#'
#' @param neutralLoss for `containsNeutralLoss()`: `numeric(1)` defining the
#'     value which should be subtracted from the spectrum's precursor m/z.
#'
#' @param normalized for `entropy()`: `logical(1)` whether the normalized
#'     entropy should be calculated (default). See also [nentropy()] for
#'     details.
#'
#' @param object A `Spectra` object.
#'
#' @param ppm For `containsMz()` and `neutralLoss()`: `numeric(1)` defining a
#'     relative, m/z-dependent, maximal accepted difference between m/z values
#'     for peaks to be matched.
#'
#' @param snr For `pickPeaks()`: `double(1)` defining the
#'     *S*ignal-to-*N*oise-*R*atio. The intensity of a local maximum has to be
#'     higher than `snr * noise` to be considered as peak.
#'
#' @param spectraVariables For `addProcessing()`: `character` with additional
#'     spectra variables that should be passed along to the function defined
#'     with `FUN`. See function description for details.
#'
#' @param threshold For `pickPeaks()`: a `numeric(1)` defining the proportion
#'     of the maximal peak intensity. Only values above the threshold are
#'     used for the weighted mean calculation.
#'     For `replaceIntensitiesBelow()`: a `numeric(1)` defining the threshold
#'     or a `function` to calculate the threshold for each spectrum on its
#'     intensity values. Defaults to `threshold = min`.
#'
#' @param tolerance For `containsMz()` and `neutralLoss()`:
#'     `numeric(1)` allowing to define a constant maximal accepted difference
#'     between m/z values for peaks to be matched.
#'
#' @param value For `replaceIntensitiesBelow()`: `numeric(1)` defining the
#'     value with which intensities should be replaced with.
#'
#' @param which For `containsMz()`: either `"any"` or `"all"` defining whether
#'     any (the default) or all provided `mz` have to be present in the
#'     spectrum.
#'
#' @param x A `Spectra`.
#'
#' @param zero.rm For `bin()`: `logical(1)` indicating whether to remove bins
#'     with zero intensity. Defaults to `TRUE`, meaning the function will
#'     discard bins created with an intensity of 0 to enhance memory
#'     efficiency.
#'
#' @param ... Additional arguments passed to internal and downstream functions.
#'
#' @return
#'
#' See the documentation of the individual functions for a description of the
#' return value.
#'
#' @md
#'
#' @seealso
#'
#' - [compareSpectra()] for calculation of spectra similarity scores.
#'
#' - [processingChunkSize()] for information on parallel and chunk-wise data
#'   processing.
#'
#' - [Spectra] for a general description of the `Spectra` object.
#'
#' @author Sebastian Gibb, Johannes Rainer, Laurent Gatto, Philippine Louail, Nir Shahaf, Mar Garcia-Aloy
#'
#' @examples
#'
#' ## Load a `Spectra` object with LC-MS/MS data.
#' fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML",
#'     package = "msdata")
#' sps_dda <- Spectra(fl)
#' sps_dda
#'
#'
#' ##  --------  FUNCTIONS RETURNING A SPECTRA  --------
#'
#' ## Replace peak intensities below 40 with a value of 1
#' sps_mod <- replaceIntensitiesBelow(sps_dda, threshold = 20, value = 1)
#' sps_mod
#'
#' ## Get the intensities of the first spectrum before and after the
#' ## operation
#' intensity(sps_dda[1])
#' intensity(sps_mod[1])
#'
#' ## Remove all peaks with an intensity below 5.
#' sps_mod <- filterIntensity(sps_dda, intensity = c(5, Inf))
#'
#' intensity(sps_mod)
#'
#' ## In addition it is possible to pass a function to `filterIntensity()`: in
#' ## the example below we want to keep only peaks that have an intensity which
#' ## is larger than one third of the maximal peak intensity in that spectrum.
#' keep_peaks <- function(x, prop = 3) {
#'     x > max(x, na.rm = TRUE) / prop
#' }
#' sps_mod <- filterIntensity(sps_dda, intensity = keep_peaks)
#' intensity(sps_mod)
#'
#' ## We can also change the proportion by simply passing the `prop` parameter
#' ## to the function. To keep only peaks that have an intensity which is
#' ## larger than half of the maximum intensity:
#' sps_mod <- filterIntensity(sps_dda, intensity = keep_peaks, prop = 2)
#' intensity(sps_mod)
#'
#' ## With the `scalePeaks()` function we can alternatively scale the
#' ## intensities of mass peaks per spectrum to relative intensities. This
#' ## is specifically useful for fragment (MS2) spectra. We below thus
#' ## scale the intensities per spectrum by the total sum of intensities
#' ## (such that the sum of all intensities per spectrum is 1).
#' ## Below we scale the intensities of all MS2 spectra in our data set.
#' sps_mod <- scalePeaks(sps_dda, msLevel = 2L)
#'
#' ## MS1 spectra were not affected
#' sps_mod |>
#'     filterMsLevel(1L) |>
#'     intensity()
#'
#' ## Intensities of MS2 spectra were scaled
#' sps_mod |>
#'     filterMsLevel(2L) |>
#'     intensity()
#'
#' ## Since data manipulation operations are by default not directly applied to
#' ## the data but only cached in the internal processing queue, it is also
#' ## possible to remove these data manipulations with the `reset()` function:
#' tmp <- reset(sps_mod)
#' tmp
#' lengths(sps_dda) |> head()
#' lengths(sps_mod) |> head()
#' lengths(tmp) |> head()
#'
#' ## Data manipulation operations cached in the processing queue can also be
#' ## applied to the mass peaks data with the `applyProcessing()` function, if
#' ## the `Spectra` uses a backend that supports that (i.e. allows replacing
#' ## the mass peaks data). Below we first change the backend to a
#' ## `MsBackendMemory()` and then use the `applyProcessing()` to modify the
#' ## mass peaks data
#' sps_dda <- setBackend(sps_dda, MsBackendMemory())
#' sps_mod <- filterIntensity(sps_dda, intensity = c(5, Inf))
#' sps_mod <- applyProcessing(sps_mod)
#' sps_mod
#'
#' ## While we can't *undo* this filtering operation now using the `reset()`
#' ## function, accessing the data would now be faster, because the operation
#' ## does no longer to be applied to the original data before returning to the
#' ## user.
#'
#'
#' ##  --------  FUNCTIONS RETURNING THE RESULT  --------
#'
#' ## With the `spectrapply()` function it is possible to apply an
#' ## arbitrary function to each spectrum in a Spectra.
#' ## In the example below we calculate the mean intensity for each spectrum
#' ## in a subset of the sciex_im data. Note that we can access all variables
#' ## of each individual spectrum either with the `$` operator or the
#' ## corresponding method.
#' res <- spectrapply(sps_dda[1:20], FUN = function(x) mean(x$intensity[[1]]))
#' head(res)
#'
#' ## As an alternative, applying a function `FUN` to a `Spectra` can be
#' ## performed *chunk-wise*. The advantage of this is, that only the data for
#' ## one chunk at a time needs to be loaded into memory reducing the memory
#' ## demand. This type of processing can be performed by specifying the size
#' ## of the chunks (i.e. number of spectra per chunk) with the `chunkSize`
#' ## parameter
#' spectrapply(sps_dda[1:20], lengths, chunkSize = 5L)
#'
#' ## Precursor intensity estimation. Some manufacturers don't report the
#' ## precursor intensity for MS2 spectra:
#' sps_dda |>
#'     filterMsLevel(2L) |>
#'     precursorIntensity()
#'
#' ## This intensity can however be estimated from the previously measured
#' ## MS1 scan with the `estimatePrecursorIntensity()` function:
#' pi <- estimatePrecursorIntensity(sps_dda)
#'
#' ## This function returned the result as a `numeric` vector with one
#' ## value per spectrum:
#' pi
#'
#' ## We can replace the precursor intensity values of the originating
#' ## object:
#' sps_dda$precursorIntensity <- pi
#' sps_dda |>
#'     filterMsLevel(2L) |>
#'     precursorIntensity()
#'
NULL

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
#' @rdname addProcessing
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

#' @rdname addProcessing
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

#' @rdname addProcessing
#'
#' @exportMethod containsMz
setMethod("containsMz", "Spectra", function(object, mz = numeric(),
                                            tolerance = 0,
                                            ppm = 20, which = c("any", "all"),
                                            BPPARAM = bpparam()) {
    if (length(object)) {
        cond_fun <- match.fun(match.arg(which))
        if (all(is.na(mz)))
            return(rep(NA, length(object)))
        mz <- unique(sort(mz))
        BPPARAM <- backendBpparam(object@backend, BPPARAM)
        unlist(.peaksapply(
            object, FUN = .peaks_contain_mz, mz = mz, tolerance = tolerance,
            ppm = ppm, condFun = cond_fun, BPPARAM = BPPARAM),
            use.names = FALSE
            )
    } else logical()
})

#' @rdname addProcessing
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

#' @rdname addProcessing
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
#' @rdname addProcessing
setMethod("entropy", "ANY", function(object, ...) {
    MsCoreUtils::entropy(object)
})

#' @rdname addProcessing
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

#' @rdname addProcessing
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

#' @rdname addProcessing
setMethod("reset", "Spectra", function(object, ...) {
    object@backend <- reset(object@backend)
    object@processingQueue <- list()
    if (!.hasSlot(object, "processingQueueVariables"))
        object <- updateObject(object, check = FALSE)
    object@processingQueueVariables <- character()
    object@processing <- .logging(object@processing, "Reset object.")
    object
})

#' @rdname addProcessing
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

#' @rdname addProcessing
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

#' @title Estimate Precursor Intensities
#'
#' @aliases estimatePrecursorIntensity
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


################################################################################
##
## Spectra similarity calculations
##
################################################################################

#' @title Spectra similarity calculations
#'
#' @name compareSpectra
#'
#' @aliases compareSpectra
#'
#' @description
#'
#' `compareSpectra()` compares each spectrum in `x` with each spectrum in `y`
#' using the function provided with `FUN` (defaults to [ndotproduct()]). If
#' `y` is missing, each spectrum in `x` is compared with each other spectrum
#' in `x`.
#' The matching/mapping of peaks between the compared spectra is done with the
#' `MAPFUN` function. The default [joinPeaks()] matches peaks of both spectra
#' and allows to keep all peaks from the first spectrum (`type = "left"`),
#' from the second (`type = "right"`), from both (`type = "outer"`) and to
#' keep only matching peaks (`type = "inner"`); see [joinPeaks()] for more
#' information and examples). The `MAPFUN` function should have parameters
#' `x`, `y`, `xPrecursorMz` and `yPrecursorMz` as these values are passed to
#' the function.
#'
#' In addition to `joinPeaks()` also [joinPeaksGnps()] is supported for
#' GNPS-like similarity score calculations. Note that `joinPeaksGnps()` should
#' only be used in combination with `FUN = MsCoreUtils::gnps`
#' (see [joinPeaksGnps()] for more information and details). Use
#' `MAPFUN = joinPeaksNone` to disable internal peak matching/mapping if a
#' similarity scoring function is used that performs the matching internally.
#'
#' `FUN` is supposed to be a function to compare intensities of (matched)
#' peaks of the two spectra that are compared. The function needs to take two
#' matrices with columns `"mz"` and `"intensity"` as input and is supposed
#' to return a single numeric as result. In addition to the two peak matrices
#' the spectra's precursor m/z values are passed to the function as parameters
#' `xPrecursorMz` (precursor m/z of the `x` peak matrix) and `yPrecursorMz`
#' (precursor m/z of the `y` peak matrix). Additional parameters to functions
#' `FUN` and `MAPFUN` can be passed with `...`. Parameters `ppm` and
#' `tolerance` are passed to both `MAPFUN` and `FUN`.
#' The function returns a `matrix` with the results of `FUN` for each
#' comparison, number of rows equal to `length(x)` and number of columns
#' equal `length(y)` (i.e. element in row 2 and column 3 is the result from
#' the comparison of `x[2]` with `y[3]`). If `SIMPLIFY = TRUE` the `matrix`
#' is *simplified* to a `numeric` if length of `x` or `y` is one. See also
#' the vignette for additional examples, such as using spectral entropy
#' similarity in the scoring.
#'
#' @param FUN function to compare intensities of peaks between two spectra.
#'     Defaults to [ndotproduct()].
#'
#' @param MAPFUN For `compareSpectra()`: function to map/match peaks between
#'     the two compared spectra. See [joinPeaks()] for more information and
#'     possible functions. Defaults to [joinPeaks()].
#'
#' @param ppm `numeric(1)` defining a relative, m/z-dependent, maximal
#'     accepted difference between m/z values for peaks to be matched. This
#'     parameter is directly passed to `MAPFUN`.
#'
#' @param tolerance `numeric(1)` allowing to define a constant maximal
#'     accepted difference between m/z values for peaks to be matched. This
#'     parameter is directly passed to `MAPFUN`.
#'
#' @param x A `Spectra` object.
#'
#' @param y A `Spectra` object.
#'
#' @param SIMPLIFY `logical(1)` defining  whether the result matrix should be
#'     *simplified* to a `numeric` if possible (i.e. if either `x` or `y` is
#'     of length 1).
#'
#' @param ... Additional arguments passed to the internal functions.
#'
#' @importFrom MsCoreUtils ndotproduct
#'
#' @importMethodsFrom ProtGenerics compareSpectra
#'
#' @exportMethod compareSpectra
#'
#' @author Sebastian Gibb, Johannes Rainer, Laurent Gatto
#'
#' @examples
#'
#' ## Load a `Spectra` object with LC-MS/MS data.
#' fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML",
#'     package = "msdata")
#' sps_dda <- Spectra(fl)
#' sps_dda
#'
#' ## Restrict to MS2 (fragment) spectra:
#' sps_ms2 <- filterMsLevel(sps_dda, msLevel = 2L)
#'
#' ## Compare spectra: comparing spectra 2 and 3 against spectra 10:20 using
#' ## the normalized dotproduct method.
#' res <- compareSpectra(sps_ms2[2:3], sps_ms2[10:20])
#' ## first row contains comparisons of spectrum 2 with spectra 10 to 20 and
#' ## the second row comparisons of spectrum 3 with spectra 10 to 20
#' res
#'
#' ## We next calculate the pairwise similarity for the first 10 spectra
#' compareSpectra(sps_ms2[1:10])
#'
#' ## Use compareSpectra to determine the number of common (matching) peaks
#' ## with a ppm of 10:
#' ## type = "inner" uses a *inner join* to match peaks, i.e. keeps only
#' ## peaks that can be mapped betwen both spectra. The provided FUN returns
#' ## simply the number of matching peaks.
#' compareSpectra(sps_ms2[2:3], sps_ms2[10:20], ppm = 10, type = "inner",
#'     FUN = function(x, y, ...) nrow(x))
#'
#' ## We repeat this calculation between all pairwise combinations
#' ## of the first 20 spectra
#' compareSpectra(sps_ms2[1:20], ppm = 10, type = "inner",
#'     FUN = function(x, y, ...) nrow(x))
NULL

#' @rdname compareSpectra
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
#' @rdname compareSpectra
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


################################################################################
##
## methods with documentation in Spectra-functions.R
##
################################################################################

#' @rdname processingChunkSize
setMethod("backendBpparam", "Spectra", function(object, BPPARAM = bpparam()) {
    backendBpparam(object@backend, BPPARAM)
})
