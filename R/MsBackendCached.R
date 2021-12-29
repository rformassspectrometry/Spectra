setClassUnion("characterOrInteger", c("character", "integer"))

#' @title Base MsBackend class providing data caching mechanism
#'
#' @aliases MsBackendCached-class MsBackendCached
#'
#' @description
#'
#' The `MsBackendCached` class is a rudimentary implementation of the
#' [MsBackend] providing a simple mechanism to cache spectra data locally. This
#' class is thought to be used as a base class for other `MsBackend`
#' implementations to reuse its caching mechanism and avoid having to
#' re-implement commonly used methods. This class is thus not thought to be
#' used directly by a user.
#'
#' The provided caching mechanism by this class allows `MsBackend` instances
#' to add or replace spectra variables even if the backend used by them does
#' not allow to alter values (e.g. if a SQL database is used as a backend). Any
#' replacement operation with `$<-` will add the specified values to a local
#' `data.frame` within the `MsBackendCached` class that allows to *cache* these
#' values (increasing obviously the memory demand of the object).
#'
#' Any data accessor functions of the extending `MsBackend` class (such as `$`
#' or `msLevel` or `spectraData`) should first use `callNextMethod` to call the
#' respective accessor of `MsBackendCached` that will evaluate if the
#' requested spectra variable(s) are in the local cache and return these. If
#' the requested spectra variables are not in the local cache, are also not in
#' listed in the `@spectraVariables` slot (which defines all spectra variables
#' that can be requested from the extending `MsBackend` class) but are *core
#' spectra variables* then missing values of the correct data type are returned.
#'
#' @section Implementation notes:
#'
#' Classes extending the `MsBackendCached` need to
#'
#' - call the `backendInitialize` method of this class in their own
#'   `backendInitialize` method and to set at least the number of spectra with
#'   the `nspectra` parameter and the `spectraVariables` that are available to
#'   the (extending) backend class.
#'
#' - implement the `spectraData` method that calls also the `spectraData`
#'   method from `MsBackendCached` to also retrieve cached values (e.g. using
#'   `res <- callNextMethod()` at the beginning of the `spectraData` function).
#'   The `spectraData,MsBackendCached` method will return `NULL` if the
#'   selected spectra variables were not cached and are not *core spectra
#'   variables* not being provided by the extending backend. Thus, the extending
#'   backend can then proceed to retrieve the respective values from its own
#'   backend/data storage.
#'
#' - implement eventually the `[` method that calls in addition the `[` from the
#'   `MsBackendCached`.
#'
#' All other methods accessing or setting spectra variables don't need to be
#' implemented by the extending backend class (the default implementations of
#' the `MsBackendCached` will then be used instead; these ensure that cached
#' values are returned first).
#' Spectra variables can be modified or added using the `$<-` method of the
#' `MsBackendCached`. Replacing or adding multiple variables using the
#' `spectraData<-` is not supported by `MsBackendCached`. The extending backend
#' might however implement such a method that internally uses `$<-` to
#' add/replace single variables.
#'
#' The `MsBackendCached` has the following slots"
#'
#' - `nspectra`: `integer(1)` defining the number of spectra of the backend.
#'   This variable needs to be set and must match the number of rows of
#'   `localData` and the actual number of spectra in the (extending) backend.
#'
#' - `localData`: `data.frame` with the cached local data. Any replacement
#'   operation with `$<-` will set/add a column with the respective values.
#'
#' - `spectraVariables`: `character` defining the spectra variables that are
#'   provided by the extending `MsBackend` class (e.g. all spectra variables
#'   that can be retrieved from the data base or original data files).
#'
#' @section Available methods:
#'
#' - `acquisitionNum`: returns the acquisition number of each
#'   spectrum. Returns an `integer` of length equal to the number of
#'   spectra (with `NA_integer_` if not available).
#'
#' - `backendInitialize`: *initializes* the backend. The method takes parameters
#'   `data` (`data.frame` with cached data), `nspectra` (`integer` defining the
#'   number of spectra) and `spectraVariables` (`character` with the spectra
#'   variables that are provided by the extending backend.
#'
#' - `centroided`, `centroided<-`: gets or sets the centroiding
#'   information of the spectra. `centroided` returns a `logical`
#'   vector of length equal to the number of spectra with `TRUE` if a
#'   spectrum is centroided, `FALSE` if it is in profile mode and `NA`
#'   if it is undefined. See also `isCentroided` for estimating from
#'   the spectrum data whether the spectrum is centroided.  `value`
#'   for `centroided<-` is either a single `logical` or a `logical` of
#'   length equal to the number of spectra in `object`.
#'
#' - `collisionEnergy`, `collisionEnergy<-`: gets or sets the
#'   collision energy for all spectra in `object`. `collisionEnergy`
#'   returns a `numeric` with length equal to the number of spectra
#'   (`NA_real_` if not present/defined), `collisionEnergy<-` takes a
#'   `numeric` of length equal to the number of spectra in `object`.
#'
#' - `dataOrigin`: gets a `character` of length equal to the number of spectra
#'   in `object` with the *data origin* of each spectrum. This could e.g. be
#'   the mzML file from which the data was read.
#'
#' - `intensity`: gets the intensity values from the spectra. Returns
#'   a [NumericList()] of `numeric` vectors (intensity values for each
#'   spectrum). The length of the `list` is equal to the number of
#'   `spectra` in `object`.
#'
#' - `ionCount`: returns a `numeric` with the sum of intensities for
#'   each spectrum. If the spectrum is empty (see `isEmpty`),
#'   `NA_real_` is returned.
#'
#' - `isEmpty`: checks whether a spectrum in `object` is empty
#'   (i.e. does not contain any peaks). Returns a `logical` vector of
#'   length equal number of spectra.
#'
#' - `isolationWindowLowerMz`, `isolationWindowLowerMz<-`: gets or sets the
#'   lower m/z boundary of the isolation window.
#'
#' - `isolationWindowTargetMz`, `isolationWindowTargetMz<-`: gets or sets the
#'   target m/z of the isolation window.
#'
#' - `isolationWindowUpperMz`, `isolationWindowUpperMz<-`: gets or sets the
#'   upper m/z boundary of the isolation window.
#'
#' - `length`: returns the number of spectra (i.e. the `@nspectra`).
#'
#' - `lengths`: gets the number of peaks (m/z-intensity values) per
#'   spectrum.  Returns an `integer` vector (length equal to the
#'   number of spectra). For empty spectra, `0` is returned.
#'
#' - `msLevel`: gets the spectra's MS level. Returns an `integer`
#'   vector (of length equal to the number of spectra) with the MS
#'   level for each spectrum (or `NA_integer_` if not available).
#'
#' - `mz`: gets the mass-to-charge ratios (m/z) from the
#'   spectra. Returns a [NumericList()] or length equal to the number of
#'   spectra, each element a `numeric` vector with the m/z values of
#'   one spectrum.
#'
#' - `polarity`, `polarity<-`: gets or sets the polarity for each
#'   spectrum.  `polarity` returns an `integer` vector (length equal
#'   to the number of spectra), with `0` and `1` representing negative
#'   and positive polarities, respectively. `polarity<-` expects an
#'   integer vector of length 1 or equal to the number of spectra.
#'
#' - `precursorCharge`, `precursorIntensity`, `precursorMz`,
#'   `precScanNum`, `precAcquisitionNum`: get the charge (`integer`),
#'   intensity (`numeric`), m/z (`numeric`), scan index (`integer`)
#'   and acquisition number (`interger`) of the precursor for MS level
#'   2 and above spectra from the object. Returns a vector of length equal to
#'   the number of spectra in `object`. `NA` are reported for MS1
#'   spectra of if no precursor information is available.
#'
#' - `rtime`, `rtime<-`: gets or sets the retention times for each
#'   spectrum (in seconds). `rtime` returns a `numeric` vector (length equal to
#'   the number of spectra) with the retention time for each spectrum.
#'   `rtime<-` expects a numeric vector with length equal to the
#'   number of spectra.
#'
#' - `scanIndex`: returns an `integer` vector with the *scan index*
#'   for each spectrum. This represents the relative index of the
#'   spectrum within each file. Note that this can be different to the
#'   `acquisitionNum` of the spectrum which is the index of the
#'   spectrum as reported in the mzML file.
#'
#' - `selectSpectraVariables`: subset the object to specified spectra variables.
#'   This will eventually remove spectra variables listed in `@spectraVariables`
#'   and will also drop columns from the local cache if not among
#'   `spectraVariables`.
#'
#' - `smoothed`,`smoothed<-`: gets or sets whether a spectrum is
#'   *smoothed*. `smoothed` returns a `logical` vector of length equal
#'   to the number of spectra. `smoothed<-` takes a `logical` vector
#'   of length 1 or equal to the number of spectra in `object`.
#'
#' - `spectraVariables`: returns the available spectra variables, i.e. the
#'   unique set of *core spectra variables*, cached spectra variables and
#'   spectra variables defined in the `@spectraVariables` slot (i.e. spectra
#'   variables thought to be provided by the extending `MsBackend` instance).
#'
#' - `spectraData`: returns a `DataFrame` with cached spectra variablers or
#'   initialized *core spectra variables*. Parameter `spectraVariables` allows
#'   to specify the variables to retrieve. The function returns `NULL` if the
#'   requested variables are not cached and are not provided by the extending
#'   backend. Note that this method **only** returns cached spectra variables
#'   or core spectra variables **not** provided by the extending backend. It is
#'   the responsibility of the extending backend to add/provide these.
#'
#' - `[`: subsets the cached data. Parameter `i` needs to be an `integer`
#'   vector.
#'
#' - `$`, `$<-`: access or set/add a single spectrum variable (column) in the
#'   backend.
#'
#'
#' @param columns For `spectraData`: `character` with the names of the spectra
#'     variables to retrieve.
#'
#' @param data For `backendInitialize`: (optional) `data.frame` with cached
#'     values. The number of rows (and their order) has to match the number of
#'     spectra.
#'
#' @param drop For `[`: not considered.
#'
#' @param i For `[`: `integer` with the indices to subset the object.
#'
#' @param j For `[`: ignored.
#'
#' @param name For `$<-`: the name of the spectra variable to set.
#'
#' @param nspectra For `backendInitialize`: `integer` with the number of
#'     spectra.
#'
#' @param object A `MsBackendCached` object.
#'
#' @param spectraVariables For `backendInitialize`: `character` with the names
#'     of the spectra variables that are provided by the extending backend.
#'     For `selectSpectraVariables`: `character` specifying the spectra
#'     variables to keep.
#'
#' @param use.names For `lengths`: whether spectrum names should be used.
#'
#' @param value replacement value for `<-` methods. See individual
#'     method description or expected data type.
#'
#' @param x A `MsBackendCached` object.
#'
#' @param ... ignored
#'
#' @author Johannes Rainer
#'
#' @name MsBackendCached
#'
#' @return See documentation of respective function.
#'
#' @seealso [MsBackend] for the documentation of MS backends.
#'
#' @exportClass MsBackendCached
NULL

setClass(
    "MsBackendCached",
    contains = "MsBackend",
    slots = c(
        spectraVariables = "character",
        localData = "data.frame",
        nspectra = "integer"),
    prototype = prototype(
        spectraVariables = character(),
        localData = data.frame(),
        nspectra = 0L,
        readonly = FALSE, version = "0.1"))

.valid_local_data <- function(x, y) {
    if (nrow(x) && nrow(x) != y)
        "Number of rows in local data and number of spectra don't match"
    else NULL
}

#' @rdname MsBackendCached
#'
#' @export
MsBackendCached <- function() {
    new("MsBackendCached")
}

setValidity("MsBackendCached", function(object) {
    msg <- c(.valid_local_data(object@localData, object@nspectra))
    ## For local data that are core spectra variables, check for valid data type
    msg <- c(msg, .valid_column_datatype(object@localData))
    if (is.null(msg)) TRUE
    else msg
})

#' @rdname MsBackendCached
setMethod("backendInitialize", "MsBackendCached",
          function(object, data = data.frame(), nspectra = 0L,
                   spectraVariables = character(), ...) {
              if (nrow(data) == 0L)
                  data <- data.frame(row.names = seq_len(nspectra))
              else nspectra <- nrow(data)
              object@nspectra <- as.integer(nspectra)
              object@localData <- data
              object@spectraVariables <- spectraVariables
              validObject(object)
              object
})

#' @rdname MsBackendCached
setMethod("dataStorage", "MsBackendCached", function(object) {
    rep("<cache>", length(object))
})

#' @rdname MsBackendCached
setMethod("length", "MsBackendCached", function(x) {
    x@nspectra
})

#' @rdname MsBackendCached
setMethod("spectraVariables", "MsBackendCached", function(object) {
    unique(c(names(.SPECTRA_DATA_COLUMNS), colnames(object@localData),
             object@spectraVariables))
})

#' @importFrom S4Vectors make_zero_col_DFrame
#'
#' @noRd
.spectra_data <- function(x, columns = spectraVariables(x)) {
    local_cols <- intersect(columns, colnames(x@localData))
    ## @spectraVariables are supposed to be added by the extending backend
    columns <- columns[!columns %in% x@spectraVariables]
    core_cols <- intersect(columns, names(.SPECTRA_DATA_COLUMNS))
    core_cols <- core_cols[!core_cols %in% c(local_cols, x@spectraVariables)]
    res <- NULL
    if (length(local_cols))
        res <- as(x@localData[, local_cols, drop = FALSE], "DataFrame")
    if (length(core_cols)) {
        mzvals <- NULL
        intvals <- NULL
        lst <- NumericList(numeric(), compress = FALSE)
        if (any(core_cols == "mz")) {
            core_cols <- core_cols[core_cols != "mz"]
            mzvals <- lst[rep(1, times = length(x))]
        }
        if (any(core_cols == "intensity")) {
            core_cols <- core_cols[core_cols != "intensity"]
            intvals <- lst[rep(1, times = length(x))]
        }
        if (length(core_cols))
            tmp <- DataFrame(lapply(.SPECTRA_DATA_COLUMNS[core_cols],
                                    function(z, n) rep(as(NA, z), n), length(x)))
        else tmp <- make_zero_col_DFrame(x@nspectra)
        tmp$mz <- mzvals
        tmp$intensity <- intvals
        if (length(res))
            res <- cbind(res, tmp)
        else res <- tmp
        if (any(core_cols == "dataStorage"))
            res$dataStorage <- dataStorage(x)
    }
    if (!all(columns %in% colnames(res)))
        stop("Column(s) ", paste0(columns[!columns %in% names(res)],
                                  collapse = ", "), " not available.",
             call. = FALSE)
    if (length(res))
        as(res, "DataFrame")
    else NULL
}

#' @rdname MsBackendCached
setMethod(
    "spectraData", "MsBackendCached",
    function(object, columns = spectraVariables(object)) {
        .spectra_data(object, columns = columns)
    })

#' @rdname MsBackendCached
setReplaceMethod("spectraData", "MsBackendCached", function(object, value) {
    stop(class(object)[1], " does not support replacing the full spectra data.")
})

#' @rdname MsBackendCached
setMethod("[", "MsBackendCached", function(x, i, j, ..., drop = FALSE) {
    if (missing(i))
        return(x)
    slot(x, "localData", check = FALSE) <- x@localData[i, , drop = FALSE]
    x@nspectra <- nrow(x@localData)
    x
})

#' @rdname MsBackendCached
setMethod("$", "MsBackendCached", function(x, name) {
    if (!any(spectraVariables(x) == name))
        stop("Spectra variable '", name, "' not available.")
    spectraData(x, name)[, 1]
})

#' @rdname MsBackendCached
setReplaceMethod("$", "MsBackendCached", function(x, name, value) {
    if (name %in% c("mz", "intensity"))
        stop("Replacing m/z and intensity values is not supported.",
             call. = FALSE)
    if (length(value) == 1)
        value <- rep(value, length(x))
    if (length(value) != length(x))
        stop("value has to be either of length 1 or length equal to the ",
             "number of spectra")
    if (length(x@localData)) {
        cn <- colnames(x@localData) == name
        if (any(cn))
            x@localData[, cn] <- value
        else {
            cn <- colnames(x@localData)
            x@localData <- cbind(x@localData, value)
            colnames(x@localData) <- c(cn, name)
        }
    } else {
        x@localData <- data.frame(value)
        colnames(x@localData) <- name
    }
    validObject(x)
    x
})

#' @rdname MsBackendCached
setMethod(
    "selectSpectraVariables", "MsBackendCached",
    function(object, spectraVariables = spectraVariables(object)) {
        if (any(!spectraVariables %in% spectraVariables(object)))
            stop("spectra variable(s) ",
                 paste(spectraVariables[!spectraVariables %in%
                                        spectraVariables(object)],
                       collapse = ", "), " not available")
        object@spectraVariables <- intersect(object@spectraVariables,
                                             spectraVariables)
        object@localData <- object@localData[, colnames(object@localData) %in%
                                               spectraVariables, drop = FALSE]
        validObject(object)
        object
    })

#' @rdname MsBackendCached
setMethod("show", "MsBackendCached", function(object) {
    n <- object@nspectra
    cat(class(object), "with", n, "spectra\n")
    if (n) {
        idx <- unique(c(1L:min(6L, n), max(1L, n-5L):n))
        spd <- spectraData(object[idx, ],
                           c("msLevel", "precursorMz", "polarity"))
        if (!length(rownames(spd)))
            rownames(spd) <- idx
        txt <- capture.output(print(spd))
        cat(txt[-1], sep = "\n")
        sp_cols <- spectraVariables(object)
        cat(" ...", length(sp_cols) - 3, "more variables/columns.\n", "Use ",
            "'spectraVariables' to list all of them.\n")
    }
})

## ------------- DEFAULT IMPLEMENTATIONS OF ACCESSORS --------------------------

#' @rdname MsBackendCached
setMethod("acquisitionNum", "MsBackendCached", function(object) {
    spectraData(object, "acquisitionNum")[, 1]
})

#' @rdname MsBackendCached
setMethod("centroided", "MsBackendCached", function(object) {
    spectraData(object, "centroided")[, 1]
})

#' @rdname MsBackendCached
setReplaceMethod("centroided", "MsBackendCached", function(object, value) {
    object$centroided <- value
    validObject(object)
    object
})

#' @rdname MsBackendCached
setMethod("collisionEnergy", "MsBackendCached", function(object) {
    spectraData(object, "collisionEnergy")[, 1]
})

#' @rdname MsBackendCached
setReplaceMethod("collisionEnergy", "MsBackendCached", function(object, value) {
    object$collisionEnergy <- value
    validObject(object)
    object
})

#' @rdname MsBackendCached
setMethod("dataOrigin", "MsBackendCached", function(object) {
    spectraData(object, "dataOrigin")[, 1]
})

#' @rdname MsBackendCached
setReplaceMethod("dataOrigin", "MsBackendCached", function(object, value) {
    object$dataOrigin <- value
    validObject(object)
    object
})

#' @rdname MsBackendCached
setMethod("msLevel", "MsBackendCached", function(object) {
    spectraData(object, "msLevel")[, 1]
})

#' @rdname MsBackendCached
setMethod("intensity", "MsBackendCached", function(object) {
    spectraData(object, "intensity")[, 1]
})

#' @importFrom MsCoreUtils vapply1d
#'
#' @rdname MsBackendCached
setMethod("ionCount", "MsBackendCached", function(object) {
    vapply1d(intensity(object), sum, na.rm = TRUE)
})

#' @importMethodsFrom S4Vectors isEmpty
#'
#' @rdname MsBackendCached
setMethod("isEmpty", "MsBackendCached", function(x) {
    lengths(intensity(x)) == 0
})

#' @rdname MsBackendCached
setMethod("isolationWindowLowerMz", "MsBackendCached", function(object) {
    spectraData(object, "isolationWindowLowerMz")[, 1]
})

#' @rdname MsBackendCached
setReplaceMethod("isolationWindowLowerMz", "MsBackendCached",
                 function(object, value) {
                     object$isolationWindowLowerMz <- value
                     validObject(object)
                     object
                 })

#' @rdname MsBackendCached
setMethod("isolationWindowTargetMz", "MsBackendCached", function(object) {
    spectraData(object, "isolationWindowTargetMz")[, 1]
})

#' @rdname MsBackendCached
setReplaceMethod("isolationWindowTargetMz", "MsBackendCached",
                 function(object, value) {
                     object$isolationWindowTargetMz <- value
                     validObject(object)
                     object
                 })

#' @rdname MsBackendCached
setMethod("isolationWindowUpperMz", "MsBackendCached", function(object) {
    spectraData(object, "isolationWindowUpperMz")[, 1]
})

#' @rdname MsBackendCached
setReplaceMethod("isolationWindowUpperMz", "MsBackendCached",
                 function(object, value) {
                     object$isolationWindowUpperMz <- value
                     validObject(object)
                     object
                 })

#' @rdname MsBackendCached
setMethod("lengths", "MsBackendCached", function(x, use.names = FALSE) {
    lengths(mz(x))
})

#' @rdname MsBackendCached
setMethod("mz", "MsBackendCached", function(object) {
    spectraData(object, "mz")[, 1]
})

#' @rdname MsBackendCached
setMethod("polarity", "MsBackendCached", function(object) {
    spectraData(object, "polarity")[, 1]
})

#' @rdname MsBackendCached
setReplaceMethod("polarity", "MsBackendCached", function(object, value) {
    if (is.numeric(value)) value <- as.integer(value)
    object$polarity <- value
    validObject(object)
    object
})

#' @rdname MsBackendCached
setMethod("precursorCharge", "MsBackendCached", function(object) {
    spectraData(object, "precursorCharge")[, 1]
})

#' @rdname MsBackendCached
setMethod("precursorIntensity", "MsBackendCached", function(object) {
    spectraData(object, "precursorIntensity")[, 1]
})

#' @rdname MsBackendCached
setMethod("precursorMz", "MsBackendCached", function(object) {
    spectraData(object, "precursorMz")[, 1]
})

#' @rdname MsBackendCached
setMethod("rtime", "MsBackendCached", function(object) {
    spectraData(object, "rtime")[, 1]
})

#' @rdname MsBackendCached
setReplaceMethod("rtime", "MsBackendCached", function(object, value) {
    object$rtime <- value
    validObject(object)
    object
})

#' @rdname MsBackendCached
setMethod("scanIndex", "MsBackendCached", function(object) {
    spectraData(object, "scanIndex")[, 1]
})

#' @rdname MsBackendCached
setMethod("smoothed", "MsBackendCached", function(object) {
    spectraData(object, "smoothed")[, 1]
})

#' @rdname MsBackendCached
setReplaceMethod("smoothed", "MsBackendCached", function(object, value) {
    object$smoothed <- value
    validObject(object)
    object
})
