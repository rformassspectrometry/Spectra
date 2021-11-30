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
#' - use the `spectraData` function in their own data accessor methods to
#'   retrieve potentially cached values. The function will return `NULL` if the
#'   selected spectra variables were not cached and are not *core spectra
#'   variables* not being provided by the extending backend. Thus, the extending
#'   backend can then proceed to retrieve the respective values from its own
#'   backend/data storage.
#'
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
#' - `backendInitialize`: *initializes* the backend. The method takes parameters
#'   `data` (`data.frame` with cached data), `nspectra` (`integer` defining the
#'   number of spectra) and `spectraVariables` (`character` with the spectra
#'   variables that are provided by the extending backend.
#'
#' - `length`: returns the number of spectra (i.e. the `@nspectra`).
#'
#' - `spectraVariables`: returns the available spectra variables, i.e. the
#'   unique set of *core spectra variables*, cached spectra variables and
#'   spectra variables defined in the `@spectraVariables` slot (i.e. spectra
#'   variables thought to be provided by the extending `MsBackend` instance).
#'
#' - `selectSpectraVariables`: subset the object to specified spectra variables.
#'   This will eventually remove spectra variables listed in `@spectraVariables`
#'   and will also drop columns from the local cache if not among
#'   `spectraVariables`.
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
#' @md
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
                  data <- data.frame(matrix(ncol = 0, nrow = nspectra))
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
setReplaceMethod("spectraData", "MsBackendCached",function(object, value) {
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
        idx <- union(1:min(6, n), max(1, n-5):n)
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
