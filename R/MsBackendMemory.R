#' @include hidden_aliases.R
NULL

#' @title In-memory MS data backend
#'
#' @description
#'
#' @name MsBackendMemory
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @noRd
#'
#' @exportClass MsBackendMemory
NULL

setClass("MsBackendMemory",
         contains = "MsBackend",
         slots = c(spectraData = "DataFrame"),
         prototype = prototype(spectraData = DataFrame(),
                               files = character(),
                               modCount = integer(),
                               readonly = FALSE,
                               version = "0.1"))

setValidity("MsBackendMemory", function(object) {
    msg <- .valid_spectra_data_required_columns(object@spectraData)
    if (length(msg))
        return(msg)
    msg <- c(
        .valid_column_datatype(object@spectraData, .SPECTRA_DATA_COLUMNS),
        .valid_intensity_column(object@spectraData),
        .valid_mz_column(object@spectraData),
        .valid_ms_backend_files(object@files),
        .valid_ms_backend_files_from_file(object@files,
                                          object@spectraData$fromFile),
        .valid_ms_backend_mod_count(object@files, object@modCount))
    if (is.null(msg)) TRUE
    else msg
})

#' @rdname hidden_aliases
setMethod("backendInitialize", signature = "MsBackendMemory",
          definition = function(object, files, spectraData, ...) {
              if (missing(files)) files <- character()
              if (missing(spectraData)) spectraData <- DataFrame()
              object@files <- files
              object@modCount <- integer(length(files))
              object@spectraData <- spectraData
              validObject(object)
              object
          })

## Data accessors

#' @rdname hidden_aliases
setMethod("acquisitionNum", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "acquisitionNum"))
        object@spectraData$acquisitionNum
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("centroided", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "centroided"))
        object@spectraData$centroided
    else rep(NA, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("centroided", "MsBackendMemory", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.logical(value) | length(value) != length(object))
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    object@spectraData$centroided <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("collisionEnergy", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "collisionEnergy"))
        object@spectraData$collisionEnergy
    else rep(NA_real_, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("collisionEnergy", "MsBackendMemory", function(object, value) {
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$collisionEnergy <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("fromFile", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "fromFile"))
        object@spectraData$fromFile
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("intensity", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "intensity"))
        object@spectraData$intensity
    else {
        lst <- list(numeric)
        lst[rep(1, times = length(object))]
    }
})

#' @rdname hidden_aliases
setMethod("ionCount", "MsBackendMemory", function(object) {
    lapply(intensity(object), sum, na.rm = TRUE)
})

#' @rdname hidden_aliases
setMethod("isCentroided", "MsBackendMemory", function(object, ...) {
    ## TODO @jo IMPLEMENT
    stop("Not implemented for ", class(object), ".")
})

#' @rdname hidden_aliases
setMethod("isEmpty", "MsBackendMemory", function(x) {
    lengths(intensity(object)) == 0
})

#' @rdname hidden_aliases
setMethod("length", "MsBackendMemory", function(x) {
    nrow(x@spectraData)
})

#' @rdname hidden_aliases
setMethod("msLevel", "MsBackendMemory", function(object, ...) {
    if (any(colnames(object@spectraData) == "msLevel"))
        object@spectraData$msLevel
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("mz", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "mz"))
        object@spectraData$mz
    else {
        lst <- list(numeric)
        lst[rep(1, times = length(object))]
    }
})

#' @rdname hidden_aliases
setMethod("peaks", "MsBackendMemory", function(object) {
    mapply(mz(object), intensity(object), FUN = function(m, i)
        data.frame(mz = mz, intensity = i))
})

#' @rdname hidden_aliases
setMethod("peaksCount", "MsBackendMemory", function(object) {
    lengths(mz(object))
})

#' @rdname hidden_aliases
setMethod("polarity", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "polarity"))
        object@spectraData$polarity
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("polarity", "MsBackendMemory", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be an 'integer' of length 1 or ", length(object))
    object@spectraData$polarity <- as.integer(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("precScanNum", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "precScanNum"))
        object@spectraData$precScanNum
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("precursorCharge", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "precursorCharge"))
        object@spectraData$precursorCharge
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("precursorIntensity", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "precursorIntensity"))
        object@spectraData$precursorIntensity
    else rep(NA_real_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("precursorMz", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "precursorMz"))
        object@spectraData$precursorMz
    else rep(NA_real_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("rtime", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "rt"))
        object@spectraData$rt
    else rep(NA_real_, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("rtime", "MsBackendMemory", function(object, value) {
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$polarity <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("scanIndex", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "scanIndex"))
        object@spectraData$scanIndex
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("smoothed", "MsBackendMemory", function(object) {
    if (any(colnames(object@spectraData) == "smoothed"))
        object@spectraData$smoothed
    else rep(NA, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("smoothed", "MsBackendMemory", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.logical(value) | length(value) != length(object))
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    object@spectraData$smoothed <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("spectraData", "MsBackendMemory", function(object, columns) {
    ## TODO @jo IMPLEMENT:
    ## allow to specify columns. Have then to check if these columns exist,
    ## if not call the respective method.
})

#' @rdname hidden_aliases
setReplaceMethod("spectraData", "MsBackendMemory", function(object, value) {
    if (!is(value, "DataFrame") || rownames(value) != length(object))
        stop("'value' has to be a 'DataFrame' with ", length(object), " rows.")
    object@spectraData <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("spectraNames", "MsBackendMemory", function(object) {
    rownames(object@spectraData)
})

#' @rdname hidden_aliases
setReplaceMethod("spectraNames", "MsBackendMemory", function(object, value) {
    rownames(object@spectraData) <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("spectraVariables", "MsBackendMemory", function(object) {
    ## TODO @jo IMPLEMENT:
    ## list the colnames of spectraData AND all columns supported by accessor
    ## methods.
})

#' @rdname hidden_aliases
setMethod("tic", "MsBackendMemory", function(object, initial = TRUE) {
    if (initial) {
        if (any(colnames(object@spectraData) == "totIonCurrent"))
            object@spectraData$totIonCurrent
        else rep(NA_real_, times = length(object))
    } else vapply(intensity(object), sum, numeric(1), na.rm = TRUE)
})
