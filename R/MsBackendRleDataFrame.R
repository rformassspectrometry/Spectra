#' @include hidden_aliases.R
NULL

#' @title In-memory MS data backend
#'
#' @description
#'
#' Columns containing exclusively `NA` values are converted to `Rle` in the
#' `spectraData` `DataFrame` to keep memory demand to a minimum.
#'
#' @name MsBackendRleDataFrame
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @noRd
#'
#' @exportClass MsBackendRleDataFrame
NULL

setClass("MsBackendRleDataFrame",
         contains = "MsBackend",
         slots = c(spectraData = "DataFrame"),
         prototype = prototype(spectraData = DataFrame(),
                               files = character(),
                               modCount = integer(),
                               readonly = FALSE,
                               version = "0.1"))

setValidity("MsBackendRleDataFrame", function(object) {
    msg <- .valid_spectra_data_required_columns(object@spectraData)
    if (length(msg))
        return(msg)
    msg <- c(
        .valid_column_rle_datatype(object@spectraData, .SPECTRA_DATA_COLUMNS),
        .valid_intensity_column(object@spectraData),
        .valid_mz_column(object@spectraData),
        .valid_ms_backend_files(object@files),
        .valid_ms_backend_files_from_file(
            object@files, as.vector(object@spectraData$fromFile)),
        .valid_ms_backend_mod_count(object@files, object@modCount),
        .valid_intensity_mz_columns(object))
    if (is.null(msg)) TRUE
    else msg
})

#' @rdname hidden_aliases
setMethod("show", "MsBackendRleDataFrame", function(object) {
    spd <- spectraData(object, c("msLevel", "rtime", "scanIndex"))
    cat(class(object), "with", nrow(spd), "spectra\n")
    if (nrow(spd)) {
        txt <- capture.output(show(spd))
        cat(txt[-1], sep = "\n")
        sp_cols <- spectraVariables(object)
        cat(" ...", length(sp_cols) - 3, "more variables/columns.\n")
    }
})

#' @rdname hidden_aliases
setMethod("backendInitialize", signature = "MsBackendRleDataFrame",
          definition = function(object, files, spectraData, ...) {
              if (missing(files)) files <- character()
              if (missing(spectraData)) spectraData <- DataFrame()
              object@files <- files
              object@modCount <- integer(length(files))
              object@spectraData <- .initialize_spectra_data(spectraData)
              validObject(object)
              object
          })

#' @importMethodsFrom S4Vectors [[
.get_rle_column <- function(x, column) {
    if (!any(colnames(x) == column))
        return(do.call(.SPECTRA_DATA_COLUMNS[column], args = list()))
    if (is(x[[column]], "Rle"))
        as.vector(x[[column]])
    else x[[column]]
}

#' @rdname hidden_aliases
setMethod("acquisitionNum", "MsBackendRleDataFrame", function(object) {
    .get_rle_column(object@spectraData, "acquisitionNum")
})

#' @rdname hidden_aliases
setMethod("centroided", "MsBackendRleDataFrame", function(object) {
    .get_rle_column(object@spectraData, "centroided")
})

#' @rdname hidden_aliases
setReplaceMethod("centroided", "MsBackendRleDataFrame", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.logical(value) | length(value) != length(object))
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    object@spectraData$centroided <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("collisionEnergy", "MsBackendRleDataFrame", function(object) {
    .get_rle_column(object@spectraData, "collisionEnergy")
})

#' @rdname hidden_aliases
setReplaceMethod("collisionEnergy", "MsBackendRleDataFrame", function(object, value) {
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$collisionEnergy <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("fromFile", "MsBackendRleDataFrame", function(object) {
    .get_rle_column(object@spectraData, "fromFile")
})

#' @rdname hidden_aliases
setMethod("intensity", "MsBackendRleDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "intensity"))
        object@spectraData$intensity
    else list()
})

#' @rdname hidden_aliases
setMethod("ionCount", "MsBackendRleDataFrame", function(object) {
    vapply(intensity(object), sum, numeric(1), na.rm = TRUE)
})

#' @rdname hidden_aliases
setMethod("isCentroided", "MsBackendRleDataFrame", function(object, ...) {
    vapply(peaks(object), .isCentroided, logical(1))
})

#' @rdname hidden_aliases
setMethod("isEmpty", "MsBackendRleDataFrame", function(x) {
    lengths(intensity(x)) == 0
})

#' @rdname hidden_aliases
setMethod("length", "MsBackendRleDataFrame", function(x) {
    nrow(x@spectraData)
})

#' @rdname hidden_aliases
setMethod("msLevel", "MsBackendRleDataFrame", function(object, ...) {
    .get_rle_column(object@spectraData, "msLevel")
})

#' @rdname hidden_aliases
setMethod("mz", "MsBackendRleDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "mz"))
        object@spectraData$mz
    else list()
})

#' @rdname hidden_aliases
setMethod("peaks", "MsBackendRleDataFrame", function(object) {
    mapply(mz(object), intensity(object), FUN = function(m, i)
        cbind(mz = m, intensity = i), SIMPLIFY = FALSE, USE.NAMES = FALSE)
})

#' @rdname hidden_aliases
setMethod("peaksCount", "MsBackendRleDataFrame", function(object) {
    lengths(mz(object))
})

#' @rdname hidden_aliases
setMethod("polarity", "MsBackendRleDataFrame", function(object) {
    .get_rle_column(object@spectraData, "polarity")
})

#' @rdname hidden_aliases
setReplaceMethod("polarity", "MsBackendRleDataFrame", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be an 'integer' of length 1 or ", length(object))
    object@spectraData$polarity <- as.integer(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("precScanNum", "MsBackendRleDataFrame", function(object) {
    .get_rle_column(object@spectraData, "precScanNum")
})

#' @rdname hidden_aliases
setMethod("precursorCharge", "MsBackendRleDataFrame", function(object) {
    .get_rle_column(object@spectraData, "precursorCharge")
})

#' @rdname hidden_aliases
setMethod("precursorIntensity", "MsBackendRleDataFrame", function(object) {
    .get_rle_column(object@spectraData, "precursorIntensity")
})

#' @rdname hidden_aliases
setMethod("precursorMz", "MsBackendRleDataFrame", function(object) {
    .get_rle_column(object@spectraData, "precursorMz")
})

#' @rdname hidden_aliases
setMethod("rtime", "MsBackendRleDataFrame", function(object) {
    .get_rle_column(object@spectraData, "rtime")
})

#' @rdname hidden_aliases
setReplaceMethod("rtime", "MsBackendRleDataFrame", function(object, value) {
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$rtime <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("scanIndex", "MsBackendRleDataFrame", function(object) {
    .get_rle_column(object@spectraData, "scanIndex")
})

#' @rdname hidden_aliases
setMethod("smoothed", "MsBackendRleDataFrame", function(object) {
    .get_rle_column(object@spectraData, "smoothed")
})

#' @rdname hidden_aliases
setReplaceMethod("smoothed", "MsBackendRleDataFrame", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.logical(value) | length(value) != length(object))
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    object@spectraData$smoothed <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
#'
#' @importFrom methods as
setMethod("spectraData", "MsBackendRleDataFrame",
          function(object, columns = spectraVariables(object)) {
              cn <- colnames(object@spectraData)
              if(!nrow(object@spectraData)) {
                  res <- lapply(.SPECTRA_DATA_COLUMNS, do.call, args = list())
                  res <- DataFrame(res)
                  res$mz <- list()
                  res$intensity <- list()
                  return(res[, columns, drop = FALSE])
              }
              not_found <- setdiff(columns, colnames(object@spectraData))
              if (length(not_found))
                  stop("Column(s) ", paste(not_found, collapse = ", "),
                       " not available")
              .uncompress_spectra_data(object@spectraData[, columns,
                                                          drop = FALSE])
          })

#' @rdname hidden_aliases
setReplaceMethod("spectraData", "MsBackendRleDataFrame", function(object,
                                                                  value) {
    if (!is(value, "DataFrame") || nrow(value) != length(object))
        stop("'value' has to be a 'DataFrame' with ", length(object), " rows.")
    object@spectraData <- .compress_spectra_data(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("spectraNames", "MsBackendRleDataFrame", function(object) {
    rownames(object@spectraData)
})

#' @rdname hidden_aliases
setReplaceMethod("spectraNames", "MsBackendRleDataFrame",
                 function(object, value) {
                     rownames(object@spectraData) <- value
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("spectraVariables", "MsBackendRleDataFrame", function(object) {
    cn <- colnames(object@spectraData)
    if (length(cn))
        cn
    else names(.SPECTRA_DATA_COLUMNS)
})

#' @rdname hidden_aliases
setMethod("tic", "MsBackendRleDataFrame", function(object, initial = TRUE) {
    if (initial) {
        if (any(colnames(object@spectraData) == "totIonCurrent"))
            .get_rle_column(object@spectraData, "totIonCurrent")
        else rep(NA_real_, times = length(object))
    } else vapply(intensity(object), sum, numeric(1), na.rm = TRUE)
})

#' @rdname hidden_aliases
setMethod("[", "MsBackendRleDataFrame", function(x, i, j, ..., drop = FALSE) {
    if (!missing(j))
        stop("Subsetting byt column ('j = ", j, "' is not supported")
    i <- .i_to_index(i, length(x), rownames(x@spectraData))
    x@spectraData <- x@spectraData[i, , drop = FALSE]
    orig_files <- x@files
    files_idx <- unique(fromFile(x))
    x@files <- orig_files[files_idx]
    x@modCount <- x@modCount[files_idx]
    x@spectraData$fromFile <- match(orig_files[fromFile(x)], x@files)
    validObject(x)
    x
})
