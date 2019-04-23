#' @include hidden_aliases.R
NULL

#' @title mzR-based backend
#'
#' @description
#'
#' The `MsBackendMzR` inherits all slots and methods from the base
#' `MsBackendDataFrame` (in-memory) backend. It overrides the base `mz` and
#' `intensity` methods as well as `peaks` to read the respective data from
#' the original raw data files.
#'
#' The validator function has to ensure that the files provided in the
#' `files` slot exist.
#'
#' The `backendInitialize` method reads the header data from the raw files and
#' hence fills the `spectraData` slot. Note that this method could be called
#' several times, e.g. also to *re-fill* `spectraData` after dropping some of
#' its columns.
#'
#' @author Johannes Rainer
#'
#' @noRd
setClass("MsBackendMzR",
         contains = "MsBackendDataFrame",
         prototype = prototype(version = "0.1", readonly = TRUE))

setValidity("MsBackendMzR", function(object) {
    msg <- .valid_spectra_data_required_columns(object@spectraData,
                                                c("fromFile", "scanIndex"))
    if (length(msg)) msg
    else TRUE
})

#' @rdname hidden_aliases
#'
#' @importFrom methods callNextMethod
#'
#' @importMethodsFrom BiocParallel bpmapply
#'
#' @importFrom BiocParallel bpparam
setMethod("backendInitialize", "MsBackendMzR",
          function(object, files, spectraData, ..., BPPARAM = bpparam()) {
              if (missing(files) || !length(files))
                  stop("Parameter 'files' is mandatory for 'MsBackendMzR'")
              files <- normalizePath(files)
              msg <- .valid_ms_backend_files(files)
              if (length(msg))
                  stop(msg)
              spectraData <- do.call(
                  rbind, bpmapply(files, seq_along(files),
                                  FUN = function(fl, index) {
                                      cbind(Spectra:::.mzR_header(fl),
                                            fromFile = index)
                                  }))
              callNextMethod(object = object, files = files,
                             spectraData = .as_rle_spectra_data(spectraData),
                             ...)
          })

#' @rdname hidden_aliases
setMethod("show", "MsBackendMzR", function(object) {
    callNextMethod()
    fls <- basename(object@files)
    if (length(fls)) {
        to <- min(3, length(fls))
        cat("\nfile(s):\n", paste(basename(fls[1:to]), collapse = "\n"),
            "\n", sep = "")
        if (length(fls) > 3)
            cat(" ...", length(fls) - 3, "more files\n")
    }
})

#' @rdname hidden_aliases
setMethod("acquisitionNum", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "acquisitionNum")
})

#' @rdname hidden_aliases
setMethod("centroided", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "centroided")
})

#' @rdname hidden_aliases
setReplaceMethod("centroided", "MsBackendMzR", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.logical(value) | length(value) != length(object))
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    object@spectraData$centroided <- .as_rle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("collisionEnergy", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "collisionEnergy")
})

#' @rdname hidden_aliases
setReplaceMethod("collisionEnergy", "MsBackendMzR", function(object, value) {
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$collisionEnergy <- .as_rle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("fromFile", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "fromFile")
})

#' @rdname hidden_aliases
setMethod("intensity", "MsBackendMzR", function(object) {
    SimpleList(lapply(peaks(object), function(z) z[, 2]))
})

#' @rdname hidden_aliases
setMethod("ionCount", "MsBackendMzR", function(object) {
    vapply(peaks(object), function(z) sum(z[, 2], na.rm = TRUE), numeric(1))
})

#' @rdname hidden_aliases
setMethod("isCentroided", "MsBackendMzR", function(object, ...) {
    vapply(peaks(object), .isCentroided, logical(1))
})

#' @rdname hidden_aliases
setMethod("isEmpty", "MsBackendMzR", function(x) {
    peaksCount(x) == 0
})

#' @rdname hidden_aliases
setMethod("msLevel", "MsBackendMzR", function(object, ...) {
    .get_rle_column(object@spectraData, "msLevel")
})

#' @rdname hidden_aliases
setMethod("mz", "MsBackendMzR", function(object) {
    SimpleList(lapply(peaks(object), function(z) z[, 1]))
})

#' @rdname hidden_aliases
setMethod("peaks", "MsBackendMzR", function(object) {
    if (!length(object))
        return(list())
    if (length(object@files) > 1) {
            unlist(mapply(FUN = .mzR_peaks, object@files,
                          split(object@spectraData$scanIndex,
                                as.factor(object@spectraData$fromFile)),
                          SIMPLIFY = FALSE, USE.NAMES = FALSE),
                   use.names = FALSE, recursive = FALSE)
    } else
        .mzR_peaks(object@files, object@spectraData$scanIndex)
})

#' @rdname hidden_aliases
setMethod("peaksCount", "MsBackendMzR", function(object) {
    vapply(peaks(object), nrow, integer(1))
})

#' @rdname hidden_aliases
setMethod("polarity", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "polarity")
})

#' @rdname hidden_aliases
setReplaceMethod("polarity", "MsBackendMzR", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be an 'integer' of length 1 or ", length(object))
    object@spectraData$polarity <- .as_rle(as.integer(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("precScanNum", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "precScanNum")
})

#' @rdname hidden_aliases
setMethod("precursorCharge", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "precursorCharge")
})

#' @rdname hidden_aliases
setMethod("precursorIntensity", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "precursorIntensity")
})

#' @rdname hidden_aliases
setMethod("precursorMz", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "precursorMz")
})

#' @rdname hidden_aliases
setMethod("rtime", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "rtime")
})

#' @rdname hidden_aliases
setReplaceMethod("rtime", "MsBackendMzR", function(object, value) {
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$rtime <- .as_rle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("scanIndex", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "scanIndex")
})

#' @rdname hidden_aliases
setMethod("smoothed", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "smoothed")
})

#' @rdname hidden_aliases
setReplaceMethod("smoothed", "MsBackendMzR", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.logical(value) | length(value) != length(object))
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    object@spectraData$smoothed <- .as_rle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
#'
#' @importFrom methods as
setMethod("spectraData", "MsBackendMzR",
          function(object, columns = spectraVariables(object)) {
              cn <- colnames(object@spectraData)
              if(!nrow(object@spectraData)) {
                  res <- lapply(.SPECTRA_DATA_COLUMNS, do.call, args = list())
                  res <- DataFrame(res)
                  res$mz <- SimpleList()
                  res$intensity <- SimpleList()
                  return(res[, columns, drop = FALSE])
              }
              not_found <- setdiff(columns, c(cn, names(.SPECTRA_DATA_COLUMNS)))
              if (length(not_found))
                  stop("Column(s) ", paste(not_found, collapse = ", "),
                       " not available")
              sp_cols <- columns[columns %in% cn]
              res <- .as_vector_spectra_data(
                  object@spectraData[, sp_cols, drop = FALSE])
              if (any(columns %in% c("mz", "intensity"))) {
                  pks <- peaks(object)
                  if (any(columns == "mz"))
                      res$mz <- SimpleList(lapply(pks, function(z) z[, 1]))
                  if (any(columns == "intensity"))
                      res$intensity <- SimpleList(lapply(pks, function(z) z[, 2]))
              }
              other_cols <- setdiff(
                  columns[!(columns %in% c("mz", "intensity"))], sp_cols)
              if (length(other_cols)) {
                  other_res <- lapply(other_cols, .get_spectra_data_column,
                                      x = object)
                  names(other_res) <- other_cols
                  res <- cbind(res, as(other_res, "DataFrame"))
              }
              res[, columns, drop = FALSE]
          })

#' @rdname hidden_aliases
setReplaceMethod("spectraData", "MsBackendMzR", function(object, value) {
    if (!is(value, "DataFrame") || nrow(value) != length(object))
        stop("'value' has to be a 'DataFrame' with ", length(object), " rows.")
    if (any(colnames(value) %in% c("mz", "intensity"))) {
        value <- value[, !(colnames(value) %in% c("mz", "intensity")),
                       drop = FALSE]
    }
    object@spectraData <- .as_rle_spectra_data(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("spectraNames", "MsBackendMzR", function(object) {
    rownames(object@spectraData)
})

#' @rdname hidden_aliases
setReplaceMethod("spectraNames", "MsBackendMzR", function(object, value) {
    rownames(object@spectraData) <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("spectraVariables", "MsBackendMzR", function(object) {
    unique(c(names(.SPECTRA_DATA_COLUMNS), colnames(object@spectraData)))
})

#' @rdname hidden_aliases
setMethod("tic", "MsBackendMzR", function(object, initial = TRUE) {
    if (initial) {
        if (any(colnames(object@spectraData) == "totIonCurrent"))
            .get_rle_column(object@spectraData, "totIonCurrent")
        else rep(NA_real_, times = length(object))
    } else vapply(intensity(object), sum, numeric(1), na.rm = TRUE)
})

## #' @rdname hidden_aliases
## setMethod("[", "MsBackendMzR", function(x, i, j, ..., drop = FALSE) {
##     if (!missing(j))
##         stop("Subsetting by column ('j = ", j, "' is not supported")
##     i <- .i_to_index(i, length(x), rownames(x@spectraData))
##     x@spectraData <- x@spectraData[i, , drop = FALSE]
##     orig_files <- x@files
##     files_idx <- unique(fromFile(x))
##     x@files <- orig_files[files_idx]
##     x@modCount <- x@modCount[files_idx]
##     x@spectraData$fromFile <- match(orig_files[fromFile(x)], x@files)
##     validObject(x)
##     x
## })
