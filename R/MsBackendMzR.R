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
    cat("validObject MsBackendMzR\n")
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
                             spectraData = .compress_spectra_data(spectraData),
                             ...)
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
    object@spectraData$centroided <- .rle_compress(value)
    validObject(object)
    object
})

## #' @rdname hidden_aliases
## setMethod("collisionEnergy", "MsBackendMzR", function(object) {
##     .get_rle_column(object@spectraData, "collisionEnergy")
## })

## #' @rdname hidden_aliases
## setReplaceMethod("collisionEnergy", "MsBackendMzR", function(object, value) {
##     if (!is.numeric(value) | length(value) != length(object))
##         stop("'value' has to be a 'numeric' of length ", length(object))
##     object@spectraData$collisionEnergy <- value
##     validObject(object)
##     object
## })

#' @rdname hidden_aliases
setMethod("fromFile", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "fromFile")
})

#' @rdname hidden_aliases
setMethod("intensity", "MsBackendMzR", function(object) {
    stop("Have to read intensity from the original files.")
})

## #' @rdname hidden_aliases
## setMethod("ionCount", "MsBackendMzR", function(object) {
##     vapply(intensity(object), sum, numeric(1), na.rm = TRUE)
## })

## #' @rdname hidden_aliases
## setMethod("isCentroided", "MsBackendMzR", function(object, ...) {
##     vapply(peaks(object), .isCentroided, logical(1))
## })

## #' @rdname hidden_aliases
## setMethod("isEmpty", "MsBackendMzR", function(x) {
##     lengths(intensity(x)) == 0
## })

## #' @rdname hidden_aliases
## setMethod("length", "MsBackendMzR", function(x) {
##     nrow(x@spectraData)
## })

#' @rdname hidden_aliases
setMethod("msLevel", "MsBackendMzR", function(object, ...) {
    .get_rle_column(object@spectraData, "msLevel")
})

#' @rdname hidden_aliases
setMethod("mz", "MsBackendMzR", function(object) {
    stop("Have to read mz from the original files.")
})

## #' @rdname hidden_aliases
## setMethod("peaks", "MsBackendMzR", function(object) {
##     mapply(mz(object), intensity(object), FUN = function(m, i)
##         cbind(mz = m, intensity = i), SIMPLIFY = FALSE, USE.NAMES = FALSE)
## })

## #' @rdname hidden_aliases
## setMethod("peaksCount", "MsBackendMzR", function(object) {
##     lengths(mz(object))
## })

#' @rdname hidden_aliases
setMethod("polarity", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "polarity")
})

## #' @rdname hidden_aliases
## setReplaceMethod("polarity", "MsBackendMzR", function(object, value) {
##     if (length(value) == 1)
##         value <- rep(value, length(object))
##     if (!is.numeric(value) | length(value) != length(object))
##         stop("'value' has to be an 'integer' of length 1 or ", length(object))
##     object@spectraData$polarity <- as.integer(value)
##     validObject(object)
##     object
## })

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

## #' @rdname hidden_aliases
## setReplaceMethod("rtime", "MsBackendMzR", function(object, value) {
##     if (!is.numeric(value) | length(value) != length(object))
##         stop("'value' has to be a 'numeric' of length ", length(object))
##     object@spectraData$rtime <- value
##     validObject(object)
##     object
## })

#' @rdname hidden_aliases
setMethod("scanIndex", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "scanIndex")
})

#' @rdname hidden_aliases
setMethod("smoothed", "MsBackendMzR", function(object) {
    .get_rle_column(object@spectraData, "smoothed")
})

## #' @rdname hidden_aliases
## setReplaceMethod("smoothed", "MsBackendMzR", function(object, value) {
##     if (length(value) == 1)
##         value <- rep(value, length(object))
##     if (!is.logical(value) | length(value) != length(object))
##         stop("'value' has to be a 'logical' of length 1 or ", length(object))
##     object@spectraData$smoothed <- value
##     validObject(object)
##     object
## })

## #' @rdname hidden_aliases
## #'
## #' @importFrom methods as
## setMethod("spectraData", "MsBackendMzR",
##           function(object, columns = spectraVariables(object)) {
##               cn <- colnames(object@spectraData)
##               if(!nrow(object@spectraData)) {
##                   res <- lapply(.SPECTRA_DATA_COLUMNS, do.call, args = list())
##                   res <- DataFrame(res)
##                   res$mz <- list()
##                   res$intensity <- list()
##                   return(res[, columns, drop = FALSE])
##               }
##               not_found <- setdiff(columns, colnames(object@spectraData))
##               if (length(not_found))
##                   stop("Column(s) ", paste(not_found, collapse = ", "),
##                        " not available")
##               .uncompress_spectra_data(object@spectraData[, columns,
##                                                           drop = FALSE])
##           })

## #' @rdname hidden_aliases
## setReplaceMethod("spectraData", "MsBackendMzR", function(object,
##                                                                   value) {
##     if (!is(value, "DataFrame") || nrow(value) != length(object))
##         stop("'value' has to be a 'DataFrame' with ", length(object), " rows.")
##     object@spectraData <- .compress_spectra_data(value)
##     validObject(object)
##     object
## })

#' @rdname hidden_aliases
setMethod("spectraNames", "MsBackendMzR", function(object) {
    rownames(object@spectraData)
})

## #' @rdname hidden_aliases
## setReplaceMethod("spectraNames", "MsBackendMzR",
##                  function(object, value) {
##                      rownames(object@spectraData) <- value
##                      validObject(object)
##                      object
##                  })

#' @rdname hidden_aliases
setMethod("spectraVariables", "MsBackendMzR", function(object) {
    unique(c(names(.SPECTRA_DATA_COLUMNS), colnames(object@spectraData)))
})

## #' @rdname hidden_aliases
## setMethod("tic", "MsBackendMzR", function(object, initial = TRUE) {
##     if (initial) {
##         if (any(colnames(object@spectraData) == "totIonCurrent"))
##             .get_rle_column(object@spectraData, "totIonCurrent")
##         else rep(NA_real_, times = length(object))
##     } else vapply(intensity(object), sum, numeric(1), na.rm = TRUE)
## })

## #' @rdname hidden_aliases
## setMethod("[", "MsBackendMzR", function(x, i, j, ..., drop = FALSE) {
##     if (!missing(j))
##         stop("Subsetting byt column ('j = ", j, "' is not supported")
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
