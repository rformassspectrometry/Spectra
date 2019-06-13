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
    msg <- c(msg, .valid_ms_backend_files_exist(object@files))
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
              if (!all(file.exists(files)))
                  stop("File(s) ", paste(files[!file.exists(files)]),
                       " not found")
              msg <- .valid_ms_backend_files(files)
              if (length(msg))
                  stop(msg)
              if (!missing(spectraData)) {
                  spectraData$mz <- NULL
                  spectraData$intensity <- NULL
              } else {
                  spectraData <- do.call(
                      rbind, bpmapply(files, seq_along(files),
                                      FUN = function(fl, index) {
                                          cbind(Spectra:::.mzR_header(fl),
                                                fromFile = index)
                                      }, BPPARAM = BPPARAM))
              }
              callNextMethod(object = object, files = files,
                             spectraData = spectraData,
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
setMethod("intensity", "MsBackendMzR", function(object) {
    NumericList(lapply(peaks(object), "[", , 2), compress = FALSE)
})

#' @rdname hidden_aliases
setReplaceMethod("intensity", "MsBackendMzR", function(object, value) {
    stop(class(object), " does not support replacing intensity values")
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
setMethod("mz", "MsBackendMzR", function(object) {
    NumericList(lapply(peaks(object), "[", , 1), compress = FALSE)
})

#' @rdname hidden_aliases
setReplaceMethod("mz", "MsBackendMzR", function(object, value) {
    stop(class(object), " does not support replacing m/z values")
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
    as.integer(lengths(peaks(object)) / 2L)
})

#' @rdname hidden_aliases
#'
#' @importFrom methods as
setMethod("spectraData", "MsBackendMzR",
          function(object, columns = spectraVariables(object)) {
              .spectra_data_mzR(object, columns)
          })

#' @rdname hidden_aliases
setReplaceMethod("spectraData", "MsBackendMzR", function(object, value) {
    if (inherits(value, "DataFrame") && any(colnames(value) %in%
                                            c("mz", "intensity"))) {
        warning("Ignoring columns \"mz\" and \"intensity\" as the ",
                "'MzBackendMzR' backend currently does not support replacing ",
                "them.")
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

#' @rdname hidden_aliases
setReplaceMethod("$", "MsBackendMzR", function(x, name, value) {
    if (name == "mz" || name == "intensity")
        stop("'MsBackendMzR' does not support replacing mz or intensity values")
    if (length(value) == 1)
        value <- rep(value, length(x))
    if (length(value) != length(x))
        stop("Length of 'value' has to be either 1 or ", length(x))
    x@spectraData[[name]] <- .as_rle(value)
    validObject(x)
    x
})
