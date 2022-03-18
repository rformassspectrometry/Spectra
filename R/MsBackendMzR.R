#' @include hidden_aliases.R
NULL

#' @title mzR-based backend
#'
#' @description
#'
#' The `MsBackendMzR` inherits all slots and methods from the base
#' `MsBackendDataFrame` (in-memory) backend. It overrides the base `mz` and
#' `intensity` methods as well as `peaksData` to read the respective data from
#' the original raw data files.
#'
#' The validator function has to ensure that the files exist and that required
#' column names are present.
#'
#' The `backendInitialize` method reads the header data from the raw files and
#' hence fills the `spectraData` slot.
#'
#' @author Johannes Rainer
#'
#' @noRd
setClass("MsBackendMzR",
         contains = "MsBackendDataFrame",
         prototype = prototype(version = "0.1", readonly = TRUE))

setValidity("MsBackendMzR", function(object) {
    msg <- .valid_spectra_data_required_columns(object@spectraData,
                                                c("dataStorage", "scanIndex"))
    msg <- c(msg, .valid_ms_backend_files_exist(
                      unique(object@spectraData$dataStorage)))
    if (length(msg)) msg
    else TRUE
})

#' @rdname hidden_aliases
#'
#' @importFrom methods callNextMethod
#'
#' @importFrom MsCoreUtils rbindFill
#'
#' @importMethodsFrom BiocParallel bpmapply
#'
#' @importFrom BiocParallel bpparam
setMethod("backendInitialize", "MsBackendMzR",
          function(object, files, ..., BPPARAM = bpparam()) {
              if (missing(files) || !length(files))
                  stop("Parameter 'files' is mandatory for 'MsBackendMzR'")
              if (!is.character(files))
                  stop("Parameter 'files' is expected to be a character vector",
                       " with the files names from where data should be",
                       " imported")
              files <- normalizePath(files, mustWork = FALSE)
              msg <- .valid_ms_backend_files_exist(files)
              if (length(msg))
                  stop(msg)
              spectraData <- do.call(
                  rbindFill, bplapply(files,
                                  FUN = function(fl) {
                                      cbind(Spectra:::.mzR_header(fl),
                                            dataStorage = fl)
                                  }, BPPARAM = BPPARAM))
              spectraData$dataOrigin <- spectraData$dataStorage
              object@spectraData <- spectraData
              validObject(object)
              object
          })

#' @rdname hidden_aliases
setMethod("show", "MsBackendMzR", function(object) {
    callNextMethod()
    fls <- unique(object@spectraData$dataStorage)
    if (length(fls)) {
        to <- min(3, length(fls))
        cat("\nfile(s):\n", paste(basename(fls[seq_len(to)]), collapse = "\n"),
            "\n", sep = "")
        if (length(fls) > 3)
            cat(" ...", length(fls) - 3, "more files\n")
    }
})

#' @rdname hidden_aliases
setMethod(
    "peaksData", "MsBackendMzR",
    function(object, columns = peaksVariables(object)) {
        if (!length(object))
            return(list())
        if (!all(columns %in% c("mz", "intensity")))
            stop("'peaksData' for 'MsBackendMzR' does only support",
                 " columns \"mz\" and \"intensity\"", call. = FALSE)
        fls <- unique(object@spectraData$dataStorage)
        read_header <- getOption("READ_HEADER", FALSE)
        if (length(fls) > 1) {
            f <- factor(dataStorage(object), levels = fls)
            unsplit(mapply(FUN = .mzR_peaks, fls, split(scanIndex(object), f),
                           MoreArgs = list(readHeader = read_header,
                                           columns = columns),
                           SIMPLIFY = FALSE, USE.NAMES = FALSE), f)
        } else
            .mzR_peaks(fls, scanIndex(object), readHeader = read_header,
                       columns = columns)
    })

#' @rdname hidden_aliases
setMethod("intensity", "MsBackendMzR", function(object) {
    NumericList(lapply(peaksData(object), "[", , 2), compress = FALSE)
})

#' @rdname hidden_aliases
setReplaceMethod("intensity", "MsBackendMzR", function(object, value) {
    stop(class(object), " does not support replacing intensity values")
})

#' @rdname hidden_aliases
setMethod("ionCount", "MsBackendMzR", function(object) {
    vapply1d(peaksData(object), function(z) sum(z[, 2], na.rm = TRUE))
})

#' @rdname hidden_aliases
setMethod("isCentroided", "MsBackendMzR", function(object, ...) {
    vapply1l(peaksData(object), .peaks_is_centroided)
})

#' @rdname hidden_aliases
setMethod("isEmpty", "MsBackendMzR", function(x) {
    lengths(x) == 0
})

#' @rdname hidden_aliases
setMethod("lengths", "MsBackendMzR", function(x, use.names = FALSE) {
    as.integer(lengths(peaksData(x)) / 2L)
})

#' @rdname hidden_aliases
setMethod("mz", "MsBackendMzR", function(object) {
    NumericList(lapply(peaksData(object), "[", , 1), compress = FALSE)
})

#' @rdname hidden_aliases
setReplaceMethod("mz", "MsBackendMzR", function(object, value) {
    stop(class(object), " does not support replacing m/z values")
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
    object@spectraData <- value
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
setReplaceMethod("$", "MsBackendMzR", function(x, name, value) {
    if (name == "mz" || name == "intensity")
        stop("'MsBackendMzR' does not support replacing mz or intensity values")
    value_len <- length(value)
    if (value_len == 1L || value_len == length(x))
        x@spectraData[[name]] <- value
    else
        stop("Length of 'value' has to be either 1 or ", length(x))
    validObject(x)
    x
})

#' @rdname hidden_aliases
setMethod("export", "MsBackendMzR", function(object, x, file = tempfile(),
                                             format = c("mzML", "mzXML"),
                                             copy = FALSE,
                                             BPPARAM = bpparam()) {
    l <- length(x)
    file <- sanitize_file_name(file)
    if (length(file) == 1)
        file <- rep_len(file, l)
    if (length(file) != l)
        stop("Parameter 'file' has to be either of length 1 or ",
             length(x), ", i.e. 'length(x)'.", call. = FALSE)
    f <- factor(file, levels = unique(file))
    tmp <- bpmapply(.write_ms_data_mzR, split(x, f), levels(f),
                    MoreArgs = list(format = format, copy = copy),
                    BPPARAM = BPPARAM)
})
