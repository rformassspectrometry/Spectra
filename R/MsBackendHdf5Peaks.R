#' @include hidden_aliases.R
NULL

#' @title Hdf5-based backend
#'
#' @description
#'
#' The `MsBackendHdf5Peaks` is a bakend that keeps general spectra variables in
#' memory while reading (writing) peak data (i.e. m/z and intensity values) from
#' and to Hdf5 files.
#'
#' @note
#'
#' For memory issues we might want to extend a MsBackendRleDataFrame instead,
#' which could be the base class for both the MsBackendMzR and the
#' MsBackendHdf5Peaks.
#'
#' @author Johannes Rainer
#'
#' @noRd
setClass("MsBackendHdf5Peaks",
         contains = "MsBackendDataFrame",
         slots = c(
             modCount = "integer"
         ),
         prototype = prototype(version = "0.1", readonly = FALSE))

setValidity("MsBackendHdf5Peaks", function(object) {
    msg <- .valid_spectra_data_required_columns(object@spectraData,
                                                c("dataStorage", "scanIndex"))
    fls <- unique(object@spectraData$dataStorage)
    msg <- c(msg, .valid_ms_backend_mod_count(object@modCount, fls))
    msg <- c(msg, .valid_ms_backend_files_exist(fls))
    msg <- c(msg, .valid_h5files(fls))
    if (is.null(msg)) TRUE
    else msg
})

#' @rdname hidden_aliases
#'
#' @importFrom tools file_path_sans_ext
setMethod("backendInitialize", "MsBackendHdf5Peaks",
          function(object, files = character(), data = DataFrame(),
                   hdf5path = character(), ..., BPPARAM = bpparam()) {
              if (!is(data, "DataFrame"))
                  stop("'data' is supposed to be a 'DataFrame' with ",
                       "spectrum data")
              if (!nrow(data))
                  return(object)
              if (length(files) != 1) {
                  if (all(colnames(data) != "dataStorage"))
                      stop("Column \"dataStorage\" is required in 'data'",
                           " if 'files' is missing or length > 1.")
                  if (!length(files))
                      files <- unique(paste0(
                          file_path_sans_ext(data$dataStorage), ".h5"))
              } else if (all(colnames(data) != "dataStorage"))
                  data$dataStorage <- files
              if (length(files) != length(unique(data$dataStorage)))
                  stop("Number of provided file names has to match unique ",
                       "elements in 'data' column \"dataStorage\" (",
                       length(unique(data$dataStorage)), "\"")
              if (length(hdf5path)) {
                  suppressWarnings(hdf5path <- normalizePath(hdf5path))
                  if (!dir.exists(hdf5path))
                      dir.create(hdf5path, recursive = TRUE)
                  files <- file.path(hdf5path, basename(files))
              }
              if (any(file.exists(files)))
                  stop("File(s) ", files[file.exists(files)],
                       " does/do already exist")
              data_storage_levels <- unique(data$dataStorage)
              file_idx <- match(data$dataStorage, data_storage_levels)
              data$dataStorage <- files[file_idx]
              if (!any(colnames(data) == "scanIndex"))
                  data$scanIndex <- seq_len(nrow(data))
              if (any(colnames(data) == "mz")) {
                  if (is.null(data$intensity))
                      data$intensity <- NA
                  peaks <- mapply(cbind,
                                  mz=data$mz,
                                  intensity=data$intensity,
                                  SIMPLIFY=FALSE, USE.NAMES=FALSE)
              } else {
                  mt <- matrix(ncol = 2, nrow = 0,
                               dimnames = list(character(),
                                               c("mz", "intensity")))
                  peaks <- replicate(nrow(data), mt)
              }
              data$mz <- NULL
              data$intensity <- NULL
              file_idx <- factor(file_idx)
              res <- bpmapply(FUN = function(pks, sidx, h5file) {
                  .initialize_h5peaks_file(h5file, modCount = 0L)
                  .h5_write_peaks(pks, scanIndex = sidx, h5file = h5file,
                                  modCount = 0L)
              },
              split(peaks, file_idx),
              split(data$scanIndex, file_idx),
              files,
              BPPARAM = BPPARAM)
              object@modCount <- rep(0L, length(files))
              object@spectraData <- data
              validObject(object)
              object
          })

#' @rdname hidden_aliases
setMethod("show", "MsBackendHdf5Peaks", function(object) {
    callNextMethod()
    fls <- unique(object@spectraData$dataStorage)
    if (length(fls)) {
        to <- min(3, length(fls))
        cat("\nfile(s):\n ", paste(basename(fls[1:to]), collapse = "\n "),
            "\n", sep = "")
        if (length(fls) > 3)
            cat(" ...", length(fls) - 3, "more files\n")
    }
})

#' @rdname hidden_aliases
setMethod("as.list", "MsBackendHdf5Peaks", function(x) {
    if (!length(x))
        return(list())
    fls <- unique(x@spectraData$dataStorage)
    if (length(fls) > 1) {
        f <- factor(dataStorage(x), levels = fls)
        unsplit(bpmapply(
            FUN = .h5_read_peaks,
            fls,
            split(scanIndex(x), f),
            x@modCount,
            SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = bpparam()),
            f)
    } else
        .h5_read_peaks(fls, scanIndex(x), x@modCount)
})

#' @rdname hidden_aliases
setMethod("intensity", "MsBackendHdf5Peaks", function(object) {
    NumericList(lapply(as.list(object), "[", , 2), compress = FALSE)
})

#' @rdname hidden_aliases
setReplaceMethod("intensity", "MsBackendHdf5Peaks", function(object, value) {
    if (!(is.list(value) || inherits(value, "NumericList")))
        stop("'value' has to be a list or NumericList")
    if (length(value) != length(object))
        stop("length of 'value' has to match the length of 'object'")
    mzs <- mz(object)
    if (!all(lengths(value) == lengths(mzs)))
        stop("lengths of 'value' has to match the number of peaks ",
             "(i.e. lengths(object))")
    pks <- mapply(cbind, mz=mzs, intensity=value,
                  SIMPLIFY = FALSE, USE.NAMES = FALSE)
    replaceList(object) <- pks
    object
})

#' @rdname hidden_aliases
setMethod("ionCount", "MsBackendHdf5Peaks", function(object) {
    vapply1d(as.list(object), function(z) sum(z[, 2], na.rm = TRUE))
})

#' @rdname hidden_aliases
setMethod("isCentroided", "MsBackendHdf5Peaks", function(object, ...) {
    vapply1l(as.list(object), .peaks_is_centroided)
})

#' @rdname hidden_aliases
setMethod("isEmpty", "MsBackendHdf5Peaks", function(x) {
    lengths(x) == 0
})

#' @rdname hidden_aliases
setMethod("lengths", "MsBackendHdf5Peaks", function(x, use.names = FALSE) {
    as.integer(lengths(as.list(x)) / 2L)
})

#' @rdname hidden_aliases
setMethod("mz", "MsBackendHdf5Peaks", function(object) {
    NumericList(lapply(as.list(object), "[", , 1), compress = FALSE)
})

#' @rdname hidden_aliases
setReplaceMethod("mz", "MsBackendHdf5Peaks", function(object, value) {
    if (!(is.list(value) || inherits(value, "NumericList")))
        stop("'value' has to be a list or NumericList")
    if (length(value) != length(object))
        stop("length of 'value' has to match the length of 'object'")
    ints <- intensity(object)
    if (!all(lengths(value) == lengths(ints)))
        stop("lengths of 'value' has to match the number of peaks ",
             "(i.e. lengths(object))")
    pks <- mapply(cbind, mz=value, intensity=ints,
                  SIMPLIFY = FALSE, USE.NAMES = FALSE)
    replaceList(object) <- pks
    object
})

#' @rdname hidden_aliases
setReplaceMethod("replaceList", "MsBackendHdf5Peaks", function(object, value) {
    if (length(value) != length(object))
        stop("Length of 'value' has to match length of 'object'")
    if (!(is.list(value) || inherits(value, "SimpleList")))
        stop("'value' has to be a list-like object")
    object@modCount <- object@modCount + 1L
    fls <- unique(object@spectraData$dataStorage)
    if (length(fls)) {
        f <- factor(dataStorage(object), levels = fls)
        res <- bpmapply(FUN = function(pks, sidx, h5file, modC) {
            .h5_write_peaks(pks, scanIndex = sidx, h5file = h5file,
                            modCount = modC)
        },
        split(value, f),
        split(scanIndex(object), f),
        fls,
        object@modCount,
        BPPARAM = bpparam())
    } else
        .h5_write_peaks(value, scanIndex = scanIndex(object), h5file = fls,
                        modCount = object@modCount)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("asDataFrame", "MsBackendHdf5Peaks",
          function(object, columns = spectraVariables(object)) {
              .spectra_data_mzR(object, columns)
          })

#' @rdname hidden_aliases
setReplaceMethod("setDataFrame", "MsBackendHdf5Peaks", function(object, value) {
    pks <- NULL
    if (!inherits(value, "DataFrame"))
        stop("'value' has to be a 'DataFrame'")
    if (nrow(value) != length(object))
        stop("Number of rows of 'value' have to match the length of 'object'")
    if (all(colnames(value) != "dataStorage"))
        value$dataStorage <- object@spectraData$dataStorage
    if (all(colnames(value) != "scanIndex"))
        if (any(colnames(object@spectraData) == "scanIndex"))
            value$scanIndex <- object@spectraData$scanIndex
        else value$scanIndex <- seq_len(nrow(value))
    any_mz <- any(colnames(value) == "mz")
    any_int <- any(colnames(value) == "intensity")
    if (!any_mz && any_int)
        stop("Column \"mz\" required if columns \"intensity\" present")
    if (any_mz) {
        if (!any_int)
            value$intensity <- NA_real_
        pks <- mapply(cbind, mz=value$mz, intensity=value$intensity,
                      SIMPLIFY = FALSE, USE.NAMES = FALSE)
        value$mz <- NULL
        value$intensity <- NULL
    }
    object <- callNextMethod(object, value = value)
    if (length(pks))
        replaceList(object) <- pks
    object
})

#' @rdname hidden_aliases
setReplaceMethod("$", "MsBackendHdf5Peaks", function(x, name, value) {
    if (name == "mz")
        mz(x) <- value
    else if (name == "intensity")
        intensity(x) <- value
    else x <- callNextMethod()
    x
})

#' @rdname hidden_aliases
#'
#' @aliases [,MsBackendHdf5Peaks-method
setMethod("[", "MsBackendHdf5Peaks", function(x, i, j, ..., drop = FALSE) {
    fls <- unique(x@spectraData$dataStorage)
    x <- .subset_backend_data_frame(x, i)
    slot(x, "modCount", check = FALSE) <-
        x@modCount[match(unique(x@spectraData$dataStorage), fls)]
    x
})

#' @rdname hidden_aliases
setMethod("backendMerge", "MsBackendHdf5Peaks", function(object, ...) {
    object <- unname(c(object, ...))
    fls <- lapply(object, function(z) unique(z@spectraData$dataStorage))
    if (anyDuplicated(unlist(fls, use.names = FALSE)))
        stop("Combining backends with the same 'dataStorage' is not supported")
    res <- .combine_backend_data_frame(object)
    res@modCount <- unlist(lapply(object, function(z) z@modCount),
                           use.names = FALSE)
    validObject(res)
    res
})
