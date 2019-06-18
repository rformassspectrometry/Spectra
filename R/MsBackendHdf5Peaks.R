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
    fls <- dataStorageLevels(object)
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
          function(object, files = character(), spectraData = DataFrame(),
                   hdf5path = character(), ..., BPPARAM = bpparam()) {
              if (!is(spectraData, "DataFrame"))
                  stop("'spectraData' is supposed to be a 'DataFrame' with ",
                       "spectrum data")
              if (!nrow(spectraData))
                  return(object)
              if (length(files) != 1) {
                  if (all(colnames(spectraData) != "dataStorage"))
                      stop("Column \"dataStorage\" is required in 'spectraData'",
                           " if 'files' is missing or length > 1.")
                  if (!length(files))
                      files <- unique(paste0(
                          file_path_sans_ext(spectraData$dataStorage), ".h5"))
              } else if (all(colnames(spectraData) != "dataStorage"))
                  spectraData$dataStorage <- files
              if (length(files) != length(unique(spectraData$dataStorage)))
                  stop("Number of provided file names has to match unique ",
                       "elements in 'spectraData' column \"dataStorage\" (",
                       length(unique(spectraData$dataStorage)), "\"")
              if (length(hdf5path)) {
                  suppressWarnings(hdf5path <- normalizePath(hdf5path))
                  if (!dir.exists(hdf5path))
                      dir.create(hdf5path, recursive = TRUE)
                  files <- file.path(hdf5path, basename(files))
              }
              if (any(file.exists(files)))
                  stop("File(s) ", files[file.exists(files)],
                       " does/do already exist")
              data_storage_levels <- unique(spectraData$dataStorage)
              file_idx <- match(spectraData$dataStorage, data_storage_levels)
              spectraData$dataStorage <- files[file_idx]
              if (!any(colnames(spectraData) == "scanIndex"))
                  spectraData$scanIndex <- seq_len(nrow(spectraData))
              if (any(colnames(spectraData) == "mz")) {
                  if (is.null(spectraData$intensity))
                      spectraData$intensity <- NA
                  peaks <- mapply(spectraData$mz, spectraData$intensity,
                                  FUN = function(mz, intensity) {
                                      cbind(mz, intensity)
                                  })
              } else {
                  mt <- matrix(ncol = 2, nrow = 0,
                               dimnames = list(character(),
                                               c("mz", "intensity")))
                  peaks <- replicate(nrow(spectraData), mt)
              }
              spectraData$mz <- NULL
              spectraData$intensity <- NULL
              file_idx <- factor(file_idx)
              res <- bpmapply(FUN = function(pks, sidx, h5file) {
                  .initialize_h5peaks_file(h5file, modCount = 0L)
                  .h5_write_peaks(pks, scanIndex = sidx, h5file = h5file,
                                  modCount = 0L)
              },
              split(peaks, file_idx),
              split(spectraData$scanIndex, file_idx),
              files,
              BPPARAM = BPPARAM)
              object@modCount <- rep(0L, length(files))
              object@spectraData <- .as_rle_spectra_data(spectraData)
              validObject(object)
              object
          })

#' @rdname hidden_aliases
setMethod("show", "MsBackendHdf5Peaks", function(object) {
    callNextMethod()
    fls <- dataStorageLevels(object)
    if (length(fls)) {
        to <- min(3, length(fls))
        cat("\nfile(s):\n ", paste(basename(fls[1:to]), collapse = "\n "),
            "\n", sep = "")
        if (length(fls) > 3)
            cat(" ...", length(fls) - 3, "more files\n")
    }
})

#' @rdname hidden_aliases
setMethod("intensity", "MsBackendHdf5Peaks", function(object) {
    NumericList(lapply(peaks(object), "[", , 2), compress = FALSE)
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
             "(i.e. peaksCount(object))")
    pks <- mapply(cbind, mz=mzs, intensity=value)
    peaks(object) <- pks
    object
})

#' @rdname hidden_aliases
setMethod("ionCount", "MsBackendHdf5Peaks", function(object) {
    vapply(peaks(object), function(z) sum(z[, 2], na.rm = TRUE), numeric(1))
})

#' @rdname hidden_aliases
setMethod("isCentroided", "MsBackendHdf5Peaks", function(object, ...) {
    vapply(peaks(object), .isCentroided, logical(1))
})

#' @rdname hidden_aliases
setMethod("isEmpty", "MsBackendHdf5Peaks", function(x) {
    peaksCount(x) == 0
})

#' @rdname hidden_aliases
setMethod("mz", "MsBackendHdf5Peaks", function(object) {
    NumericList(lapply(peaks(object), "[", , 1), compress = FALSE)
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
             "(i.e. peaksCount(object))")
    pks <- mapply(cbind, mz=value, intensity=ints)
    peaks(object) <- pks
    object
})

#' @rdname hidden_aliases
setMethod("peaks", "MsBackendHdf5Peaks", function(object) {
    if (!length(object))
        return(list())
    fls <- dataStorageLevels(object)
    if (length(fls) > 1) {
        f <- factor(dataStorage(object), levels = fls)
        unsplit(bpmapply(
            FUN = .h5_read_peaks,
            fls,
            split(scanIndex(object), f),
            object@modCount,
            SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = bpparam()),
            f)
    } else
        .h5_read_peaks(fls, scanIndex(object), object@modCount)
})

#' @rdname hidden_aliases
setReplaceMethod("peaks", "MsBackendHdf5Peaks", function(object, value) {
    if (length(value) != length(object))
        stop("Length of 'value' has to match length of 'object'")
    if (!(is.list(value) || inherits(value, "SimpleList")))
        stop("'value' has to be a list-like object")
    object@modCount <- object@modCount + 1L
    fls <- dataStorageLevels(object)
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
setMethod("peaksCount", "MsBackendHdf5Peaks", function(object) {
    as.integer(lengths(peaks(object)) / 2L)
})

#' @rdname hidden_aliases
setMethod("spectraData", "MsBackendHdf5Peaks",
          function(object, columns = spectraVariables(object)) {
              .spectra_data_mzR(object, columns)
          })

#' @rdname hidden_aliases
setReplaceMethod("spectraData", "MsBackendHdf5Peaks", function(object, value) {
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
        pks <- mapply(function(mz, intensity) cbind(mz, intensity),
                      value$mz, value$intensity, SIMPLIFY = FALSE,
                      USE.NAMES = FALSE)
        value$mz <- NULL
        value$intensity <- NULL
    }
    object <- callNextMethod(object, value = value)
    if (length(pks))
        peaks(object) <- pks
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
setMethod("[", "MsBackendHdf5Peaks", function(x, i, j, ..., drop = FALSE) {
    fls <- dataStorageLevels(x)
    if (!missing(j))
        stop("Subsetting by column ('j = ", j, "' is not supported")
    i <- .i_to_index(i, length(x), rownames(x@spectraData))
    x@spectraData <- x@spectraData[i, , drop = FALSE]
    x@modCount <- x@modCount[match(dataStorageLevels(x), fls)]
    validObject(x)
    x
})

#' @rdname hidden_aliases
setMethod("backendMerge", "MsBackendHdf5Peaks", function(object, ...) {
    object <- unname(c(object, ...))
    fls <- lapply(object, dataStorageLevels)
    if (any(duplicated(unlist(fls, use.names = FALSE))))
        stop("Combining backends with the same 'dataStorage' is not supported")
    res <- .combine_backend_data_frame(object)
    res@modCount <- unlist(lapply(object, function(z) z@modCount),
                           use.names = FALSE)
    validObject(res)
    res
})
