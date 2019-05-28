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
         prototype = prototype(version = "0.1", readonly = FALSE,
                               h4files = character()))

setValidity("MsBackendHdf5Peaks", function(object) {
    msg <- .valid_spectra_data_required_columns(object@spectraData,
                                                c("fromFile", "scanIndex"))
    msg <- c(msg, .valid_h5files(object@files))
    if (is.null(msg)) TRUE
    else msg
})

#' @rdname hidden_aliases
#'
#' @importFrom tools file_ext
setMethod("backendInitialize", "MsBackendHdf5Peaks",
          function(object, files = character(), spectraData,
                   hdf5path = character(), modCount = 0L, ...,
                   BPPARAM = bpparam()) {
              if (missing(spectraData))
                  spectraData <- DataFrame()
              else if (!is(spectraData, "DataFrame"))
                  stop("'spectraData' is supposed to be a 'DataFrame' with ",
                       "spectrum data")
              if (!nrow(spectraData))
                  return(object)
              if (!length(files) || is.na(files))
                  files <- "spectra_peaks.h5"
              if (length(hdf5path)) {
                  suppressWarnings(hdf5path <- normalizePath(hdf5path))
                  if (!dir.exists(hdf5path))
                      dir.create(hdf5path, recursive = TRUE)
                  files <- file.path(hdf5path, files)
              }
              if (length(files) != length(modCount))
                  modCount <- rep(modCount[1], length(files))
              if (any(file.exists(files)))
                  stop("File(s) ", files[file.exists(files)],
                       " does/do already exist")
              if (!any(colnames(spectraData) == "fromFile"))
                  spectraData$fromFile <- 1L
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
                  mt <- matrix(ncol = 2, nrow = 0)
                  peaks <- replicate(nrow(spectraData), mt)
              }
              spectraData$mz <- NULL
              spectraData$intensity <- NULL
              if (length(files) != length(unique(spectraData$fromFile)))
                  stop("Number of files does not match with indices in",
                       " 'spectraData' column 'fromFile'")
              fromF <- factor(spectraData$fromFile,
                              levels = unique(spectraData$fromFile))
              res <- bpmapply(FUN = function(pks, sidx, h5file, modC) {
                  .initialize_h5peaks_file(h5file, modCount = modC)
                  .h5_write_peaks(pks, scanIndex = sidx, h5file = h5file,
                                  modCount = modC)
              },
              split(peaks, fromF),
              split(spectraData$scanIndex, fromF),
              files[as.integer(levels(fromF))],
              modCount[as.integer(levels(fromF))],
              BPPARAM = BPPARAM)
              object@modCount <- modCount
              callNextMethod(object = object, files = files,
                             spectraData = spectraData, ...)
          })

#' @rdname hidden_aliases
setMethod("show", "MsBackendHdf5Peaks", function(object) {
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
setMethod("intensity", "MsBackendHdf5Peaks", function(object) {
    SimpleList(lapply(peaks(object), function(z) z[, 2]))
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
    SimpleList(lapply(peaks(object), function(z) z[, 1]))
})

#' @rdname hidden_aliases
setMethod("peaks", "MsBackendHdf5Peaks", function(object) {
    if (!length(object))
        return(list())
    fromF <- factor(object@spectraData$fromFile,
                    levels = unique(object@spectraData$fromFile))
    if (length(object@files) > 1) {
        unsplit(mapply(
            FUN = .h5_read_peaks,
            object@files[as.integer(levels(fromF))],
            split(scanIndex(object), fromF),
            object@modCount[as.integer(levels(fromF))],
            SIMPLIFY = FALSE, USE.NAMES = FALSE),
            fromF)
    } else
        .h5_read_peaks(object@files, object@spectraData$scanIndex,
                       object@modCount)
})

#' @rdname hidden_aliases
setMethod("peaksCount", "MsBackendHdf5Peaks", function(object) {
    vapply(peaks(object), nrow, integer(1))
})

#' @rdname hidden_aliases
setMethod("spectraData", "MsBackendHdf5Peaks",
          function(object, columns = spectraVariables(object)) {
              .spectra_data_mzR(object, columns)
          })

## spectraData<- (support m/z and intensity replacement) -> increase modCount
## $<- (support m/z and intensity replacement) -> increase modCount
