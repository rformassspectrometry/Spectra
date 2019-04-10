#' @include hidden_aliases.R
NULL

#' @title mzR-based backend
#'
#' @description
#'
#' The `MsBackendMzR` inherits all slots and methods from the base
#' `MsBackendMemory` (in-memory) backend. It overrides the base `mz` and
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
         contains = "MsBackendMemory",
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
              ## Read the header information in parallel
              callNextMethod(object = object, files = files,
                             spectraData = spectraData, ...)
          })

## Need to implement the mz, intensity and peaks methods that will read the
## data from the original files