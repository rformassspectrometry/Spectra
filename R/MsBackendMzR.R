#' @include hidden_aliases.R
NULL

#' @title mzR-based backend
#'
#' @description
#'
#' The `MsBackendMzR` inherits all slots and methods from the base `MsBackend`
#' (in-memory) backend. It overrides the base `mz` and `intensity` methods as
#' well as `peaks` to read the respective data from the original raw data
#' files.
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
         contains = "MsBackend",
         prototype = prototype(version = "0.1"))

setValidity("MsBackendMzR", function(object) {
    cat("MsBackendMzR validity\n")
    msg <- .valid_ms_backend_files_exist(object@files)
})

#' @description
#'
#' `backendInitialize` reads the header information from all raw data files
#' and stores this information to the object's `spectraData` slot.
#'
#' @importFrom BiocParallel bpparam bpmapply
#'
#' @noRd
setMethod("backendInitialize", "MsBackendMzR",
          function(object, files, spectraData, ..., BPPARAM = bpparam()) {
              if (missing(files) || !length(files))
                  stop("Parameter 'files' is mandatory for 'MsBackendMzR'")
              files <- normalizePath(files)
              msg <- .valid_ms_backend_files_exist(files)
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

## backendInitialize: read header from file and add/put that into spectraData
## Need also a column spIdx, which is the index of the spectrum within the file.
## scanIndex corresponds to seqNum which is reported by the header call.
## Seems to be the index of the spectrum/scan within the file (also by looking
## at the pwiz code. Subsetting needs the `scanIndex` column.

## Need to implement the mz, intensity and peaks methods that will read the
## data from the original files