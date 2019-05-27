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
         slots = c(h5files = "character"),
         prototype = prototype(version = "0.1", readonly = FALSE,
                               h4files = character()))

setValidity("MsBackendHdf5Peaks", function(object) {
    msg <- .valid_h5files(object@h5files, object@files)
    if (is.null(msg)) TRUE
    else msg
})

setMethod("backendInitialize", "MsBackendHdf5Peaks",
          function(object, files, spectraData, hdf5path = ".", ...,
                   BPPARAM = bpparam()) {
              ## create hdf5 file names.
              ## - if spectraData contains m/z and intensity -> save to hdf5 and
              ##   remove from spectraData
          })
