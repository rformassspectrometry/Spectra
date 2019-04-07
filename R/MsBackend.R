#' @include hidden_aliases.R
NULL

#' @title Mass spectrometry data backends
#'
#' @description
#'
#' [MsBackend-class] objects provide access to mass spectrometry data. Such
#' backends can be generally classidied into *in-memory* and *on-disk* backends.
#' In-memory backends, such as the base `MsBackend`, keep all (spectra) data in
#' memory ensuring fast data access. On-disk backends keep only part of the
#' data in memory retrieving the remaining data (mostly m/z and intensity
#' values) on-demand from disk.
#'
#' Available backends for [Spectra()] objects are listed below. The last section
#' provides information for the development of new backends.
#'
#' @section `MsBackend`:
#'
#' The base backend class, `MsBackend`, keeps all spectra data in a `DataFrame`.
#'
#' @section Implementation notes for new backend classes:
#'
#' New backend classes should extend the base `MsBackend` class.
#'
#' TODO: add description of required methods etc.
#'
#' @name Backend
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @md
#'
#' @exportClass MsBackend
NULL

setClass("MsBackend",
         slots = c(spectraData = "DataFrame",
                   version = "character"),
         prototype = prototype(spectraData = DataFrame(),
                               version = "0.1"))

setValidity("MsBackend", function(object) {
    msg <- c(.valid_column_datatype(object@spectraData, .SPECTRA_DATA_COLUMNS),
             .valid_intensity_column(object@spectraData),
             .valid_mz_column(object@spectraData))
    if (is.null(msg)) TRUE
    else msg
})
