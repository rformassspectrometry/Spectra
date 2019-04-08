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
#' @name MsBackend
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @md
#'
#' @exportClass MsBackend
NULL

setClass("MsBackend",
         slots = c(spectraData = "DataFrame",
                   files = "character",
                   modCount = "integer",
                   version = "character"),
         prototype = prototype(spectraData = DataFrame(),
                               files = character(),
                               modCount = integer(),
                               version = "0.1"))

setValidity("MsBackend", function(object) {
    msg <- .valid_spectra_data_required_columns(object@spectraData)
    if (length(msg))
        return(msg)
    msg <- c(
        .valid_column_datatype(object@spectraData, .SPECTRA_DATA_COLUMNS),
        .valid_intensity_column(object@spectraData),
        .valid_mz_column(object@spectraData),
        .valid_ms_backend_files(object@files),
        .valid_ms_backend_files_from_file(object@files,
                                          object@spectraData$fromFile),
        .valid_ms_backend_mod_count(object@files, object@modCount))
    if (is.null(msg)) TRUE
    else msg
})

#' Initialize a backend
#'
#' This generic is used to setup and initialize a backend.
#'
#' It should be reimplemented if the backend needs specific requirements,
#' e.g. for *HDF5* the creation of .h5 files; for *SQL* the creation of a
#' database or tables or establish and check a database connection.
#'
#' @param object An object inheriting from [Backend-class],
#' i.e. [BackendHdf5-class]
#'
#' @param files The path to the source (generally .mzML) files.
#'
#' @param spectraData A [S4Vectors::DataFrame-class]
#'
#' @param ... Other arguments passed to the methods.
#'
#' @return A fully operational [Backend-class] derivate.
#'
#' @family Backend generics
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#'
#' @noRd
setGeneric(
    "backendInitialize",
    def = function(object, files, spectraData, ...)
        standardGeneric("backendInitialize"),
    valueClass = "MsBackend"
)
setMethod("backendInitialize", signature = "MsBackend",
          definition = function(object, files, spectraData, ...) {
              if (missing(files)) files <- character()
              if (missing(spectraData)) spectraData <- DataFrame()
              object@files <- files
              object@modCount <- integer(length(files))
              object@spectraData <- spectraData
              validObject(object)
              object
          })

#' @exportMethod length
#'
#' @noRd
setMethod("length", "MsBackend", function(x) {
    nrow(x@spectraData)
})

#' @exportMethod msLevel
#'
#' @noRd
setMethod("msLevel", "MsBackend", function(object, ...) {
    if (any(colnames(object@spectraData) == "msLevel"))
        object@spectraData$msLevel
    else rep(NA_integer_, times = length(object))
})
