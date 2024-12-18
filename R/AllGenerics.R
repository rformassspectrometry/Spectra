#' @include hidden_aliases.R
NULL

setGeneric("backendRequiredSpectraVariables", function(object, ...)
          standardGeneric("backendRequiredSpectraVariables"))
#' @rdname hidden_aliases
setMethod("bin", "numeric", MsCoreUtils::bin)
setGeneric("combinePeaks", function(object, ...)
           standardGeneric("combinePeaks"))
setGeneric("containsMz", function(object, ...)
    standardGeneric("containsMz"))
setGeneric("containsNeutralLoss", function(object, ...)
    standardGeneric("containsNeutralLoss"))
setGeneric("dataStorageBasePath", function(object, ...)
    standardGeneric("dataStorageBasePath"))
setGeneric("dataStorageBasePath<-", function(object, ..., value)
    standardGeneric("dataStorageBasePath<-"))
setGeneric("dropNaSpectraVariables", function(object, ...)
    standardGeneric("dropNaSpectraVariables"))
setGeneric("entropy", function(object, ...)
  standardGeneric("entropy"))
setGeneric("export", function(object, ...)
    standardGeneric("export"))
setGeneric("filterFourierTransformArtefacts", function(object, ...)
    standardGeneric("filterFourierTransformArtefacts"))
setGeneric("neutralLoss", function(object, param, ...)
    standardGeneric("neutralLoss"))
setGeneric("pickPeaks", function(object, ...)
    standardGeneric("pickPeaks"))
setGeneric("plotSpectraMirror", function(x, y, ...)
    standardGeneric("plotSpectraMirror"))
setGeneric("replaceIntensitiesBelow", function(object, threshold = min, ...)
    standardGeneric("replaceIntensitiesBelow"))
setGeneric("reset", function(object, ...)
    standardGeneric("reset"))
setGeneric("selectSpectraVariables", function(object, ...)
    standardGeneric("selectSpectraVariables"))
setGeneric("Spectra", function(object, ...) standardGeneric("Spectra"))

#' @title Mapping between spectra variables and data file fields
#'
#' @description
#'
#' The `spectraVariableMapping` function provides the mapping
#' between *spectra variables* of a [Spectra()] object with data fields from a
#' data file. Such name mapping is expected to enable an easier import of data
#' files with specific *dialects*, e.g. files in MGF format that use a
#' different naming convention for core spectra variables.
#'
#' [MsBackend()] implementations are expected to implement this function
#' (if needed) to enable import of data from file formats with non-standardized
#' data fields.
#'
#' @param object An instance of an object extending [MsBackend()].
#'
#' @param ... Optional parameters.
#'
#' @return A named `character` with names being spectra variable names (use
#'     [spectraVariables()] for a list of supported names) and values being the
#'     data field names.
#'
#' @author Johannes Rainer
#'
#' @export
setGeneric("spectraVariableMapping", function(object, ...)
    standardGeneric("spectraVariableMapping"))
