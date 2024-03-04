#' @include hidden_aliases.R
NULL

setGeneric("backendBpparam", def = function(object, ...)
    standardGeneric("backendBpparam"))
#' @rdname hidden_aliases
setGeneric("backendParallelFactor", def = function(object, ...)
    standardGeneric("backendParallelFactor"))
#' @rdname hidden_aliases
setMethod("bin", "numeric", MsCoreUtils::bin)
setGeneric("combinePeaks", function(object, ...)
           standardGeneric("combinePeaks"))
#' @rdname hidden_aliases
setGeneric("containsMz", function(object, ...)
    standardGeneric("containsMz"))
#' @rdname hidden_aliases
setGeneric("containsNeutralLoss", function(object, ...)
    standardGeneric("containsNeutralLoss"))
#' @rdname hidden_aliases
setGeneric("dropNaSpectraVariables", function(object, ...)
    standardGeneric("dropNaSpectraVariables"))
#' @rdname hidden_aliases
setGeneric("entropy", function(object, ...)
  standardGeneric("entropy"))
#' @rdname hidden_aliases
setGeneric("export", function(object, ...)
    standardGeneric("export"))
setGeneric("filterFourierTransformArtefacts", function(object, ...)
    standardGeneric("filterFourierTransformArtefacts"))
#' @rdname hidden_aliases
setGeneric("filterMzRange", function(object, ...)
    standardGeneric("filterMzRange"))
#' @rdname hidden_aliases
setGeneric("filterMzValues", function(object, ...)
    standardGeneric("filterMzValues"))
#' @rdname hidden_aliases
setGeneric("filterPrecursorMzValues", function(object, ...)
    standardGeneric("filterPrecursorMzValues"))
#' @rdname hidden_aliases
setGeneric("filterPrecursorMzRange", function(object, ...)
    standardGeneric("filterPrecursorMzRange"))
#' @rdname hidden_aliases
setGeneric("filterRanges", function(object, ...)
    standardGeneric("filterRanges"))
#' @rdname hidden_aliases
setGeneric("filterValues", function(object, ...)
    standardGeneric("filterValues"))
#' @rdname hidden_aliases
setGeneric("isReadOnly", function(object, ...)
    standardGeneric("isReadOnly"))
#' @rdname neutralLoss
setGeneric("neutralLoss", function(object, param, ...)
    standardGeneric("neutralLoss"))
#' @rdname hidden_aliases
setGeneric("pickPeaks", function(object, ...)
    standardGeneric("pickPeaks"))
setGeneric("plotSpectraMirror", function(x, y, ...)
    standardGeneric("plotSpectraMirror"))
#' @rdname hidden_aliases
setGeneric("replaceIntensitiesBelow", function(object, threshold = min, ...)
    standardGeneric("replaceIntensitiesBelow"))
#' @rdname hidden_aliases
setGeneric("reset", function(object, ...)
    standardGeneric("reset"))
#' @rdname hidden_aliases
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

setGeneric("supportsSetBackend", function(object, ...)
    standardGeneric("supportsSetBackend"))
