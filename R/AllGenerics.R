#' @include hidden_aliases.R
NULL

#' @rdname hidden_aliases
setGeneric("replaceList<-", function(object, value)
    standardGeneric("replaceList<-"))
#' @rdname hidden_aliases
setGeneric("backendInitialize", def = function(object, ...)
    standardGeneric("backendInitialize"),
    valueClass = "MsBackend")
#' @rdname hidden_aliases
setGeneric("backendMerge", def = function(object, ...)
    standardGeneric("backendMerge"),
    valueClass = "MsBackend")
#' @rdname hidden_aliases
setGeneric("bin", function(x, ...) standardGeneric("bin"))
#' @rdname hidden_aliases
setMethod("bin", "numeric", MsCoreUtils::bin)
#' @rdname hidden_aliases
setGeneric("compareSpectra", function(x, y, ...)
    standardGeneric("compareSpectra"))
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
setGeneric("export", function(object, ...)
    standardGeneric("export"))
#' @rdname hidden_aliases
setGeneric("filterIntensity", function(object, ...)
    standardGeneric("filterIntensity"))
#' @rdname hidden_aliases
setGeneric("isReadOnly", function(object, ...)
    standardGeneric("isReadOnly"))
#' @rdname hidden_aliases
setGeneric("peaksData", function(object, ...) standardGeneric("peaksData"))
#' @rdname hidden_aliases
setGeneric("peaksData<-", function(object, value)
    standardGeneric("peaksData<-"))
#' @rdname hidden_aliases
setGeneric("pickPeaks", function(object, ...)
    standardGeneric("pickPeaks"))
#' @rdname hidden_aliases
setGeneric("replaceIntensitiesBelow", function(object, threshold = min, ...)
    standardGeneric("replaceIntensitiesBelow"))
#' @rdname hidden_aliases
setGeneric("reset", function(object, ...)
    standardGeneric("reset"))
#' @rdname hidden_aliases
setGeneric("selectSpectraVariables", function(object, ...)
    standardGeneric("selectSpectraVariables"))
#' @rdname hidden_aliases
setGeneric("setBackend", function(object, backend, ...)
    standardGeneric("setBackend"))
setGeneric("Spectra", function(object, ...) standardGeneric("Spectra"))
