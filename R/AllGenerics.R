#' @include hidden_aliases.R
NULL

#' @rdname hidden_aliases
setGeneric("asDataFrame", function(object, ...) standardGeneric("asDataFrame"))
#' @rdname hidden_aliases
setGeneric("asDataFrame<-", function(object, value)
    standardGeneric("asDataFrame<-"))
#' @rdname hidden_aliases
setGeneric("backendInitialize",
           def = function(object, ...)
               standardGeneric("backendInitialize"),
           valueClass = "MsBackend"
           )
#' @rdname hidden_aliases
setGeneric("backendMerge", def = function(object, ...)
    standardGeneric("backendMerge"),
    valueClass = "MsBackend")
#' @rdname hidden_aliases
setGeneric("bin", function(x, ...) standardGeneric("bin"))
#' @rdname hidden_aliases
setMethod("bin", "numeric", MsCoreUtils::bin)
#' @rdname hidden_aliases
setGeneric("clean", function(object, ...) standardGeneric("clean"))
#' @rdname hidden_aliases
setGeneric("compareSpectra", function(x, y, ...)
    standardGeneric("compareSpectra"))
#' @rdname hidden_aliases
setGeneric("dropNaSpectraVariables", function(object, ...)
    standardGeneric("dropNaSpectraVariables"))
#' @rdname hidden_aliases
setGeneric("isReadOnly", function(object, ...)
    standardGeneric("isReadOnly"))
#' @rdname hidden_aliases
setGeneric("pickPeaks", function(object, ...)
    standardGeneric("pickPeaks"))
#' @rdname hidden_aliases
setGeneric("removePeaks", function(object, t="min", ...)
    standardGeneric("removePeaks"))
#' @rdname hidden_aliases
setGeneric("replaceList<-", function(object, value)
    standardGeneric("replaceList<-"))
#' @rdname hidden_aliases
setGeneric("selectSpectraVariables", function(object, ...)
    standardGeneric("selectSpectraVariables"))
#' @rdname hidden_aliases
setGeneric("setBackend", function(object, backend, ...)
    standardGeneric("setBackend"))
setGeneric("Spectra", function(object, ...) standardGeneric("Spectra"))
