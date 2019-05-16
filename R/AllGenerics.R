#' @include hidden_aliases.R
NULL

#' @rdname hidden_aliases
setGeneric("backendInitialize",
           def = function(object, files, spectraData, ...)
               standardGeneric("backendInitialize"),
           valueClass = "MsBackend"
)
#' @rdname hidden_aliases
setGeneric("clean", function(object, ...) standardGeneric("clean"))
#' @rdname hidden_aliases
setGeneric("fileNames", function(object, ...) standardGeneric("fileNames"))
#' @rdname hidden_aliases
setGeneric("fromFile", function(object, ...)
    standardGeneric("fromFile"))
#' @rdname hidden_aliases
setGeneric("ionCount", function(object, ...)
    standardGeneric("ionCount"))
#' @rdname hidden_aliases
setGeneric("isReadOnly", function(object, ...)
    standardGeneric("isReadOnly"))
#' @rdname hidden_aliases
setGeneric("peaksCount", function(object, ...)
    standardGeneric("peaksCount"))
#' @rdname hidden_aliases
setGeneric("removePeaks", function(object, t="min", ...)
    standardGeneric("removePeaks"))
#' @rdname hidden_aliases
setGeneric("selectSpectraVariables", function(object, ...)
    standardGeneric("selectSpectraVariables"))
setGeneric("Spectra", function(object, ...) standardGeneric("Spectra"))
