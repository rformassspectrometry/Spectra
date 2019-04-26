#' @include hidden_aliases.R
NULL

#' @rdname hidden_aliases
setGeneric("acquisitionNum", function(object, ...)
    standardGeneric("acquisitionNum"))
#' @rdname hidden_aliases
setGeneric("backendInitialize",
           def = function(object, files, spectraData, ...)
               standardGeneric("backendInitialize"),
           valueClass = "MsBackend"
)
#' @rdname hidden_aliases
setGeneric("centroided", function(object, ...)
    standardGeneric("centroided"))
#' @rdname hidden_aliases
setGeneric("centroided<-", function(object, value)
    standardGeneric("centroided<-"))
#' @rdname hidden_aliases
setGeneric("clean", function(object, ...) standardGeneric("clean"))
#' @rdname hidden_aliases
setGeneric("collisionEnergy", function(object, ...)
    standardGeneric("collisionEnergy"))
#' @rdname hidden_aliases
setGeneric("collisionEnergy<-", function(object, value)
    standardGeneric("collisionEnergy<-"))
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
setGeneric("msLevel", function(object, ...)
    standardGeneric("msLevel"))
#' @rdname hidden_aliases
setGeneric("peaksCount", function(object, ...)
    standardGeneric("peaksCount"))
#' @rdname hidden_aliases
setGeneric("polarity", function(object, ...)
    standardGeneric("polarity"))
#' @rdname hidden_aliases
setGeneric("polarity<-", function(object, value)
    standardGeneric("polarity<-"))
#' @rdname hidden_aliases
setGeneric("precScanNum", function(object, ...)
    standardGeneric("precScanNum"))
#' @rdname hidden_aliases
setGeneric("precursorCharge", function(object, ...)
    standardGeneric("precursorCharge"))
#' @rdname hidden_aliases
setGeneric("precursorIntensity", function(object, ...)
    standardGeneric("precursorIntensity"))
#' @rdname hidden_aliases
setGeneric("precursorMz", function(object, ...)
    standardGeneric("precursorMz"))
#' @rdname hidden_aliases
setGeneric("removePeaks", function(object, t="min", ...)
    standardGeneric("removePeaks"))
#' @rdname hidden_aliases
setGeneric("rtime<-", function(object, value)
    standardGeneric("rtime<-"))
#' @rdname hidden_aliases
setGeneric("scanIndex", function(object, ...)
    standardGeneric("scanIndex"))
#' @rdname hidden_aliases
setGeneric("selectSpectraVariables", function(object, ...)
    standardGeneric("selectSpectraVariables"))
#' @rdname hidden_aliases
setGeneric("smoothed", function(object, ...)
    standardGeneric("smoothed"))
#' @rdname hidden_aliases
setGeneric("smoothed<-", function(object, value)
    standardGeneric("smoothed<-"))
#' @rdname hidden_aliases
setGeneric("Spectra", function(object, ...) standardGeneric("Spectra"))
#' @rdname hidden_aliases
setGeneric("spectraData", function(object, ...)
    standardGeneric("spectraData"))
#' @rdname hidden_aliases
setGeneric("spectraData<-", function(object, value)
    standardGeneric("spectraData<-"))
#' @rdname hidden_aliases
setGeneric("spectraNames", function(object, ...)
    standardGeneric("spectraNames"))
#' @rdname hidden_aliases
setGeneric("spectraNames<-", function(object, value)
    standardGeneric("spectraNames<-"))
#' @rdname hidden_aliases
setGeneric("spectraVariables", function(object, ...)
    standardGeneric("spectraVariables"))
