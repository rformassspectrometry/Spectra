#' @include hidden_aliases.R
NULL

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
setGeneric("dataOrigin", function(object, ...) standardGeneric("dataOrigin"))
#' @rdname hidden_aliases
setGeneric("dataOrigin<-", function(object, value)
    standardGeneric("dataOrigin<-"))
#' @rdname hidden_aliases
setGeneric("dataStorage", function(object, ...) standardGeneric("dataStorage"))
#' @rdname hidden_aliases
setGeneric("dataStorage<-", function(object, value)
    standardGeneric("dataStorage<-"))
#' @rdname hidden_aliases
setGeneric("filterAcquisitionNum", function(object, ...)
    standardGeneric("filterAcquisitionNum"))
#' @rdname hidden_aliases
setGeneric("filterDataOrigin", function(object, ...)
    standardGeneric("filterDataOrigin"))
#' @rdname hidden_aliases
setGeneric("filterDataStorage", function(object, ...)
    standardGeneric("filterDataStorage"))
#' @rdname hidden_aliases
setGeneric("filterEmptySpectra", function(object, ...)
    standardGeneric("filterEmptySpectra"))
#' @rdname hidden_aliases
setGeneric("filterIsolationWindow", function(object, ...)
    standardGeneric("filterIsolationWindow"))
#' @rdname hidden_aliases
setGeneric("filterMsLevel", function(object, ...)
    standardGeneric("filterMsLevel"))
#' @rdname hidden_aliases
setGeneric("filterPolarity", function(object, ...)
    standardGeneric("filterPolarity"))
#' @rdname hidden_aliases
setGeneric("filterPrecursorMz", function(object, ...)
    standardGeneric("filterPrecursorMz"))
#' @rdname hidden_aliases
setGeneric("filterPrecursorScan", function(object, ...)
    standardGeneric("filterPrecursorScan"))
#' @rdname hidden_aliases
setGeneric("filterRt", function(object, ...) standardGeneric("filterRt"))
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
#' @rdname hidden_aliases
setGeneric("setBackend", function(object, backend, ...)
    standardGeneric("setBackend"))
setGeneric("Spectra", function(object, ...) standardGeneric("Spectra"))
