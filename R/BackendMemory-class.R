#' @include hidden_aliases.R
NULL

setClass("BackendMemory",
    contains = "Backend",
    slots = c(
        ## list of data.frames with columns mz and intensity.
        spectra = "list"
    )
)

.valid.BackendMemory.spectra.names <- function(x) {
    n <- length(x)

    if (n) {
        nms <- names(x)
        if (any(is.null(nms)))
            return("Spectra names should not be NULL.")
        if (anyNA(nms))
            return("Spectra names should not contain NA.")
        if (!all(nchar(nms)))
            return("Spectra names should not be missing.")
        if (anyDuplicated(nms))
            return("Duplicated spectra names found.")
    }
    NULL
}

.valid_spectrum_data_frame <- function(x) {
    if (is.data.frame(x)) {
        if (!all(c("mz", "intensity") %in% colnames(x)))
            return("Spectra data.frame must have columns \"mz\" and \"intensity\"")
        if (is.unsorted(x$mz))
            return("m/z values have to be increasingly ordered")
        if (!is.numeric(x$mz) || !is.numeric(x$intensity))
            return("columns \"mz\" and \"intensity\" need to be of type numeric")
    } else return("Spectra should be a data.frame")
    NULL
}

setValidity("BackendMemory", function(object) {
    lapply(object@spectra, function(z) {
        res <- .valid_spectrum_data_frame
        if (length(character))
            stop("Validity: first error: ", res, call. = FALSE)
    })
    msg <- .valid.BackendMemory.spectra.names(object@spectra)
    if (is.null(msg)) { TRUE } else { msg }
})

#' @rdname Backend
#'
#' @export
BackendMemory <- function() { new("BackendMemory") }

#' @rdname hidden_aliases
setMethod("backendSubset", "BackendMemory", function(object, spectraData) {
    fidx <- unique(spectraData$fileIdx)
    ## Update also `@fromFile` in the spectra.
    object@spectra <- object@spectra[rownames(spectraData)]
    callNextMethod()
})

#' @rdname hidden_aliases
setReplaceMethod(
    "backendSplitByFile",
    "BackendMemory",
    function(object, spectraData, ..., value) {
    rn <- split(rownames(spectraData), spectraData$fileIdx)
    fidx <- as.integer(sort.int(unique(spectraData$fileIdx)))
    for (i in seq_along(value)) {
        object@spectra[rn[[i]]] <- value[[i]]@spectra
    }
    callNextMethod()
})

#' @rdname hidden_aliases
setMethod(
    "backendInitialize",
    signature = "BackendMemory",
    definition = function(object, files, spectraData, ..., BPPARAM = bpparam()) {

    object@spectra <- vector(mode = "list", length = nrow(spectraData))
    names(object@spectra) <- rownames(spectraData)
    object <- callNextMethod()
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod(
    "backendReadSpectra",
    signature = "BackendMemory",
    definition = function(object, spectraData, ..., BPPARAM = bpparam()) {
    object@spectra[rownames(spectraData)]
})

#' @rdname hidden_aliases
setMethod(
    "backendWriteSpectra",
    signature = "BackendMemory",
    definition = function(object, spectra, spectraData, updateModCount, ...) {
        object@spectra[rownames(spectraData)] <- spectra
        if (updateModCount) {
            idx <- unique(vapply(spectra, fromFile, integer(1L)))
            object@modCount[idx] <- object@modCount[idx] + 1L
        }
        validObject(object)
        object
})

## #' @rdname hidden_aliases
## setMethod("backendUpdateMetadata", "BackendMemory", function(object,
##                                                              spectraData) {
##     object@spectra <- object@spectra[rownames(spectraData)]
##     object@spectra <- mapply(object@spectra,
##                              split(spectraData, seq_len(nrow(spectraData))),
##                              FUN = .spectrum_set_header)
##     validObject(object)
##     object
## })
