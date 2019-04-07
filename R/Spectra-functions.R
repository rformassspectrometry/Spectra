#' @include hidden_aliases.R
NULL

.valid_processing_queue <- function(x) {
    if (length(x))
        if (!all(vapply(x, inherits, logical(1), "ProcessingStep")))
            return("'processingQueue' should only contain ProcessingStep objects.")
    NULL
}

.valid_spectra <- function(object) {
    msg <- c(.valid_processing_queue(x@processingQueue))
    if (length(msg)) msg
    else NULL
}
