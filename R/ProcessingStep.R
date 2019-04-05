#' @include hidden_aliases.R
NULL

setClassUnion("characterOrFunction", c("character", "function"))

#' @title Processing step
#'
#' @description
#'
#' Class containing the function and arguments to be applied in a lazy-execution
#' framework.
#'
#' Objects of this class are created using the `ProcessingStep()` function. The
#' processing step is executed with the `executeProcessingStep()` function.
#'
#' @details
#'
#' This object contains all relevant information of a data analysis processing
#' step, i.e. the function and all of its arguments to be applied to the data.
#' This object is mainly used to record possible processing steps of a
#' [Spectra()] object.
#'
#' @author Johannes Rainer
#'
#' @exportClass ProcessingStep
#'
#' @md
#'
#' @examples
#'
#' ## Create a simple processing step object
#' ps <- ProcessingStep(sum)
#'
#' executeProcessingStep(ps, 1:10)
#'
#' @name ProcessingStep
NULL

#' @rdname hidden_aliases
setClass("ProcessingStep",
         representation = representation(
             FUN = "characterOrFunction",
             ARGS = "list"
         ),
         prototype = prototype(
             ARGS = list(),
             FUN = character()
         ),
         validity = function(object) {
             msg <- character()
             if (length(object@FUN)) {
                 if (!is.function(object@FUN)) {
                     res <- try(match.fun(object@FUN), silent = TRUE)
                     if (is(res, "try-error"))
                         msg <- c(msg, paste0("Function '", object@FUN,
                                              "' not found."))
                 }
             }
             if (length(msg))
                 msg
             else TRUE
         })

#' @param FUN `function` or `character` representing a function name.
#'
#' @param ARGS `list` of arguments to be passed along to `FUN`.
#'
#' @rdname ProcessingStep
#'
#' @md
#'
#' @export
ProcessingStep <- function(FUN = character(), ARGS = list())  {
    if (missing(FUN))
        FUN <- character()
    new("ProcessingStep", FUN = FUN, ARGS = ARGS)
}

#' @rdname hidden_aliases
#'
#' @exportMethod show
setMethod("show", "ProcessingStep", function(object) {
    cat("Object of class \"", class(object), "\"\n", sep = "")
    cat(" Function: ", object@FUN, "\n", sep = "")
    args <- object@ARGS
    if (length(args) > 0) {
        cat(" Arguments:\n")
        for (i in 1:length(args)) {
            cat("  o ", names(args)[i], " = ",
                paste(args[[i]], collapse = ", "), "\n", sep = "")
        }
    }
})

#' @param object `ProcessingStep` object.
#'
#' @param ... optional additional arguments to be passed along.
#'
#' @rdname ProcessingStep
#'
#' @md
#'
#' @export
executeProcessingStep <- function(object, ...) {
    if (!is(object, "ProcessingStep"))
        stop("'object' is supposed to be a 'ProcessingStep' object!")
    do.call(object@FUN, args = c(list(...), object@ARGS))
}

#' Internal function to apply the lazy processing queue to each spectrum
#' in the provided list.
#'
#' @param x `list` of `Spectrum` objects.
#'
#' @param queue `list` (or `NULL`) of `ProcessingStep` objects.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.apply_processing_queue <- function(x, queue = NULL) {
    if (length(queue)) {
        if (!is.list(x))
            x <- list(x)
        x <- lapply(x, function(z, q) {
            for (pStep in q) {
                z <- executeProcessingStep(pStep, z)
            }
            z
        }, q = queue)
    }
    x
}
