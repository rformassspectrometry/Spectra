#' @include hidden_aliases.R
NULL

.valid_processing_queue <- function(x) {
    if (length(x))
        if (!all(vapply(x, inherits, logical(1), "ProcessingStep")))
            return("'processingQueue' should only contain ProcessingStep objects.")
    NULL
}

#' @description
#'
#' Combine two `DataFrame`s. The resulting `DataFrame` contains all
#' columns from both `x` and `y`. For columns present in both `DataFrame`s those
#' in `x` will be used. Also, the resulting `DataFrame` uses the row names of
#' `x` unless `x` has no row names.
#'
#' @param x `DataFrame`
#'
#' @param y `DataFrame`
#'
#' @return `DataFrame` with the merged columns.
#'
#' @author Johannes Rainer
#'
#' @importMethodsFrom S4Vectors nrow rownames colnames cbind
#'
#' @noRd
.combine_data_frame <- function(x, y) {
    if (nrow(y) == 0)
        return(x)
    if (nrow(x) == 0)
        return(y)
    if (nrow(x) != nrow(y))
        stop("'x' and 'y' have to have the same number of rows")
    if (is.null(rownames(x)) & !is.null(rownames(y)))
        rownames(x) <- rownames(y)
    cols_y <- !(colnames(y) %in% colnames(x))
    if (any(cols_y))
        x <- cbind(x, y[, cols_y, drop = FALSE])
    x
}
