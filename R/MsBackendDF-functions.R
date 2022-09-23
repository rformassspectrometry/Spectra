.df_pdata_column <- function(x, column) {
    idx <- which(colnames(x[[1L]]) == column)
    if (length(idx)) {
        lapply(x, `[`, , j = idx)
    } else stop("No peaks variable \"", column, "\" available")
}

#' Check columns of a DataFrame or data.frame for their data type. If they
#' are `list` or inherit `List`, and they have the same lengths than a column
#' mz or intensity return them.
#'
#' Note that the function will not return any column if there is not at least
#' one column called mz or intensity.
#'
#' @noRd
.get_peaks_columns_data_frame <- function(x) {
    lns <- integer()
    if (any(colnames(x) == "mz"))
        lns <- lengths(x$mz)
    if (!length(lns) && any(colnames(x) == "intensity"))
        lns <- lengths(x$intensity)
    if (length(lns)) {
        colnames(x)[vapply(x, function(z) {
            (inherits(z, "NumericList") || inherits(z, "SimpleList") ||
             is.list(z)) && all(lengths(z) == lns)
        }, logical(1))]
    } else character()
}
