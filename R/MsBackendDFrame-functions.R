#' @include hidden_aliases.R
NULL

#' @rdname MsBackend
#'
#' @export MsBackendDFrame
MsBackendDFrame <- function() {
    new("MsBackendDFrame")
}

#' Helper function to return a column from the (spectra data) `DataFrame`. If
#' the column `column` is an `Rle` `as.vector` is called on it. If column is
#' the name of a mandatory variable but it is not available it is created on
#' the fly.
#'
#' @note This function is equivalent to the `get_rle_column` function in
#'     the `Chromatograms` package (in *ChromBackendDataFrame-functions.R*).
#'
#' @param x `DataFrame`
#'
#' @param column `character(1)` with the name of the column to return.
#'
#' @importClassesFrom S4Vectors Rle
#'
#' @importFrom S4Vectors Rle
#'
#' @importMethodsFrom S4Vectors [[
#'
#' @author Johannes Rainer
#'
#' @noRd
.get_rle_column <- function(x, column) {
    if (any(colnames(x) == column)) {
        if (is(x[[column]], "Rle"))
            as.vector(x[[column]])
        else x[[column]]
    } else if (any(names(.SPECTRA_DATA_COLUMNS) == column)) {
        nr_x <- nrow(x)
        if (nr_x)
            as(rep(NA, nr_x), .SPECTRA_DATA_COLUMNS[column])
        else
            do.call(.SPECTRA_DATA_COLUMNS[column], args = list())
    } else stop("column '", column, "' not available")
}

#' Helper function to combine backends that base on [MsBackendDFrame()].
#'
#' @param objects `list` of `MsBackend` objects.
#'
#' @return [MsBackend()] object with combined content.
#'
#' @author Johannes Rainer
#'
#' @importFrom MsCoreUtils vapply1c asRleDataFrame rbindFill
#'
#' @noRd
.combine_backend_dframe <- function(objects) {
    if (length(objects) == 1)
        return(objects[[1]])
    if (!all(vapply1c(objects, class) == class(objects[[1]])))
        stop("Can only merge backends of the same type: ", class(objects[[1]]))
    res <- new(class(objects[[1]]))
    suppressWarnings(
        res@spectraData <- asRleDataFrame(do.call(
            rbindFill, lapply(objects, function(z) z@spectraData)),
            columns = c("dataStorage", "dataOrigin"))
    )
    res
}
