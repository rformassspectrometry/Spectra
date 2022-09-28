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
.df_peaks_columns_data_frame <- function(x) {
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

.df_spectra_data <- function(object, columns = spectraVariables(object)) {
        if (!all(columns %in% spectraVariables(object)))
            stop("Some of the requested spectra variables are not available")
        p_vars <- peaksVariables(object)
        sp_vars <- setdiff(columns, p_vars)
        df_columns <- intersect(sp_vars, colnames(object@spectraData))
        res <- object@spectraData[, df_columns, drop = FALSE]
        ## Get missing core variables.
        other_columns <- setdiff(sp_vars, colnames(object@spectraData))
        if (length(other_columns)) {
            other_res <- lapply(other_columns, .get_column,
                                x = object@spectraData)
            names(other_res) <- other_columns
            res <- cbind(res, as.data.frame(other_res))
        }
        if (any(columns == "mz"))
            res$mz <- mz(object)
        if (any(columns == "intensity"))
            res$intensity <- intensity(object)
        p_columns <- setdiff(p_vars, c("mz", "intensity"))
        for (p_column in p_columns) {
            res <- do.call(
                `$<-`, list(res, name = p_column,
                            value = .df_pdata_column(object@peaksDataFrame,
                                                     p_column)))
        }
        res[, columns, drop = FALSE]
}

.df_subset <- function(x, i) {
    if (missing(i))
        return(x)
    i <- i2index(i, length(x), rownames(x@spectraData))
    slot(x, "spectraData", check = FALSE) <- x@spectraData[i, , drop = FALSE]
    if (length(x@peaksData))
        slot(x, "peaksData", check = FALSE) <- x@peaksData[i]
    if (length(x@peaksDataFrame))
        slot(x, "peaksDataFrame", check = FALSE) <- x@peaksDataFrame[i]
    x
}
