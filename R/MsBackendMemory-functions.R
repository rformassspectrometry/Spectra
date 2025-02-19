#' @rdname MsBackend
#'
#' @export MsBackendMemory
MsBackendMemory <- function() {
    new("MsBackendMemory")
}

.df_pdata_column <- function(x, column) {
    idx <- which(colnames(x[[1L]]) == column)
    if (length(idx))
        lapply(x, `[`, , j = idx)
    else {
        if (column %in% c("mz", "intensity"))
            NumericList(vector("list", length(x)), compress = FALSE)
        else stop("No peaks variable \"", column, "\" available")
    }
}

.df_empty_peaks_data <- function(x) {
    emat <- matrix(numeric(), ncol = 2, nrow = 0,
                   dimnames = list(character(), c("mz", "intensity")))
    replicate(x, emat, simplify = FALSE)
}

#' Check columns of a DataFrame or data.frame for their data type. If they
#' are `list` or inherit `List`, and they have the same lengths as a column
#' mz or intensity return them.
#'
#' Note that the function will not return any column if there is not at least
#' one column called mz or intensity.
#'
#' @noRd
.df_peaks_columns_data_frame <- function(x) {
    lns <- lengths(x$mz)
    if (!length(lns) && any(colnames(x) == "intensity"))
        lns <- lengths(x$intensity)
    if (length(lns)) {
        colnames(x)[vapply(x, function(z) {
            (is.list(z) || inherits(z, "NumericList") ||
             inherits(z, "SimpleList")) && all(lengths(z) == lns)
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
            res <- fillCoreSpectraVariables(res, other_columns)
            ## other_res <- lapply(other_columns, .get_column,
            ##                     x = object@spectraData)
            ## names(other_res) <- other_columns
            ## res <- cbind(res, as.data.frame(other_res))
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

#' Helper function to combine/concatenate backends that base on
#' [MsBackendMemory()].
#'
#' @param objects `list` of `MsBackend` objects. Note that this function does
#'     not expect **empty** objects!
#'
#' @return [MsBackend()] object with combined content.
#'
#' @author Johannes Rainer
#'
#' @importFrom MsCoreUtils vapply1c rbindFill
#'
#' @noRd
.df_combine <- function(objects) {
    if (length(objects) == 1)
        return(objects[[1]])
    if (!all(vapply1c(objects, class) == class(objects[[1]])))
        stop("Can only merge backends of the same type: ", class(objects[[1]]))
    res <- objects[[1]]
    pv <- peaksVariables(res)
    for (i in 2:length(objects)) {
        res@spectraData <- rbindFill(res@spectraData, objects[[i]]@spectraData)
        pv2 <- peaksVariables(objects[[i]])
        if (length(pv) == length(pv2) && all(pv == pv2)) {
            res@peaksData <- c(res@peaksData, objects[[i]]@peaksData)
            res@peaksDataFrame <- c(res@peaksDataFrame,
                                    objects[[i]]@peaksDataFrame)
        } else
            stop("Provided objects have different sets of peak variables. ",
                 "Combining such objects is currently not supported.")
    }
    res
}

.check_peaks_data_value <- function(x, lo) {
    if (!(is.list(x) || inherits(x, "SimpleList")))
        stop("'value' has to be a list-like object")
    if (length(x) != lo)
        stop("Length of 'value' has to match length of 'object'")
    if (!(is.matrix(x[[1L]]) | is.data.frame(x[[1L]])))
        stop("'value' is expected to be a 'list' of 'matrix' ",
             "or 'data.frame'")
    cn <- colnames(x[[1L]])
    lcn <- length(cn)
    lapply(x, function(z) {
        cur_cn <- colnames(z)
        if (lcn != length(cur_cn) || !all(cn == cur_cn))
            stop("provided matrices don't have the same column names")
    })
}
