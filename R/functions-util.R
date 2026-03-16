#' @param acquisitionNum `integer` with the acquisition numbers of all scans.
#'
#' @param precursorScanNum `integer` with the precursor scan numbers.
#'
#' @param an integer, acquisitionNum of spectrum of interest (parent and
#' children will be selected)
#'
#' @author Sebastian Gibb
#'
#' @noRd
.filterSpectraHierarchy <- function(acquisitionNum = integer(),
                                    precursorScanNum = integer(), an) {
    if (length(acquisitionNum) != length(precursorScanNum))
        stop("length of 'acquisitionNum' and 'precursorScanNum' have to be ",
             "the same")
    ## we could use recursion which is slow in R
    ## or reformat the adjacency list into a nested tree
    ## list model but most ms data are limited to at most 3 levels and the
    ## filtering isn't done very often, so we use for loops here

    parents <- logical(length(acquisitionNum))

    ## find current scan
    parents[acquisitionNum %in% an] <- TRUE
    children <- parents

    ## find parent scan
    nLastParents <- 0L
    nParents <- 1L
    while (nLastParents < nParents) {
        parents[acquisitionNum %in% precursorScanNum[parents]] <- TRUE
        nLastParents <- nParents
        nParents <- sum(parents)
    }

    ## find children scans
    nLastChildren <- 0L
    nChildren <- 1L
    while (nLastChildren < nChildren) {
        children[precursorScanNum %in% acquisitionNum[children]] <- TRUE
        nLastChildren <- nChildren
        nChildren <- sum(children)
    }
    parents | children
}

.logging <- function(x, ...) {
    c(x, paste0(..., " [", date(), "]"))
}

setAs("logical", "factor", function(from, to) factor(from))

#' @description
#'
#' Simple helper function to sanitize the file name. The function uses the
#' [path_sanitize()] function only on the file name but not on the path.
#'
#' @param x `character` with file names/paths.
#'
#' @noRd
#'
#' @importFrom fs path_sanitize
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' a <- c(tempfile(), tempfile(), ">hello")
#' sanitize_file_name(a)
sanitize_file_name <- function(x) {
    file.path(normalizePath(dirname(x)), path_sanitize(basename(x)))
}

#' Helper function that matches `x` against `mz` (using the `closest` function)
#' and returns the indices of `x` that match any of the values in `mz`. The
#' function takes care of sorting `x` and `mz` and deals also with missing
#' values.
#'
#' @return `integer` with the indices of values in `x` that are not `NA` and
#'     are matching any of the values in `mz` given `ppm` and `tolerance`.
#'
#' @noRd
#'
#' @author Sebastian Gibb, Johannes Rainer
.values_match_mz <- function(x, mz, ppm = 20, tolerance = 0) {
    o <- order(x, na.last = NA)
    cmn <- common(x[o], sort(mz), tolerance = tolerance, ppm = ppm,
                  duplicates = "keep", .check = FALSE)
    sort(o[cmn])
}

#' Convert a `spectraData` data.frame to long format.
#'
#' @param spectraData `data.frame` with the spectra and peaks data as e.g.,
#'     returned by a call to `spectraData()` on a `MsBackend`. For performance
#'     reasons it should ideally be a `data.frame`, not `DataFrame`.
#'
#' @param peaksVariables `character` with the names of the columns in
#'     `spectraData` containing the peaks data (as `list` of `numeric`). If
#'     `length(peaksVariables) == 0` `spectraData` is returned *as-is*.
#'
#' @return `data.frame` in long format.
#'
#' @noRd
.long_spectra_data2 <- function(spectraData,
                                peaksVariables = c("mz", "intensity")) {
    peaksVariables <- intersect(colnames(spectraData), peaksVariables)
    if (!length(peaksVariables))
        return(spectraData)
    ls <- lengths(spectraData[[peaksVariables[1]]])
    res <- lapply(spectraData[, !colnames(spectraData) %in% peaksVariables,
                              drop = FALSE], rep, times = ls)
    for (i in peaksVariables)
        res[[i]] <- unlist(spectraData[[i]], use.names = FALSE,
                           recursive = FALSE)
    base::as.data.frame(res[colnames(spectraData)])
}

#' Convert `spectraData` and `peaksData` to a `data.frame` in long format.
#' Ideally, `spectraData` and `peaksData` should be base R data types
#' (`data.frame` and `list` and not `DataFrame` and `List`).
#'
#' @param spectraData `data.frame` with the data on spectra variables.
#'     Should **not** include peaks variables.
#'
#' @param peaksData `list` of `numeric` matrices with peak data.
#'
#' @param peaksVariables character with the names of the peaks variables. These
#'     **have** to match the column names of the peak matrices in `peaksData`.
#'
#' @return `data.frame` in long format.
#'
#' @noRd
.long_spectra_data3 <- function(spectraData, peaksData,
                                peaksVariables = colnames(peaksData[[1L]])) {
    ls <- lengths(peaksData) / length(peaksVariables)
    if (ncol(spectraData))
        cbind.data.frame(
            base::as.data.frame(lapply(spectraData, rep, times = ls)),
            base::as.data.frame(do.call(base::rbind, peaksData)))
    else base::as.data.frame(do.call(base::rbind, peaksData))
}

#' @title Fast *rbind-ing* `data.frame`s preserving row names
#'
#' @description
#'
#' The `rbindlistWithRownames()` function uses the [data.table::rbindlist()]
#' function for a fast concatenation (row-binding) of a `list` of provided
#' `data.frame`s and in addition also preserves their row names (if set).
#'
#' The parameters are the same as for `rbindlist()`.
#'
#' @param l A `list` containing `data.frame`s that should be combined.
#'
#' @param use.names If `TRUE` binds by matching column names, if `FALSE` by
#'     position. `"check"` (default) warns if not all items have the same names
#'     in the same order. See [data.table::rbindlist()] for details.
#'
#' @param fill `logical(1)`: if `TRUE` fills missing columns with NAs or `NULL`
#'     for missing list columns. By defalt `FALSE`.
#'
#' @param idcol Creates a column in the result showing which list item those
#'     rows cam from. See [data.table::rbindlist()] for details.
#'
#' @param ignore.attr `logical(1)`, default `FALSE`. When `TRUE`, allows
#'     binding columns with different attributes (e.g. class).
#'
#' @note
#'
#' Row names are dropped if duplicated row names are present, or if row
#' names are not defined for all `data.frame`s in `l`.
#'
#' The function uses `.row_names_info()` to guess wheter row names are *really*
#' set or just the default `seq_len(nrow(x))` are used.
#'
#' @return
#'
#' A combined `data.frame` with row names (if present in all `data.frame`s
#' provided).
#'
#' @seealso [data.table::rbindlist()]
#'
#' @author Johannes Rainer
#'
#' @export
rbindlistWithRownames <- function(l, use.names = TRUE, fill = FALSE,
                                  idcol = NULL, ignore.attr = FALSE) {
    rn <- unlist(lapply(l, function(z) {
        if (.row_names_info(z) < 0) NULL else attr(z, "row.names")
    }), use.names = FALSE, recursive = FALSE)
    if (anyDuplicated(rn))
        rn <- rn[1L]
    lrn <- length(rn)
    r <- as.data.frame(rbindlist(l, use.names = use.names, fill = fill,
                                 idcol = idcol, ignore.attr = ignore.attr))
    if (lrn) {
        if (lrn == nrow(r))
            attr(r, "row.names") <- rn
        else warning("Dropping rownames: duplicated rownames present or ",
                     "rownames not available for all data.frames")
    }
    r
}
