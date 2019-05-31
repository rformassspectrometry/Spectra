#' @include hidden_aliases.R
NULL

#' @rdname MsBackend
#'
#' @export MsBackendMzR
MsBackendMzR <- function() {
    if (!requireNamespace("mzR", quietly = TRUE))
        stop("The use of 'MsBackendMzR' requires package 'mzR'. Please ",
             "install with 'Biobase::install(\"mzR\")'")
    new("MsBackendMzR")
}

#' Read the header for each spectrum from the MS file `x`
#'
#' @author Johannes Rainer
#'
#' @return `DataFrame` with the header.
#'
#' @noRd
.mzR_header <- function(x = character()) {
    if (length(x) != 1)
        stop("'x' should have length 1")
    requireNamespace("mzR", quietly = TRUE)
    msd <- mzR::openMSfile(x)
    on.exit(mzR::close(msd))
    hdr <- mzR::header(msd)
    colnames(hdr)[colnames(hdr) == "seqNum"] <- "scanIndex"
    colnames(hdr)[colnames(hdr) == "precursorScanNum"] <- "precScanNum"
    colnames(hdr)[colnames(hdr) == "precursorMZ"] <- "precursorMz"
    colnames(hdr)[colnames(hdr) == "retentionTime"] <- "rtime"
    colnames(hdr)[colnames(hdr) == "isolationWindowTargetMZ"] <-
        "isolationWindowTargetMz"
    hdr$isolationWindowLowerMz <- hdr$isolationWindowTargetMz -
        hdr$isolationWindowLowerOffset
    hdr$isolationWindowUpperMz <- hdr$isolationWindowTargetMz +
        hdr$isolationWindowUpperOffset
    hdr$isolationWindowUpperOffset <- NULL
    hdr$isolationWindowLowerOffset <- NULL
    S4Vectors::DataFrame(hdr)
}

#' Read peaks from a single mzML file.
#'
#' @param x `character(1)` with the file to read from.
#'
#' @param scanIndex (required) indices of spectra from which the data should be
#'     retrieved.
#'
#' @author Johannes Rainer
#'
#' @return list of `matrix`
#'
#' @noRd
.mzR_peaks <- function(x = character(), scanIndex = integer()) {
    if (length(x) != 1)
        stop("'x' should have length 1")
    msd <- mzR::openMSfile(x)
    on.exit(mzR::close(msd))
    ## hd_spectra <- mzR::header(msd, max(scanIndex))
    pks <- mzR::peaks(msd, scanIndex)
    if (is.matrix(pks))
        pks <- list(pks)
    lapply(pks, function(z) {
        colnames(z) <- c("mz", "intensity")
        z
    })
}

#' Utility function to convert columns in the `x` `DataFrame` that have only
#' a single element to `Rle`. Also columns specified with parameter `columns`
#' will be converted (if present).
#'
#' @param x `DataFrame`
#'
#' @param columns `character` of column names that should be converted to `Rle`
#'
#' @return `DataFrame`
#'
#' @author Johannes Rainer
#'
#' @importClassesFrom S4Vectors Rle
#'
#' @importFrom S4Vectors Rle
#'
#' @noRd
.as_rle_spectra_data <- function(x, columns = c("fromFile")) {
    if (nrow(x) <= 1)
        return(x)
    for (col in colnames(x)) {
        x[[col]] <- .as_rle(x[[col]])
    }
    columns <- intersect(columns, colnames(x))
    for (col in columns) {
        if (!is(x[[col]], "Rle"))
            x[[col]] <- Rle(x[[col]])
    }
    x
}


#' *Uncompress* a `DataFrame` by converting all of its columns that contain
#' an `Rle` with `as.vector`.
#'
#' @param x `DataFrame`
#'
#' @return `DataFrame` with all `Rle` columns converted to vectors.
#'
#' @author Johannes Rainer
#'
#' @noRd
.as_vector_spectra_data <- function(x) {
    cols <- colnames(x)[vapply(x, is, logical(1), "Rle")]
    for (col in cols) {
        x[[col]] <- as.vector(x[[col]])
    }
    x
}

#' Helper function to return a column from the (spectra data) `DataFrame`. If
#' the column `column` is an `Rle` `as.vector` is called on it. If column is
#' the name of a mandatory variable but it is not available it is created on
#' the fly.
#'
#' @param x `DataFrame`
#'
#' @param column `character(1)` with the name of the column to return.
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

#' @importFrom IRanges NumericList
#'
#' @description
#'
#' Helper to build the spectraData `DataFrame` for `MsBackendMzR` backend.
#'
#' @noRd
.spectra_data_mzR <- function(x, columns = spectraVariables(x)) {
    cn <- colnames(x@spectraData)
    if(!nrow(x@spectraData)) {
        res <- lapply(
            .SPECTRA_DATA_COLUMNS[!(names(.SPECTRA_DATA_COLUMNS) %in%
                                    c("mz", "intensity"))], do.call,
            args = list())
        res <- DataFrame(res)
        res$mz <- NumericList(compress = FALSE)
        res$intensity <- NumericList(compress = FALSE)
        return(res[, columns, drop = FALSE])
    }
    not_found <- setdiff(columns, c(cn, names(.SPECTRA_DATA_COLUMNS)))
    if (length(not_found))
        stop("Column(s) ", paste(not_found, collapse = ", "),
             " not available")
    sp_cols <- columns[columns %in% cn]
    res <- .as_vector_spectra_data(
        x@spectraData[, sp_cols, drop = FALSE])
    if (any(columns %in% c("mz", "intensity"))) {
        pks <- peaks(x)
        if (any(columns == "mz"))
            res$mz <- NumericList(lapply(pks, function(z) z[, 1]),
                                  compress = FALSE)
        if (any(columns == "intensity"))
            res$intensity <- NumericList(lapply(pks,
                                                function(z) z[, 2]),
                                         compress = FALSE)
    }
    other_cols <- setdiff(
        columns[!(columns %in% c("mz", "intensity"))], sp_cols)
    if (length(other_cols)) {
        other_res <- lapply(other_cols, .get_spectra_data_column,
                            x = x)
        names(other_res) <- other_cols
        res <- cbind(res, as(other_res, "DataFrame"))
    }
    res[, columns, drop = FALSE]
}
