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
.compress_spectra_data <- function(x, columns = c("fromFile")) {
    if (nrow(x) <= 1)
        return(x)
    for (col in colnames(x)) {
        x[[col]] <- .rle_compress(x[[col]])
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
.uncompress_spectra_data <- function(x) {
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
