#' @include hidden_aliases.R
NULL

.valid_ms_backend_mod_count <- function(x, y) {
    if (length(x) != length(y))
        "Different number of source files and modification counters."
    else
        NULL
}

.valid_ms_backend_files_from_file <- function(x, y) {
    if (length(x) && !length(y))
        return("'fromFile' can not be empty if 'files' are defined.")
    if (length(y) && !length(x))
        return("'files' can not be empty if 'fromFile' is defined.")
    if (length(x) && !all(y %in% seq_along(x)))
        return("Index in 'fromFile' outside of the number of files")
    else NULL
}

#' @description
#'
#' Check if spectraData has all required columns.
#'
#' @noRd
#'
#' @param x spectraData `DataFrame`
.valid_spectra_data_required_columns <- function(x, columns = c("fromFile")) {
    if (nrow(x)) {
        missing_cn <- setdiff(columns, colnames(x))
        if (length(missing_cn))
            return(paste0("Required column(s): ",
                          paste(missing_cn, collapse = ", "),
                          " is/are missing"))
    }
    NULL
}

#' Function to check data types of selected columns in the provided `DataFrame`.
#'
#' @param x `DataFrame` to validate.
#'
#' @param datatypes named `character`, names being column names and elements
#'     expected data types.
#'
#' @author Johannes Rainer
#'
#' @noRd
.valid_column_datatype <- function(x, datatypes = .SPECTRA_DATA_COLUMNS) {
    datatypes <- datatypes[names(datatypes) %in% colnames(x)]
    x_class <- vapply(x, class, character(1))[names(datatypes)]
    if (!all(datatypes == x_class))
        paste0("The following columns have a wrong data type: ",
               paste(names(datatypes)[datatypes != x_class], collapse = ", "),
               ". The expected data type(s) is/are: ",
               paste(datatypes[datatypes != x_class], collapse = ", "), ".")
    else NULL
}

.valid_mz_column <- function(x) {
    if (length(x$mz)) {
        if (!all(vapply(x$mz, is.numeric, logical(1))))
            return("mz column should contain a list of numeric")
        if (any(vapply(x$mz, is.unsorted, logical(1))))
            return("mz values have to be sorted increasingly")
    }
    NULL
}

.valid_intensity_column <- function(x) {
    if (length(x$intensity))
        if (!all(vapply(x$intensity, is.numeric, logical(1))))
            return("intensity column should contain a list of numeric")
    NULL
}

.valid_intensity_mz_columns <- function(x) {
    if (any(lengths(mz(x)) != lengths(intensity(x))))
        "Length of mz and intensity values differ for some spectra"
    else NULL
}

#' data types of spectraData columns
#'
#' @noRd
.SPECTRA_DATA_COLUMNS <- c(
    msLevel = "integer",
    rtime = "numeric",
    acquisitionNum = "integer",
    scanIndex = "integer",
    mz = "list",
    intensity = "list",
    fromFile = "integer",
    centroided = "logical",
    smoothed = "logical",
    polarity = "integer",
    precScanNum = "integer",
    precursorMz = "numeric",
    precursorIntensity = "numeric",
    precursorCharge = "integer",
    collisionEnergy = "numeric"
)

#' accessor methods for spectraData columns.
#'
#' @noRd
.SPECTRA_DATA_COLUMN_METHODS <- c(
    msLevel = "msLevel",
    rtime = "rtime",
    acquisitionNum = "acquisitionNum",
    scanIndex = "scanIndex",
    mz = "mz",
    intensity = "intensity",
    fromFile = "fromFile",
    centroided = "centroided",
    smoothed = "smoothed",
    polarity = "polarity",
    precScanNum = "precScanNum",
    precursorMz = "precursorMz",
    precursorIntensity = "precursorIntensity",
    precursorCharge = "precursorCharge",
    collisionEnergy = "collisionEnergy"
)

#' @rdname MsBackend
#'
#' @export MsBackendDataFrame
MsBackendDataFrame <- function() {
    new("MsBackendDataFrame")
}

#' Helper function to extract a certain column from the spectraData data frame.
#' If the data frame has no such column it will use the accessor method to
#' retrieve the corresponding data.
#'
#' @param x object with a `@spectraData` slot containing a `DataFrame`.
#'
#' @param column `character(1)` with the column name.
#'
#' @author Johannes Rainer
#'
#' @noRd
.get_spectra_data_column <- function(x, column) {
    if (missing(column) || length(column) != 1)
        stop("'column' should be a 'character' of length 1.")
    if (any(colnames(x@spectraData) == column))
        x@spectraData[, column]
    else {
        if (any(names(.SPECTRA_DATA_COLUMN_METHODS) == column))
            do.call(.SPECTRA_DATA_COLUMN_METHODS[column], args = list(x))
        else stop("No column '", column, "' available")
    }
}
