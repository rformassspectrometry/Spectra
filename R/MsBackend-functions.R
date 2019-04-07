#' @include hidden_aliases.R
NULL

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

.SPECTRA_DATA_COLUMNS <- c(
    msLevel = "integer",
    rt = "numeric",
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

## MsBackend <- function(spectraData, msLevel, rt, mz, intensity, acquisitionNum,
##                       scanIndex, fromFile, centroided, smoothed, polarity,
##                       precScanNum, precursorMz, precursorIntensity,
##                       precursorCharge, collisionEnergy) {

## }
