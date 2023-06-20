#' @include hidden_aliases.R
NULL

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
.valid_spectra_data_required_columns <-
    function(x, columns = c("dataStorage")) {
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
    res <- mapply(FUN = is, x[, names(datatypes), drop = FALSE],
                  datatypes)
    if (!all(res))
        paste0("The following columns have a wrong data type: ",
               paste(names(res[!res]), collapse = ", "),
               ". The expected data type(s) is/are: ",
               paste(datatypes[names(res)[!res]], collapse = ", "), ".")
    else NULL
}

.valid_mz_column <- function(x) {
    if (length(x$mz)) {
        if (!inherits(x$mz, "NumericList"))
            return("mz column should be of type NumericList")
        if (any(vapply1l(x$mz, is.unsorted)))
            return("mz values have to be sorted increasingly")
    }
    NULL
}

.valid_intensity_column <- function(x) {
    if (length(x$intensity))
        if (!inherits(x$intensity, "NumericList"))
            return("intensity column should be of type NumericList")
    NULL
}

.peaks_variables <- function(x) {
    if (.hasSlot(x, "peaksVariables")) {
        x@peaksVariables
    } else c("mz", "intensity")
}

.valid_peaks_variable_columns <- function(x, pvars) {
    lens <- lapply(pvars, function(z) lengths(x[[z]]))
    names(lens) <- pvars
    lens <- lens[lengths(lens) > 0]
    if (length(lens) > 1) {
        for (i in 2:length(lens))
            if (any(lens[[1L]] != lens[[i]]))
                return(paste0("Number of values per spectra differ for peak ",
                              "variables \"", names(lens)[1L], "\" and \"",
                              names(lens)[i], "\"."))
    }
    NULL
}

.valid_intensity_mz_columns <- function(x) {
    ## Don't want to have that tested on all on-disk objects.
    if (length(x$intensity) && length(x$mz))
        if (any(lengths(x$mz) != lengths(x$intensity)))
            return("Length of mz and intensity values differ for some spectra")
    NULL
}

#' data types of spectraData columns
#'
#' @noRd
.SPECTRA_DATA_COLUMNS <- c(
    msLevel = "integer",
    rtime = "numeric",
    acquisitionNum = "integer",
    scanIndex = "integer",
    mz = "NumericList",
    intensity = "NumericList",
    dataStorage = "character",
    dataOrigin = "character",
    centroided = "logical",
    smoothed = "logical",
    polarity = "integer",
    precScanNum = "integer",
    precursorMz = "numeric",
    precursorIntensity = "numeric",
    precursorCharge = "integer",
    collisionEnergy = "numeric",
    isolationWindowLowerMz = "numeric",
    isolationWindowTargetMz = "numeric",
    isolationWindowUpperMz = "numeric"
)

#' @rdname MsBackend
#'
#' @export MsBackendDataFrame
MsBackendDataFrame <- function() {
    new("MsBackendDataFrame")
}

#' Helper function to return a column from the (spectra data) `DataFrame`.
#' If column is the name of a mandatory variable but it is not available
#' it is created on the fly.
#'
#' @note This function is equivalent to the `get_column` function in
#'     the `Chromatograms` package (in *ChromBackendDataFrame-functions.R*).
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
.get_column <- function(x, column) {
    if (any(colnames(x) == column)) {
        x[[column]]
    } else if (any(names(.SPECTRA_DATA_COLUMNS) == column)) {
        nr_x <- nrow(x)
        if (nr_x)
            as(rep(NA, nr_x), .SPECTRA_DATA_COLUMNS[column])
        else
            do.call(.SPECTRA_DATA_COLUMNS[column], args = list())
    } else stop("column '", column, "' not available")
}

#' @description
#'
#' Helper to be used in the filter functions to select the file/origin in
#' which the filtering should be performed.
#'
#' @param object `MsBackend`
#'
#' @param dataStorage `character` or `integer` with either the names of the
#'     `dataStorage` or their index (in `unique(object$dataStorage)`) in which
#'     the filtering should be performed.
#'
#' @param dataOrigin same as `dataStorage`, but for the `dataOrigin` spectra
#'     variable.
#'
#' @return `logical` of length equal to the number of spectra in `object`.
#'
#' @noRd
.sel_file <- function(object, dataStorage = integer(), dataOrigin = integer()) {
    if (length(dataStorage)) {
        lvls <- unique(object@spectraData$dataStorage)
        if (!(is.numeric(dataStorage) || is.character(dataStorage)))
            stop("'dataStorage' has to be either an integer with the index of",
                 " the data storage, or its name")
        if (is.numeric(dataStorage)) {
            if (dataStorage < 1 || dataStorage > length(lvls))
                stop("'dataStorage' should be an integer between 1 and ",
                     length(lvls))
            dataStorage <- lvls[dataStorage]
        }
        dataStorage(object) %in% dataStorage
    } else if (length(dataOrigin)) {
        lvls <- unique(object@spectraData$dataOrigin)
        if (!(is.numeric(dataOrigin) || is.character(dataOrigin)))
            stop("'dataOrigin' has to be either an integer with the index of",
                 " the data origin, or its name")
        if (is.numeric(dataOrigin)) {
            if (dataOrigin < 1 || dataOrigin > length(lvls))
                stop("'dataOrigin' should be an integer between 1 and ",
                     length(lvls))
            dataOrigin <- lvls[dataOrigin]
        }
        dataOrigin(object) %in% dataOrigin
    } else rep(TRUE, length(object))
}

#' Helper function to combine backends that base on [MsBackendDataFrame()].
#'
#' @param objects `list` of `MsBackend` objects.
#'
#' @return [MsBackend()] object with combined content.
#'
#' @author Johannes Rainer
#'
#' @importFrom MsCoreUtils vapply1c rbindFill
#'
#' @noRd
.combine_backend_data_frame <- function(objects) {
    if (length(objects) == 1)
        return(objects[[1]])
    if (!all(vapply1c(objects, class) == class(objects[[1]])))
        stop("Can only merge backends of the same type: ", class(objects[[1]]))
    pvars <- lapply(objects, peaksVariables)
    for (i in 2:length(pvars))
        if (length(pvars[[i]]) != length(pvars[[1L]]) ||
            any(pvars[[i]] != pvars[[1L]]))
            stop("Provided backends have different peaks variables. Can only ",
                 "merge backends with the same set of peaks variables.")
    res <- objects[[1L]]
    suppressWarnings(
        res@spectraData <- do.call(
            rbindFill, lapply(objects, function(z) z@spectraData))
    )
    res
}

#' @description
#'
#' Subset the `DataFrame` of a `MsBackendDataFrame` *by rows*.
#'
#' @param x `MsBackendDataFrame
#'
#' @param i `integer`, `character` or `logical`.
#'
#' @return Subsetted `x`
#'
#' @author Johannes Rainer
#'
#' @importMethodsFrom S4Vectors extractROWS
#'
#' @importFrom methods slot<-
#'
#' @noRd
.subset_backend_data_frame <- function(x, i) {
    if (missing(i))
        return(x)
    i <- i2index(i, length(x), rownames(x@spectraData))
    slot(x, "spectraData", check = FALSE) <- extractROWS(x@spectraData, i)
    x
}
