#' @include hidden_aliases.R
NULL

## .valid_ms_backend_mod_count <- function(x, y) {
##     if (length(x) != length(y))
##         "Different number of source files and modification counters."
##     else
##         NULL
## }

## .valid_ms_backend_files_from_file <- function(x, y) {
##     if (length(x) && !length(y))
##         return("'fromFile' can not be empty if 'files' are defined.")
##     if (length(y) && !length(x))
##         return("'files' can not be empty if 'fromFile' is defined.")
##     if (length(x) && !all(y %in% seq_along(x)))
##         return("Index in 'fromFile' outside of the number of files")
##     else NULL
## }

## #' @description
## #'
## #' Check if spectraData has all required columns.
## #'
## #' @noRd
## #'
## #' @param x spectraData `DataFrame`
## .valid_spectra_data_required_columns <- function(x, columns = c("fromFile")) {
##     if (nrow(x)) {
##         missing_cn <- setdiff(columns, colnames(x))
##         if (length(missing_cn))
##             return(paste0("Required column(s): ",
##                           paste(missing_cn, collapse = ", "),
##                           " is/are missing"))
##     }
##     NULL
## }

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
.valid_column_rle_datatype <- function(x, datatypes = .SPECTRA_DATA_COLUMNS) {
    datatypes <- datatypes[names(datatypes) %in% colnames(x)]
    for (i in seq_along(datatypes)) {
        if (!.is_class(x[[names(datatypes)[i]]], datatypes[i]))
            return(paste0("Column ", names(datatypes)[i], " has the wrong ",
                          "data type (", class(x[[names(datatypes)[i]]]),
                          ") expected: ", datatypes[i]))
    }
    NULL
}

## .valid_mz_column <- function(x) {
##     if (length(x$mz)) {
##         if (!all(vapply(x$mz, is.numeric, logical(1))))
##             return("mz column should contain a list of numeric")
##         if (any(vapply(x$mz, is.unsorted, logical(1))))
##             return("mz values have to be sorted increasingly")
##     }
##     NULL
## }

## .valid_intensity_column <- function(x) {
##     if (length(x$intensity))
##         if (!all(vapply(x$intensity, is.numeric, logical(1))))
##             return("intensity column should contain a list of numeric")
##     NULL
## }

## .valid_intensity_mz_columns <- function(x) {
##     if (any(lengths(mz(x)) != lengths(intensity(x))))
##         "Length of mz and intensity values differ for some spectra"
##     else NULL
## }

## #' data types of spectraData columns
## #'
## #' @noRd
## .SPECTRA_DATA_COLUMNS <- c(
##     msLevel = "integer",
##     rtime = "numeric",
##     acquisitionNum = "integer",
##     scanIndex = "integer",
##     mz = "list",
##     intensity = "list",
##     fromFile = "integer",
##     centroided = "logical",
##     smoothed = "logical",
##     polarity = "integer",
##     precScanNum = "integer",
##     precursorMz = "numeric",
##     precursorIntensity = "numeric",
##     precursorCharge = "integer",
##     collisionEnergy = "numeric"
## )

## #' accessor methods for spectraData columns.
## #'
## #' @noRd
## .SPECTRA_DATA_COLUMN_METHODS <- c(
##     msLevel = "msLevel",
##     rtime = "rtime",
##     acquisitionNum = "acquisitionNum",
##     scanIndex = "scanIndex",
##     mz = "mz",
##     intensity = "intensity",
##     fromFile = "fromFile",
##     centroided = "centroided",
##     smoothed = "smoothed",
##     polarity = "polarity",
##     precScanNum = "precScanNum",
##     precursorMz = "precursorMz",
##     precursorIntensity = "precursorIntensity",
##     precursorCharge = "precursorCharge",
##     collisionEnergy = "collisionEnergy"
## )

#' @rdname MsBackend
#'
#' @export MsBackendRleDataFrame
MsBackendRleDataFrame <- function() {
    new("MsBackendRleDataFrame")
}

## #' Helper function to extract a certain column from the spectraData data frame.
## #' If the data frame has no such column it will use the accessor method to
## #' retrieve the corresponding data.
## #'
## #' @param x object with a `@spectraData` slot containing a `DataFrame`.
## #'
## #' @param column `character(1)` with the column name.
## #'
## #' @author Johannes Rainer
## #'
## #' @noRd
## .get_spectra_data_column <- function(x, column) {
##     if (missing(column) || length(column) != 1)
##         stop("'column' should be a 'character' of length 1.")
##     if (any(colnames(x@spectraData) == column))
##         x@spectraData[, column]
##     else {
##         if (any(names(.SPECTRA_DATA_COLUMN_METHODS) == column))
##             do.call(.SPECTRA_DATA_COLUMN_METHODS[column], args = list(x))
##         else stop("No column '", column, "' available")
##     }
## }

#' Utility function to convert columns in the `x` `DataFrame` that have only
#' `NA`s to `Rle`. Also columns specified with parameter `columns` will be
#' converted (if present).
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
    all_na <- vapply(x, function(z) all(is.na(z)), logical(1))
    all_na <- names(all_na)[all_na]
    nrx <- nrow(x)
    for (col in all_na) {
        if (!is(x[[col]], "Rle"))
            x[[col]] <- Rle(x[[col]][1], nrx)
    }
    columns <- intersect(columns, colnames(x))
    for (col in columns) {
        if (!is(x[[col]], "Rle"))
            x[[col]] <- Rle(x[[col]])
    }    
    x
}

.uncompress_spectra_data <- function(x) {
    rle_col <- vapply(x, is, logical(1), "Rle")
    cols <- names(rle_col[intersect(colnames(x)[rle_col],
                                    names(.SPECTRA_DATA_COLUMNS))])
    for (col in cols) {
        x[[col]] <- as.vector(x[[col]])
    }
    x
}

#' Initialize the `spectraData` `DataFrame` adding eventually missing columns.
#'
#' @param x `DataFrame`
#'
#' @author Johannes Rainer
#' 
#' @noRd
#' 
#' @examples
#'
#' df <- DataFrame(fromFile = rep(1L, 10))
#' .initialize_spectra_data(df)
.initialize_spectra_data <- function(x) {
    nrx <- nrow(x)
    if (nrx) {
        mis_cols <- setdiff(names(.SPECTRA_DATA_COLUMNS),
                            c("mz", "intensity", colnames(x)))
        cols <- lapply(mis_cols, function(z) {
            Rle(as(NA, .SPECTRA_DATA_COLUMNS[z]), nrx)
        })
        names(cols) <- mis_cols
        x <- cbind(x, DataFrame(cols))
        if (!any(colnames(x) == "mz"))
            x$mz <- rep(list(numeric()), nrx)
        if (!any(colnames(x) == "intensity"))
            x$intensity <- rep(list(numeric()), nrx)
    } else {
        x <- DataFrame(lapply(.SPECTRA_DATA_COLUMNS, do.call, args = list()))        
        colnames(x) <- names(
            .SPECTRA_DATA_COLUMNS)[!(names(.SPECTRA_DATA_COLUMNS) %in%
                                     c("mz", "intensity"))]
        x$mz <- list()
        x$intensity <- list()
    }
    x
}

#' Helper function to check the data type even if `x` is an `Rle`.
#'
#' @param x
#'
#' @param class2 `character(1)` specifying the class for which should be tested.
#' 
#' @author Johannes Rainer
#' 
#' @noRd
.is_class <- function(x, class2) {
    if (is(x, "Rle"))
        is(x@values, class2)
    else is(x, class2)
}
