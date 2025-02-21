#' @include hidden_aliases.R
NULL

.valid_ms_backend_data_storage <- function(x) {
    if (anyNA(x))
        return("'NA' values in dataStorage are not allowed.")
    NULL
}

.valid_ms_backend_files_exist <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) && !all(file.exists(x)))
        return(paste0("File(s) ", paste(x[!file.exists(x)], collapse = ", "),
                      " not found"))
    NULL
}

#' @title Fill spectra data  with columns for missing core variables
#'
#' @description
#'
#' `fillCoreSpectraVariables()` fills a provided `data.frame`
#' with columns for eventually missing *core* spectra variables.
#' The missing core variables are added as new columns with missing values
#' (`NA`) of the correct data type.
#' Use [coreSpectraVariables()] to list the set of core variables and their
#' data types.
#'
#' @param x `data.frame` or `DataFrame` with potentially present core
#'     variable columns.
#'
#' @param columns `character` with the names of the (core) spectra variables
#'     that should be added if not already present in `x`. Defaults to
#'     `columns = names(coreSpectraVariables())`.
#'
#' @return input data frame `x` with missing core variables added (with the
#'     correct data type).
#'
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#'
#' ## Define a data frame
#' a <- data.frame(msLevel = c(1L, 1L, 2L), other_column = "b")
#'
#' ## Add missing core chromatogram variables to this data frame
#' fillCoreSpectraVariables(a)
#'
#' ## The data.frame thus contains columns for all core spectra
#' ## variables in the respective expected data type (but filled with
#' ## missing values).
fillCoreSpectraVariables <- function(x = data.frame(),
                                     columns = names(coreSpectraVariables())) {
    nr <- nrow(x)
    cv <- .SPECTRA_DATA_COLUMNS[names(.SPECTRA_DATA_COLUMNS) %in% columns]
    miss <- cv[setdiff(names(cv), c(colnames(x), c("mz", "intensity")))]
    if (length(miss))
        x <- cbind(x, lapply(miss, function(z, n) rep(as(NA, z), n), nr))
    if (any(columns == "mz")) {
        a <- numeric()
        x$mz <- NumericList(replicate(nr, a), compress = FALSE)
    }
    if (any(columns == "intensity")) {
        a <- numeric()
        x$intensity <- NumericList(replicate(nr, a), compress = FALSE)
    }
    x
}
