#' @include hidden_aliases.R
NULL

.valid_ms_backend_data_storage <- function(x) {
    if (any(is.na(x)))
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
