#' @include hidden_aliases.R
NULL

.valid_ms_backend_files <- function(x) {
    x <- x[!is.na(x)]
    if (length(x)) {
        if (anyDuplicated(x))
            return("Duplicated file names found.")
        if (!all(file.exists(x)))
            return(paste0("File(s) ", paste(x[!file.exists(x)], collapse = ", "),
                          " not found."))
    }
    NULL
}

.valid_ms_backend_mod_count <- function(x, y) {
    if (length(x) != length(y))
        "Different number of source files and modification counters."
    else
        NULL
}
