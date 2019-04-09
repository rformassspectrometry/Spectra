#' @include hidden_aliases.R
NULL

#' Helper function checking that files exist
#'
#' @param x `character` of file names.
#'
#' @return `character` or `NULL`
#'
#' @author Johannes Rainer
#'
#' @noRd
.valid_ms_backend_files_exist <- function(x) {
    if (!all(file.exists(x)))
        paste0("File(s) ", paste(x[!file.exists(x)], collapse = ", "),
               " not found.")
    else NULL
}

#' @rdname MsBackend
#'
#' @export MsBackendMzR
MsBackendMzR <- function() {
    if (!require("mzR", character.only = TRUE, quietly = TRUE))
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
    suppressPackageStartupMessages(
        require("mzR", quietly = TRUE, character.only = TRUE)
    )
    msd <- mzR::openMSfile(x)
    on.exit(mzR::close(msd))
    hdr <- mzR::header(msd)
    colnames(hdr)[colnames(hdr) == "seqNum"] <- "scanIndex"
    colnames(hdr)[colnames(hdr) == "precursorScanNum"] <- "precScanNum"
    colnames(hdr)[colnames(hdr) == "precursorMZ"] <- "precursorMz"
    S4Vectors::DataFrame(hdr)
}
