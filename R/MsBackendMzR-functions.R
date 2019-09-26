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
    ## Remove core spectra variables that contain only `NA`
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
    res <- asVectorDataFrame(x@spectraData[, sp_cols, drop = FALSE])
    any_mz <- any(columns == "mz")
    any_int <- any(columns == "intensity")
    if (any_mz || any_int) {
        pks <- peaks(x)
        if (any_mz)
            res$mz <- NumericList(lapply(pks, "[", , 1), compress = FALSE)
        if (any_int)
            res$intensity <- NumericList(lapply(pks, "[", , 2),
                                         compress = FALSE)
    }
    other_cols <- setdiff(columns, c(sp_cols, "mz", "intensity"))
    if (length(other_cols)) {
        other_res <- lapply(other_cols, .get_rle_column, x = x@spectraData)
        names(other_res) <- other_cols
        res <- cbind(res, as(other_res, "DataFrame"))
    }
    res[, columns, drop = FALSE]
}
