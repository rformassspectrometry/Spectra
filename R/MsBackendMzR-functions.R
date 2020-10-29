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
#' @importFrom MsCoreUtils vapply1l
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
    S4Vectors::DataFrame(hdr[, !(vapply1l(hdr, function(z) all(is.na(z))) &
                                 colnames(hdr) %in%
                                 names(.SPECTRA_DATA_COLUMNS))
                             ])
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
    res <- x@spectraData[, sp_cols, drop = FALSE]
    any_mz <- any(columns == "mz")
    any_int <- any(columns == "intensity")
    if (any_mz || any_int) {
        pks <- peaksData(x)
        if (any_mz)
            res$mz <- NumericList(lapply(pks, "[", , 1), compress = FALSE)
        if (any_int)
            res$intensity <- NumericList(lapply(pks, "[", , 2),
                                         compress = FALSE)
    }
    other_cols <- setdiff(columns, c(sp_cols, "mz", "intensity"))
    if (length(other_cols)) {
        other_res <- lapply(other_cols, .get_column, x = x@spectraData)
        names(other_res) <- other_cols
        res <- cbind(res, as(other_res, "DataFrame"))
    }
    res[, columns, drop = FALSE]
}

#' Function to write the MS data (for a single file/sample) to an mzML/mzXML
#' file.
#'
#' @param x `Spectra` object with the data
#'
#' @param file `character(1)` with the (full) output file name.
#'
#' @param format `character(1)` with the format.
#'
#' @param copy `logical(1)` whether general file information should be copied
#'     from the original MS file. Note this only works on `Spectra` with a
#'     `MsBackendMzR` backend as it assumes `dataStorage` contains the original
#'     MS data file.
#'
#' @noRd
#'
#' @author Johannes Rainer
.write_ms_data_mzR <- function(x, file, format = c("mzML", "mzXML"),
                               software_processing = NULL, copy = FALSE) {
    match.arg(tolower(format), c("mzml", "mzxml"))
    message("Writing file ", basename(file), "...", appendLF = FALSE)
    file <- path.expand(file)
    if (copy && length(unique(dataOrigin(x))) == 1) {
        original_file <- as.character(dataOrigin(x))[1]
        if (!file.exists(original_file)) {
            warning("Original data file not found, setting 'copy = FALSE'.",
                    call. = FALSE)
            copy <- FALSE
        }
    } else copy <- FALSE
    hdr <- as.data.frame(spectraData(x))
    pks <- as(x, "list")
    updated_values <- t(vapply(pks, function(z) {
        if (!nrow(z))
            return(c(peaksCount = 0, totIonCurrent = 0,
                     basePeakMZ = NA_real_, basePeakIntensity = NA_real_,
                     lowMZ = NA_real_, highMZ = NA_real_))
        z <- unname(z)
        max_peak <- which.max(z[, 2L])[1L]
        npk <- nrow(z)
        c(peaksCount = npk,
          totIonCurrent = sum(z[, 2L], na.rm = TRUE),
          basePeakMZ = z[max_peak, 1L],
          basePeakIntensity = z[max_peak, 2L],
          lowMZ = z[1L, 1L],
          highMZ = z[npk, 1L])
    }, numeric(6)))
    hdr$peaksCount <- updated_values[, "peaksCount"]
    hdr$totIonCurrent <- updated_values[, "totIonCurrent"]
    hdr$basePeakMZ <- updated_values[, "basePeakMZ"]
    hdr$basePeakIntensity <- updated_values[, "basePeakIntensity"]
    hdr$lowMZ <- updated_values[, "lowMZ"]
    hdr$highMZ <- updated_values[, "highMZ"]
    if (all(c("isolationWindowTargetMz", "isolationWindowLowerMz") %in%
            colnames(hdr))) {
        hdr$isolationWindowLowerOffset <- hdr$isolationWindowTargetMz -
            hdr$isolationWindowLowerMz
        hdr$isolationWindowLowerMz <- NULL
    }
    if (all(c("isolationWindowTargetMz", "isolationWindowUpperMz") %in%
            colnames(hdr))) {
        hdr$isolationWindowUpperOffset <- hdr$isolationWindowUpperMz -
            hdr$isolationWindowTargetMz
        hdr$isolationWindowUpperMz <- NULL
    }
    colnames(hdr)[colnames(hdr) == "rtime"] <- "retentionTime"
    colnames(hdr)[colnames(hdr) == "scanIndex"] <- "seqNum"
    colnames(hdr)[colnames(hdr) == "precursorMz"] <- "precursorMZ"
    colnames(hdr)[colnames(hdr) == "isolationWindowTargetMz"] <-
        "isolationWindowTargetMZ"
    req_cols <- c(acquisitionNum = "numeric",
                  msLevel = "integer",
                  polarity = "integer",
                  retentionTime = "numeric",
                  collisionEnergy = "numeric",
                  ionisationEnergy = "numeric",
                  precursorScanNum = "numeric",
                  precursorMZ = "numeric",
                  precursorCharge = "numeric",
                  precursorIntensity = "numeric",
                  mergedScan = "numeric",
                  mergedResultScanNum = "numeric",
                  mergedResultStartScanNum = "numeric",
                  mergedResultEndScanNum = "numeric",
                  injectionTime = "numeric",
                  filterString = "character",
                  centroided = "logical",
                  ionMobilityDriftTime = "numeric",
                  isolationWindowTargetMZ = "numeric",
                  isolationWindowLowerOffset = "numeric",
                  isolationWindowUpperOffset = "numeric",
                  scanWindowLowerLimit = "numeric",
                  scanWindowUpperLimit = "numeric"
                  )
    miss_cols <- req_cols[!(names(req_cols) %in% colnames(hdr))]
    if (length(miss_cols)) {
        nas <- rep(NA, nrow(hdr))
        miss_hdr <- do.call(cbind.data.frame,
                            lapply(miss_cols, function(z) {
                                as(nas, z)
                            }))
        hdr <- cbind(hdr, miss_hdr)
    }
    if (any(is.na(hdr$acquisitionNum)))
        hdr$acquisitionNum <- seq_len(nrow(hdr))
    ## seqNum has to be from 1...nrow
    hdr$seqNum <- seq_len(nrow(hdr))
    soft_proc <- .guess_software_processing(x)
    if (copy)
        mzR::copyWriteMSData(object = pks, file = file, header = hdr,
                             original_file = original_file, outformat = format,
                             software_processing = soft_proc)
    else
        mzR::writeMSData(object = pks, file = file, header = hdr,
                         outformat = format, software_processing = soft_proc)
    message("OK")
}

#' @description Determine data processing steps based on the `@processing`.
#'
#' @param x `Spectra`.
#'
#' @return `list` with the software processing(s).
#'
#' @author Johannes Rainer
#'
#' @importFrom utils packageVersion
#'
#' @noRd
.guess_software_processing <- function(x) {
    res <- vapply(x@processing, .pattern_to_cv, character(1))
    list(c("Spectra", as.character(packageVersion("Spectra")),
           "MS:-1", res[!is.na(res)]))
}

#' @description `.pattern_to_cv` performs a mapping of pattern to PSI-MS terms.
#'
#' @details The mapping bases on a manually curated list of pattern-CV-term
#'    pairs.
#'
#' @param pattern `character(1)` with the pattern for which a CV term should be
#'     returned.
#'
#' @param ifnotfound value to be returned if none of the cv terms
#'     matches the pattern.
#'
#' @param return `character` with the PSI-MS term or the value of `ifnotfound`
#'     if it was not found. Note that the length of the character can be > 1 if
#'     multiple terms match.
#'
#' @author Johannes Rainer
#'
#' @noRd
.pattern_to_cv <- function(pattern, ifnotfound = NA_character_) {
    ## character vector mapping pattern to PSI-MS term.
    .PATTERN.TO.CV <- c(
        filter = "MS:1001486",
        normali = "MS:1001484",
        calibration = "MS:1001485",
        pick = "MS:1000035",
        centroid = "MS:1000035",
        smooth = "MS:1000592",
        baseline = "MS:1000593",
        alignment = "MS:1000745"
    )
    ## Note: we are matching each names(.PATTERN.TO.CV) to `pattern`, not the
    ## other way round. So we can provide a string describing the processing
    ## step and return the best matching CV terms.
    matches <- vapply(names(.PATTERN.TO.CV), FUN = function(z) {
        any(grepl(pattern = z, x = pattern, ignore.case = TRUE))
    }, FUN.VALUE = logical(1))
    if (any(matches))
        unname(unique(.PATTERN.TO.CV[matches]))
    else
        ifnotfound
}
