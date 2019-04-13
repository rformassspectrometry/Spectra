#' @include hidden_aliases.R
NULL

.valid_processing_queue <- function(x) {
    if (length(x))
        if (!all(vapply(x, inherits, logical(1), "ProcessingStep")))
            return("'processingQueue' should only contain ProcessingStep objects.")
    NULL
}

#' @description
#'
#' Combine two `DataFrame`s. The resulting `DataFrame` contains all
#' columns from both `x` and `y`. For columns present in both `DataFrame`s those
#' in `x` will be used. Also, the resulting `DataFrame` uses the row names of
#' `x` unless `x` has no row names.
#'
#' @param x `DataFrame`
#'
#' @param y `DataFrame`
#'
#' @return `DataFrame` with the merged columns.
#'
#' @author Johannes Rainer
#'
#' @importMethodsFrom S4Vectors nrow rownames colnames cbind
#'
#' @noRd
.combine_data_frame <- function(x, y) {
    if (nrow(y) == 0)
        return(x)
    if (nrow(x) == 0)
        return(y)
    if (nrow(x) != nrow(y))
        stop("'x' and 'y' have to have the same number of rows")
    if (is.null(rownames(x)) & !is.null(rownames(y)))
        rownames(x) <- rownames(y)
    cols_y <- !(colnames(y) %in% colnames(x))
    if (any(cols_y))
        x <- cbind(x, y[, cols_y, drop = FALSE])
    x
}

#' @importMethodsFrom S4Vectors $ $<-
#'
#' @rdname Spectra
#'
#' @export
Spectra <- function(spectraData = DataFrame(), msLevel, rtime, acquisitionNum,
                    scanIndex, fromFile, centroided, smoothed, polarity,
                    precScanNum, precursorMz, precursorIntensity,
                    precursorCharge, collisionEnergy, mz, intensity,
                    backend = MsBackendDataFrame(), processingQueue = list(),
                    metadata = list()) {
    spData <- list()
    if (!missing(msLevel))
        spData <- c(spData, list(msLevel = .as_integer(msLevel)))
    if (!missing(rtime))
        spData <- c(spData, list(rtime = rtime))
    if (!missing(acquisitionNum))
        spData <- c(spData, list(acquisitionNum = .as_integer(acquisitionNum)))
    if (!missing(scanIndex))
        spData <- c(spData, list(scanIndex = .as_integer(scanIndex)))
    if (!missing(fromFile))
        spData <- c(spData, list(fromFile = .as_integer(fromFile)))
    if (!missing(centroided))
        spData <- c(spData, list(centroided = as.logical(centroided)))
    if (!missing(smoothed))
        spData <- c(spData, list(smoothed = as.logical(smoothed)))
    if (!missing(polarity))
        spData <- c(spData, list(polarity = .as_integer(polarity)))
    if (!missing(precScanNum))
        spData <- c(spData, list(precScanNum = .as_integer(precScanNum)))
    if (!missing(precursorMz))
        spData <- c(spData, list(precursorMz = precursorMz))
    if (!missing(precursorIntensity))
        spData <- c(spData, list(precursorIntensity = precursorIntensity))
    if (!missing(precursorCharge))
        spData <- c(spData, list(precursorCharge = .as_integer(precursorCharge)))
    if (!missing(collisionEnergy))
        spData <- c(spData, list(collisionEnergy = collisionEnergy))
    ## Check lengths and build spectraData
    if (length(unique(lengths(spData))) > 1)
        stop("Different lengths of parameters. Please ensure that the length",
             " of all passed parameters matches the number of spectra.")
    spData <- DataFrame(spData)
    spectraData <- .combine_data_frame(spData, spectraData)
    if (!missing(mz)) {
        if (nrow(spectraData) && length(mz) != nrow(spectraData))
            stop("Length of 'mz' does not match the length of other parameters")
        if (!nrow(spectraData))
            spectraData <- DataFrame(msLevel = rep(NA_integer_, length(mz)))
        spectraData$mz <- mz
    }
    if (!missing(intensity)) {
        if (nrow(spectraData) && length(intensity) != nrow(spectraData))
            stop("Length of 'intensity' does not match the length of other",
                 " parameters")
        if (!nrow(spectraData))
            spectraData <- DataFrame(msLevel = rep(NA_integer_, length(mz)))
        spectraData$intensity <- intensity
    }
    if (nrow(spectraData))
        spectraData$fromFile <- 1L
    fl <- if (nrow(spectraData)) NA_character_ else character()
    backend <- backendInitialize(backend, files = fl, spectraData)
    new("Spectra", backend = backend, processingQueue = processingQueue,
        metadata = metadata)
}
