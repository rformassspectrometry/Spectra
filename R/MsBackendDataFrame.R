#' @include hidden_aliases.R
NULL

#' @title In-memory MS data backend
#'
#' @description
#'
#' @name MsBackendDataFrame
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @noRd
#'
#' @exportClass MsBackendDataFrame
NULL

setClass("MsBackendDataFrame",
         contains = "MsBackend",
         slots = c(spectraData = "DataFrame"),
         prototype = prototype(spectraData = DataFrame(),
                               readonly = FALSE,
                               version = "0.1"))

setValidity("MsBackendDataFrame", function(object) {
    msg <- .valid_spectra_data_required_columns(object@spectraData)
    if (length(msg))
        return(msg)
    msg <- c(
        .valid_column_datatype(object@spectraData, .SPECTRA_DATA_COLUMNS),
        .valid_intensity_column(object@spectraData),
        .valid_mz_column(object@spectraData),
        .valid_intensity_mz_columns(object@spectraData))
    if (is.null(msg)) TRUE
    else msg
})

#' @rdname hidden_aliases
setMethod("show", "MsBackendDataFrame", function(object) {
    spd <- spectraData(object, c("msLevel", "rtime", "scanIndex"))
    cat(class(object), "with", nrow(spd), "spectra\n")
    if (nrow(spd)) {
        txt <- capture.output(show(spd))
        cat(txt[-1], sep = "\n")
        sp_cols <- spectraVariables(object)
        cat(" ...", length(sp_cols) - 3, "more variables/columns.\n")
    }
})

#' @importMethodsFrom S4Vectors $ $<-
#'
#' @importClassesFrom IRanges NumericList
#'
#' @importFrom IRanges NumericList
#'
#' @rdname hidden_aliases
setMethod("backendInitialize", signature = "MsBackendDataFrame",
          function(object, spectraData, ...) {
              if (missing(spectraData)) spectraData <- DataFrame()
              if (is.data.frame(spectraData))
                  spectraData <- DataFrame(spectraData)
              if (!is(spectraData, "DataFrame"))
                  stop("'spectraData' has to be a 'DataFrame'")
              if (nrow(spectraData) && is.null(spectraData$dataStorage))
                  spectraData$dataStorage <- "<memory>"
              spectraData <- .as_rle_spectra_data(spectraData)
              if (nrow(spectraData) && !is(spectraData$mz, "NumericList"))
                  spectraData$mz <- NumericList(spectraData$mz, compress = FALSE)
              if (nrow(spectraData) && !is(spectraData$intensity, "NumericList"))
                  spectraData$intensity <- NumericList(spectraData$intensity,
                                                       compress = FALSE)
              object@spectraData <- spectraData
              validObject(object)
              object
          })

#' @rdname hidden_aliases
setMethod("backendMerge", "MsBackendDataFrame", function(object, ...) {
    object <- unname(c(object, ...))
    object <- object[lengths(object) > 0]
    res <- .combine_backend_data_frame(object)
    validObject(res)
    res
})

## Data accessors

#' @rdname hidden_aliases
setMethod("acquisitionNum", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "acquisitionNum")
})

#' @rdname hidden_aliases
setMethod("centroided", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "centroided")
})

#' @rdname hidden_aliases
setReplaceMethod("centroided", "MsBackendDataFrame", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.logical(value) | length(value) != length(object))
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    object@spectraData$centroided <- .as_rle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("collisionEnergy", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "collisionEnergy")
})

#' @rdname hidden_aliases
setReplaceMethod("collisionEnergy", "MsBackendDataFrame", function(object, value) {
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$collisionEnergy <- .as_rle(as.numeric(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("dataOrigin", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "dataOrigin")
})

#' @rdname hidden_aliases
setReplaceMethod("dataOrigin", "MsBackendDataFrame", function(object, value) {
    if (!is.character(value) | length(value) != length(object))
        stop("'value' has to be a 'character' of length ", length(object))
    object@spectraData$dataOrigin <- .as_rle(as.character(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("dataStorage", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "dataStorage")
})

#' @rdname hidden_aliases
setReplaceMethod("dataStorage", "MsBackendDataFrame", function(object, value) {
    if (!is.character(value) | length(value) != length(object))
        stop("'value' has to be a 'character' of length ", length(object))
    object@spectraData$dataStorage <- .as_rle(as.character(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("dataStorageNames", "MsBackendDataFrame", function(object) {
    if (is(object@spectraData$dataStorage, "Rle"))
        unique(object@spectraData$dataStorage@values)
    else unique(dataStorage(object))
})

#' @rdname hidden_aliases
setMethod("intensity", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "intensity"))
        object@spectraData$intensity
    else {
        lst <- NumericList(numeric(), compress = FALSE)
        lst[rep(1, times = length(object))]
    }
})

#' @rdname hidden_aliases
setReplaceMethod("intensity", "MsBackendDataFrame", function(object, value) {
    if (!(is.list(value) || inherits(value, "NumericList")))
        stop("'value' has to be a list or NumericList")
    if (length(value) != length(object))
        stop("length of 'value' has to match the length of 'object'")
    if (!all(lengths(value) == peaksCount(object)))
        stop("lengths of 'value' has to match the number of peaks ",
             "(i.e. peaksCount(object))")
    if (!is(value, "NumericList"))
        value <- NumericList(value, compress = FALSE)
    object@spectraData$intensity <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("ionCount", "MsBackendDataFrame", function(object) {
    vapply(intensity(object), sum, numeric(1), na.rm = TRUE)
})

#' @rdname hidden_aliases
setMethod("isCentroided", "MsBackendDataFrame", function(object, ...) {
    vapply(peaks(object), .isCentroided, logical(1))
})

#' @rdname hidden_aliases
setMethod("isEmpty", "MsBackendDataFrame", function(x) {
    lengths(intensity(x)) == 0
})

#' @rdname hidden_aliases
setMethod("isolationWindowLowerMz", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "isolationWindowLowerMz")
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowLowerMz", "MsBackendDataFrame",
                 function(object, value) {
                     if (!is.numeric(value) | length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$isolationWindowLowerMz <-
                         .as_rle(as.numeric(value))
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("isolationWindowTargetMz", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "isolationWindowTargetMz")
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowTargetMz", "MsBackendDataFrame",
                 function(object, value) {
                     if (!is.numeric(value) || length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$isolationWindowTargetMz <-
                         .as_rle(as.numeric(value))
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("isolationWindowUpperMz", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "isolationWindowUpperMz")
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowUpperMz", "MsBackendDataFrame",
                 function(object, value) {
                     if (!is.numeric(value) || length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$isolationWindowUpperMz <-
                         .as_rle(as.numeric(value))
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("length", "MsBackendDataFrame", function(x) {
    nrow(x@spectraData)
})

#' @rdname hidden_aliases
setMethod("msLevel", "MsBackendDataFrame", function(object, ...) {
    .get_rle_column(object@spectraData, "msLevel")
})

#' @rdname hidden_aliases
setMethod("mz", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "mz"))
        object@spectraData$mz
    else {
        lst <- NumericList(numeric(), compress = FALSE)
        lst[rep(1, times = length(object))]
    }
})

#' @rdname hidden_aliases
setReplaceMethod("mz", "MsBackendDataFrame", function(object, value) {
    if (!(is.list(value) || inherits(value, "NumericList")))
        stop("'value' has to be a list or NumericList")
    if (length(value) != length(object))
        stop("length of 'value' has to match the length of 'object'")
    if (!all(lengths(value) == peaksCount(object)))
        stop("lengths of 'value' has to match the number of peaks ",
             "(i.e. peaksCount(object))")
    if (!is(value, "NumericList"))
        value <- NumericList(value, compress = FALSE)
    object@spectraData$mz <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("peaks", "MsBackendDataFrame", function(object) {
    mapply(mz(object), intensity(object), FUN = function(m, i)
        cbind(mz = m, intensity = i), SIMPLIFY = FALSE, USE.NAMES = FALSE)
})

#' @rdname hidden_aliases
setReplaceMethod("peaks", "MsBackendDataFrame", function(object, value) {
    if (!(is.list(value) || inherits(value, "SimpleList")))
        stop("'value' has to be a list-like object")
    if (length(value) != length(object))
        stop("Length of 'value' has to match length of 'object'")
    object@modCount <- object@modCount + 1L
    vals <- lapply(value, "[", , 1L)
    if (!is(vals, "NumericList"))
        vals <- NumericList(vals, compress = FALSE)
    object@spectraData$mz <- vals
    vals <- lapply(value, "[", , 2L)
    if (!is(vals, "NumericList"))
        vals <- NumericList(vals, compress = FALSE)
    object@spectraData$intensity <- vals
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("peaksCount", "MsBackendDataFrame", function(object) {
    lengths(mz(object))
})

#' @rdname hidden_aliases
setMethod("polarity", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "polarity")
})

#' @rdname hidden_aliases
setReplaceMethod("polarity", "MsBackendDataFrame", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be an 'integer' of length 1 or ", length(object))
    object@spectraData$polarity <- .as_rle(as.integer(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("precScanNum", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "precScanNum")
})

#' @rdname hidden_aliases
setMethod("precursorCharge", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "precursorCharge")
})

#' @rdname hidden_aliases
setMethod("precursorIntensity", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "precursorIntensity")
})

#' @rdname hidden_aliases
setMethod("precursorMz", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "precursorMz")
})

#' @rdname hidden_aliases
setMethod("rtime", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "rtime")
})

#' @rdname hidden_aliases
setReplaceMethod("rtime", "MsBackendDataFrame", function(object, value) {
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$rtime <- .as_rle(as.numeric(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("scanIndex", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "scanIndex")
})

#' @rdname hidden_aliases
setMethod("selectSpectraVariables", "MsBackendDataFrame",
          function(object, spectraVariables = spectraVariables(object)) {
              if (!all(spectraVariables %in% spectraVariables(object)))
                  stop("Spectra variables ",
                       paste(spectraVariables[!(spectraVariables %in%
                                                spectraVariables(object))],
                             collapse = ", "), " not available")
              keep <- spectraVariables[spectraVariables %in%
                                            colnames(object@spectraData)]
              if (length(keep))
                  object@spectraData <- object@spectraData[, keep,
                                                           drop = FALSE]
              msg <- .valid_spectra_data_required_columns(object@spectraData)
              if (length(msg))
                  stop(msg)
              validObject(object)
              object
})

#' @rdname hidden_aliases
setMethod("smoothed", "MsBackendDataFrame", function(object) {
    .get_rle_column(object@spectraData, "smoothed")
})

#' @rdname hidden_aliases
setReplaceMethod("smoothed", "MsBackendDataFrame", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.logical(value) || length(value) != length(object))
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    object@spectraData$smoothed <- .as_rle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
#'
#' @importFrom methods as
#'
#' @importFrom S4Vectors SimpleList
#'
#' @importMethodsFrom S4Vectors lapply
setMethod("spectraData", "MsBackendDataFrame",
          function(object, columns = spectraVariables(object)) {
              df_columns <- intersect(columns,colnames(object@spectraData))
              res <- object@spectraData[, df_columns, drop = FALSE]
              other_columns <- setdiff(columns,colnames(object@spectraData))
              if (length(other_columns)) {
                  other_res <- lapply(other_columns, .get_spectra_data_column,
                                      x = object)
                  names(other_res) <- other_columns
                  is_mz_int <- names(other_res) %in% c("mz", "intensity")
                  if (!all(is_mz_int))
                      res <- cbind(res, as(other_res[!is_mz_int], "DataFrame"))
                  if (any(names(other_res) == "mz"))
                      res$mz <- if (length(other_res$mz)) other_res$mz
                                else NumericList(compress = FALSE)
                  if (any(names(other_res) == "intensity"))
                      res$intensity <- if (length(other_res$intensity)) other_res$intensity
                                       else NumericList(compress = FALSE)
              }
              .as_vector_spectra_data(res[, columns, drop = FALSE])
          })

#' @rdname hidden_aliases
setReplaceMethod("spectraData", "MsBackendDataFrame", function(object, value) {
    if (inherits(value, "DataFrame")) {
        if (length(object) && nrow(value) != length(object))
            stop("'value' has to be a 'DataFrame' with ", length(object), " rows.")
        if (!is(value$mz, "NumericList"))
            value$mz <- NumericList(value$mz, compress = FALSE)
        if (!is(value$intensity, "NumericList"))
            value$intensity <- NumericList(value$intensity, compress = FALSE)
        if (is.null(value$dataStorage))
            value$dataStorage <- "<memory>"
    } else {
        if (length(value) == 1)
            value <- rep(value, length(object))
        if (length(value) != length(object))
            stop("length of 'value' has to be ", length(object))
    }
    object@spectraData <- .as_rle_spectra_data(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("spectraNames", "MsBackendDataFrame", function(object) {
    rownames(object@spectraData)
})

#' @rdname hidden_aliases
setReplaceMethod("spectraNames", "MsBackendDataFrame", function(object, value) {
    rownames(object@spectraData) <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("spectraVariables", "MsBackendDataFrame", function(object) {
    unique(c(names(.SPECTRA_DATA_COLUMNS), colnames(object@spectraData)))
})

#' @rdname hidden_aliases
setMethod("tic", "MsBackendDataFrame", function(object, initial = TRUE) {
    if (initial) {
        if (any(colnames(object@spectraData) == "totIonCurrent"))
            .get_rle_column(object@spectraData, "totIonCurrent")
        else rep(NA_real_, times = length(object))
    } else vapply(intensity(object), sum, numeric(1), na.rm = TRUE)
})

#' @rdname hidden_aliases
setMethod("$", "MsBackendDataFrame", function(x, name) {
    if (!any(spectraVariables(x) == name))
        stop("spectra variable '", name, "' not available")
    spectraData(x, name)[, 1]
})

#' @rdname hidden_aliases
setReplaceMethod("$", "MsBackendDataFrame", function(x, name, value) {
    if (is.list(value) && any(c("mz", "intensity") == name))
        value <- NumericList(value, compress = FALSE)
    if (name == "dataStorage")
        value <- Rle(value)
    x@spectraData[[name]] <- .as_rle(value)
    validObject(x)
    x
})

#### ---------------------------------------------------------------------------
##
##                      FILTERING AND SUBSETTING
##
#### ---------------------------------------------------------------------------

#' @importMethodsFrom S4Vectors [
#'
#' @rdname hidden_aliases
setMethod("[", "MsBackendDataFrame", function(x, i, j, ..., drop = FALSE) {
    if (!missing(j))
        stop("Subsetting by column ('j = ", j, "' is not supported")
    i <- .i_to_index(i, length(x), rownames(x@spectraData))
    x@spectraData <- x@spectraData[i, , drop = FALSE]
    validObject(x)
    x
})

#' @rdname hidden_aliases
setMethod("filterAcquisitionNum", "MsBackendDataFrame",
          function(object, n = integer(), dataStorage = integer(),
                   dataOrigin = integer()) {
    if (!length(n) | !length(object)) return(object)
    if (!is.integer(n)) stop("'n' has to be an integer representing the ",
                             "acquisition number(s) for sub-setting")
    sel_file <- .sel_file(object, dataStorage, dataOrigin)
    sel_acq <- acquisitionNum(object) %in% n & sel_file
    object[sel_acq | !sel_file]
})

#' @rdname hidden_aliases
setMethod("filterDataOrigin", "MsBackendDataFrame",
          function(object, dataOrigin = integer()) {
              if (length(dataOrigin)) {
                  lvls <- unique(dataOrigin(object))
                  dataOrigin <- .i_to_index(dataOrigin, length(lvls), lvls)
                  object <- object[dataOrigin(object) %in% lvls[dataOrigin]]
                  if (is.unsorted(dataOrigin))
                      object[order(match(dataOrigin(object),
                                         lvls[dataOrigin]))]
                  else object
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterDataStorage", "MsBackendDataFrame",
          function(object, dataStorage = integer()) {
              if (length(dataStorage)) {
                  lvls <- dataStorageNames(object)
                  dataStorage <- .i_to_index(dataStorage, length(lvls), lvls)
                  object <- object[dataStorage(object) %in% lvls[dataStorage]]
                  if (is.unsorted(dataStorage))
                      object[order(match(dataStorage(object),
                                         lvls[dataStorage]))]
                  else object
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterEmptySpectra", "MsBackendDataFrame", function(object) {
    if (!length(object)) return(object)
    object[as.logical(peaksCount(object))]
})

#' @rdname hidden_aliases
setMethod("filterIsolationWindow", "MsBackendDataFrame", {
    function(object, mz = numeric(), ...) {
        if (length(mz)) {
            if (length(mz) > 1)
                stop("'mz' is expected to be a single m/z value", call. = FALSE)
            keep <- which(isolationWindowLowerMz(object) <= mz &
                          isolationWindowUpperMz(object) >= mz)
            object[keep]
        } else object
    }
})

#' @rdname hidden_aliases
setMethod("filterMsLevel", "MsBackendDataFrame",
          function(object, msLevel = integer()) {
              if (length(msLevel)) {
                  object[msLevel(object) %in% msLevel]
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterPolarity", "MsBackendDataFrame",
          function(object, polarity = integer()) {
              if (length(polarity))
                  object[polarity(object) %in% polarity]
              else object
          })

#' @rdname hidden_aliases
setMethod("filterPrecursorMz", "MsBackendDataFrame",
          function(object, mz = numeric(), ppm = 0) {
              if (length(mz)) {
                  if (length(mz) > 1)
                      stop("'mz' is expected to be a single m/z value",
                           call. = FALSE)
                  mz_ppm <- ppm * mz / 1e6
                  keep <- which(precursorMz(object) >= (mz - mz_ppm) &
                                precursorMz(object) <= (mz + mz_ppm))
                  object[keep]
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterPrecursorScan", "MsBackendDataFrame",
          function(object, acquisitionNum = integer()) {
              if (length(acquisitionNum)) {
                  object[.filterSpectraHierarchy(acquisitionNum(object),
                                                 precScanNum(object),
                                                 acquisitionNum)]
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterRt", "MsBackendDataFrame",
          function(object, rt = numeric(), msLevel = integer()) {
              if (length(rt)) {
                  rt <- range(rt)
                  if (!is.numeric(rt))
                      stop("'rt' must be a numeric")
                  if (!length(msLevel))
                      msLevel <- unique(msLevel(object))
                  sel_ms <- msLevel(object) %in% msLevel
                  sel_rt <- rtime(object) >= rt[1] &
                      rtime(object) <= rt[2] & sel_ms
                  object[sel_rt | !sel_ms]
              } else object
          })
