#' @include hidden_aliases.R
NULL

#' @title In-memory MS data backend
#'
#' @name MsBackendDFrame
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @noRd
#'
#' @exportClass MsBackendDFrame
NULL

setClass("MsBackendDFrame",
         contains = "MsBackend",
         slots = c(spectraData = "DataFrame"),
         prototype = prototype(spectraData = DataFrame(),
                               readonly = FALSE,
                               version = "0.1"))

setValidity("MsBackendDFrame", function(object) {
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
setMethod("show", "MsBackendDFrame", function(object) {
    spd <- spectraData(object, c("msLevel", "rtime", "scanIndex"))
    cat(class(object), "with", nrow(spd), "spectra\n")
    if (nrow(spd)) {
        txt <- capture.output(print(spd))
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
setMethod("backendInitialize", signature = "MsBackendDFrame",
          function(object, spectraData, ...) {
              if (missing(spectraData)) spectraData <- DataFrame()
              if (is.data.frame(spectraData))
                  spectraData <- DataFrame(spectraData)
              if (!is(spectraData, "DataFrame"))
                  stop("'spectraData' has to be a 'DataFrame'")
              if (!nrow(spectraData))
                  return(object)
              spectraData$dataStorage <- "<memory>"
              spectraData <- asRleDataFrame(
                  spectraData, columns = c("dataStorage", "dataOrigin"))
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
setMethod("backendMerge", "MsBackendDFrame", function(object, ...) {
    object <- unname(c(object, ...))
    object <- object[lengths(object) > 0]
    res <- .combine_backend_dframe(object)
    validObject(res)
    res
})

## Data accessors

#' @rdname hidden_aliases
setMethod("acquisitionNum", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "acquisitionNum")
})

#' @rdname hidden_aliases
setMethod("as.list", "MsBackendDFrame", function(x) {
    mapply(cbind, mz = mz(x), intensity = intensity(x),
           SIMPLIFY = FALSE, USE.NAMES = FALSE)
})

#' @rdname hidden_aliases
setMethod("centroided", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "centroided")
})

#' @rdname hidden_aliases
#'
#' @aliases centroided<-,MsBackendDFrame-method
#'
#' @importFrom MsCoreUtils asRle
setReplaceMethod("centroided", "MsBackendDFrame", function(object, value) {
    value_len <- length(value)
    value_type <- is.logical(value)
    if (value_len == 1 && value_type)
        object@spectraData$centroided <- Rle(value, length(object))
    else if (value_len == length(object) && value_type)
        object@spectraData$centroided <- asRle(value)
    else
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("collisionEnergy", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "collisionEnergy")
})

#' @rdname hidden_aliases
setReplaceMethod("collisionEnergy", "MsBackendDFrame", function(object, value) {
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$collisionEnergy <- asRle(as.numeric(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("dataOrigin", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "dataOrigin")
})

#' @rdname hidden_aliases
setReplaceMethod("dataOrigin", "MsBackendDFrame", function(object, value) {
    if (!is.character(value) || length(value) != length(object))
        stop("'value' has to be a 'character' of length ", length(object))
    object@spectraData$dataOrigin <- asRle(as.character(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("dataStorage", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "dataStorage")
})

#' @rdname hidden_aliases
setReplaceMethod("dataStorage", "MsBackendDFrame", function(object, value) {
    if (!is.character(value) || length(value) != length(object))
        stop("'value' has to be a 'character' of length ", length(object))
    object@spectraData$dataStorage <- asRle(as.character(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("intensity", "MsBackendDFrame", function(object) {
    if (any(colnames(object@spectraData) == "intensity"))
        object@spectraData$intensity
    else {
        lst <- NumericList(numeric(), compress = FALSE)
        lst[rep(1, times = length(object))]
    }
})

#' @rdname hidden_aliases
setReplaceMethod("intensity", "MsBackendDFrame", function(object, value) {
    if (!(is.list(value) || inherits(value, "NumericList")))
        stop("'value' has to be a list or NumericList")
    if (length(value) != length(object))
        stop("length of 'value' has to match the length of 'object'")
    if (!all(lengths(value) == lengths(object)))
        stop("lengths of 'value' has to match the number of peaks ",
             "(i.e. lengths(object))")
    if (!is(value, "NumericList"))
        value <- NumericList(value, compress = FALSE)
    object@spectraData$intensity <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
#' @importFrom MsCoreUtils vapply1d
setMethod("ionCount", "MsBackendDFrame", function(object) {
    vapply1d(intensity(object), sum, na.rm = TRUE)
})

#' @rdname hidden_aliases
#' @importFrom MsCoreUtils vapply1l
setMethod("isCentroided", "MsBackendDFrame", function(object, ...) {
    vapply1l(as.list(object), .peaks_is_centroided)
})

#' @rdname hidden_aliases
setMethod("isEmpty", "MsBackendDFrame", function(x) {
    lengths(intensity(x)) == 0
})

#' @rdname hidden_aliases
setMethod("isolationWindowLowerMz", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "isolationWindowLowerMz")
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowLowerMz", "MsBackendDFrame",
                 function(object, value) {
                     if (!is.numeric(value) | length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$isolationWindowLowerMz <-
                         asRle(as.numeric(value))
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("isolationWindowTargetMz", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "isolationWindowTargetMz")
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowTargetMz", "MsBackendDFrame",
                 function(object, value) {
                     if (!is.numeric(value) || length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$isolationWindowTargetMz <-
                         asRle(as.numeric(value))
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("isolationWindowUpperMz", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "isolationWindowUpperMz")
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowUpperMz", "MsBackendDFrame",
                 function(object, value) {
                     if (!is.numeric(value) || length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$isolationWindowUpperMz <-
                         asRle(as.numeric(value))
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("length", "MsBackendDFrame", function(x) {
    nrow(x@spectraData)
})

#' @rdname hidden_aliases
setMethod("lengths", "MsBackendDFrame", function(x, use.names = FALSE) {
    lengths(mz(x))
})

#' @rdname hidden_aliases
setMethod("msLevel", "MsBackendDFrame", function(object, ...) {
    .get_rle_column(object@spectraData, "msLevel")
})

#' @rdname hidden_aliases
setMethod("mz", "MsBackendDFrame", function(object) {
    if (any(colnames(object@spectraData) == "mz"))
        object@spectraData$mz
    else {
        lst <- NumericList(numeric(), compress = FALSE)
        lst[rep(1, times = length(object))]
    }
})

#' @rdname hidden_aliases
setReplaceMethod("mz", "MsBackendDFrame", function(object, value) {
    if (!(is.list(value) || inherits(value, "NumericList")))
        stop("'value' has to be a list or NumericList")
    if (length(value) != length(object))
        stop("length of 'value' has to match the length of 'object'")
    if (!all(lengths(value) == lengths(object)))
        stop("lengths of 'value' has to match the number of peaks ",
             "(i.e. lengths(object))")
    if (!is(value, "NumericList"))
        value <- NumericList(value, compress = FALSE)
    object@spectraData$mz <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("polarity", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "polarity")
})

#' @rdname hidden_aliases
setReplaceMethod("polarity", "MsBackendDFrame", function(object, value) {
    value_len <- length(value)
    value_type <- is.numeric(value)
    if (value_len == 1 && value_type)
        object@spectraData$polarity <- Rle(as.integer(value), length(object))
    else if (value_len == length(object) && value_type)
        object@spectraData$polarity <- asRle(as.integer(value))
    else
        stop("'value' has to be an 'integer' of length 1 or ", length(object))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("precScanNum", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "precScanNum")
})

#' @rdname hidden_aliases
setMethod("precursorCharge", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "precursorCharge")
})

#' @rdname hidden_aliases
setMethod("precursorIntensity", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "precursorIntensity")
})

#' @rdname hidden_aliases
setMethod("precursorMz", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "precursorMz")
})

#' @rdname hidden_aliases
setReplaceMethod("replaceList", "MsBackendDFrame", function(object, value) {
    if (!(is.list(value) || inherits(value, "SimpleList")))
        stop("'value' has to be a list-like object")
    if (length(value) != length(object))
        stop("Length of 'value' has to match length of 'object'")
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
setMethod("rtime", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "rtime")
})

#' @rdname hidden_aliases
setReplaceMethod("rtime", "MsBackendDFrame", function(object, value) {
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$rtime <- asRle(as.numeric(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("scanIndex", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "scanIndex")
})

#' @rdname hidden_aliases
setMethod("selectSpectraVariables", "MsBackendDFrame",
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
setMethod("smoothed", "MsBackendDFrame", function(object) {
    .get_rle_column(object@spectraData, "smoothed")
})

#' @rdname hidden_aliases
#'
#' @aliases smoothed<-,MsBackendDFrame-method
setReplaceMethod("smoothed", "MsBackendDFrame", function(object, value) {
    value_len <- length(value)
    value_type <- is.logical(value)
    if (value_len == 1 && value_type)
        object@spectraData$smoothed <- Rle(value, length(object))
    else if (value_len == length(object) && value_type)
        object@spectraData$smoothed <- asRle(value)
    else
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    validObject(object)
    object
})

#' @rdname hidden_aliases
#'
#' @importFrom methods as
#'
#' @importFrom S4Vectors SimpleList
#'
#' @importFrom MsCoreUtils asVectorDataFrame
#'
#' @importMethodsFrom S4Vectors lapply
setMethod("spectraData", "MsBackendDFrame",
          function(object, columns = spectraVariables(object)) {
              df_columns <- intersect(columns,colnames(object@spectraData))
              res <- object@spectraData[, df_columns, drop = FALSE]
              other_columns <- setdiff(columns,colnames(object@spectraData))
              if (length(other_columns)) {
                  other_res <- lapply(other_columns, .get_rle_column,
                                      x = object@spectraData)
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
              asVectorDataFrame(res[, columns, drop = FALSE])
          })

#' @rdname hidden_aliases
setReplaceMethod("spectraData", "MsBackendDFrame", function(object, value) {
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
    object@spectraData <- asRleDataFrame(
        value, columns = c("dataStorage", "dataOrigin"))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("spectraNames", "MsBackendDFrame", function(object) {
    rownames(object@spectraData)
})

#' @rdname hidden_aliases
setReplaceMethod("spectraNames", "MsBackendDFrame", function(object, value) {
    rownames(object@spectraData) <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("spectraVariables", "MsBackendDFrame", function(object) {
    unique(c(names(.SPECTRA_DATA_COLUMNS), colnames(object@spectraData)))
})

#' @rdname hidden_aliases
setMethod("tic", "MsBackendDFrame", function(object, initial = TRUE) {
    if (initial) {
        if (any(colnames(object@spectraData) == "totIonCurrent"))
            .get_rle_column(object@spectraData, "totIonCurrent")
        else rep(NA_real_, times = length(object))
    } else vapply1d(intensity(object), sum, na.rm = TRUE)
})

#' @rdname hidden_aliases
setMethod("$", "MsBackendDFrame", function(x, name) {
    if (!any(spectraVariables(x) == name))
        stop("spectra variable '", name, "' not available")
    spectraData(x, name)[, 1]
})

#' @rdname hidden_aliases
setReplaceMethod("$", "MsBackendDFrame", function(x, name, value) {
    if (is.list(value) && any(c("mz", "intensity") == name))
        value <- NumericList(value, compress = FALSE)
    if (name == "dataStorage")
        value <- Rle(value)
    x@spectraData[[name]] <- asRle(value)
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
#' @importFrom MsCoreUtils i2index
#'
#' @rdname hidden_aliases
setMethod("[", "MsBackendDFrame", function(x, i, j, ..., drop = FALSE) {
    .subset_backend_data_frame(x, i)
})

#' @rdname hidden_aliases
setMethod("split", "MsBackendDFrame", function(x, f, drop = FALSE, ...) {
    if (!is.factor(f))
        f <- as.factor(f)
    if (length(levels(f)) > (nrow(x@spectraData) / 10))
        slot(x, "spectraData", check = FALSE) <-
            asVectorDataFrame(x@spectraData)
    lapply(split(seq_along(x), f, ...), function(i) x[i, ])
})

#' @rdname hidden_aliases
setMethod("filterAcquisitionNum", "MsBackendDFrame",
          function(object, n = integer(), dataStorage = character(),
                   dataOrigin = character()) {
    if (!length(n) || !length(object)) return(object)
    if (!is.integer(n)) stop("'n' has to be an integer representing the ",
                             "acquisition number(s) for sub-setting")
    sel_file <- .sel_file(object, dataStorage, dataOrigin)
    sel_acq <- acquisitionNum(object) %in% n & sel_file
    object[sel_acq | !sel_file]
})

#' @rdname hidden_aliases
setMethod("filterDataOrigin", "MsBackendDFrame",
          function(object, dataOrigin = character()) {
              if (length(dataOrigin)) {
                  object <- object[dataOrigin(object) %in% dataOrigin]
                  if (is.unsorted(dataOrigin))
                      object[order(match(dataOrigin(object), dataOrigin))]
                  else object
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterDataStorage", "MsBackendDFrame",
          function(object, dataStorage = character()) {
              if (length(dataStorage)) {
                  object <- object[dataStorage(object) %in% dataStorage]
                  if (is.unsorted(dataStorage))
                      object[order(match(dataStorage(object), dataStorage))]
                  else object
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterEmptySpectra", "MsBackendDFrame", function(object) {
    if (!length(object)) return(object)
    object[as.logical(lengths(object))]
})

#' @rdname hidden_aliases
setMethod("filterIsolationWindow", "MsBackendDFrame", {
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
setMethod("filterMsLevel", "MsBackendDFrame",
          function(object, msLevel = integer()) {
              if (length(msLevel)) {
                  object[msLevel(object) %in% msLevel]
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterPolarity", "MsBackendDFrame",
          function(object, polarity = integer()) {
              if (length(polarity))
                  object[polarity(object) %in% polarity]
              else object
          })

#' @rdname hidden_aliases
setMethod("filterPrecursorMz", "MsBackendDFrame",
          function(object, mz = numeric()) {
              if (length(mz)) {
                  mz <- range(mz)
                  keep <- which(precursorMz(object) >= mz[1] &
                                precursorMz(object) <= mz[2])
                  object[keep]
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterPrecursorScan", "MsBackendDFrame",
          function(object, acquisitionNum = integer()) {
              if (length(acquisitionNum)) {
                  object[.filterSpectraHierarchy(acquisitionNum(object),
                                                 precScanNum(object),
                                                 acquisitionNum)]
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterRt", "MsBackendDFrame",
          function(object, rt = numeric(), msLevel. = unique(msLevel(object))) {
              if (length(rt)) {
                  rt <- range(rt)
                  sel_ms <- msLevel(object) %in% msLevel.
                  sel_rt <- rtime(object) >= rt[1] &
                      rtime(object) <= rt[2] & sel_ms
                  object[sel_rt | !sel_ms]
              } else object
          })
