#' @include hidden_aliases.R
NULL

#' @title Improved in-memory MS data backend
#'
#' @name MsBackendDF
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @noRd
#'
#' @exportClass MsBackendDF
NULL

#' - spectraData: `data.frame` with arbitrary spectra annotations.
#' - peaksData: `list` of `numeric` matrices with m/z and intensity values.
#' - peaksDataFrame: `list` of `data.frame`s with potentially arbitrary peak
#'     annotations.
#'
#' @noRd
setClass("MsBackendDF",
         contains = "MsBackend",
         slots = c(spectraData = "data.frame",
                   peaksData = "list",
                   peaksDataFrame = "list"),
         prototype = prototype(spectraData = data.frame(),
                               peaksData = list(),
                               peaksDataFrame = list(),
                               readonly = FALSE,
                               version = "0.1"))

setValidity("MsBackendDF", function(object) {
    msg <- .valid_spectra_data_required_columns(object@spectraData)
    if (length(msg))
        return(msg)
    msg <-
        .valid_column_datatype(object@spectraData, .SPECTRA_DATA_COLUMNS)
    if (length(object@peaksDataFrame) && length(object@peaksData) !=
        length(object@peaksDataFrame))
        msg <- c(msg, "peaksData and peaksDataFrame have different length")
    if (is.null(msg)) TRUE
    else msg
})

#' @rdname hidden_aliases
setMethod("show", "MsBackendDF", function(object) {
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
#' @rdname MsBackend
setMethod("backendInitialize", signature = "MsBackendDF",
          function(object, data, ...) {
              if (missing(data)) data <- data.frame()
              if (is(data, "DataFrame"))
                  data <- as(data, "data.frame")
              if (!is.data.frame(data))
                  stop("'data' has to be a 'DataFrame' or 'data.frame'")
              if (nrow(data)) {
                  data$dataStorage <- "<memory>"
                  ## Check for *peaks data* columns
                  peaks_cols <- .df_peaks_columns_data_frame(data)
                  ## Get m/z and intensity and put into peaksData
                  p_cols <- intersect(peaks_cols, c("mz", "intensity"))
                  if (length(p_cols)) {
                      object@peaksData <- do.call(
                          mapply, c(list(FUN = cbind, SIMPLIFY = FALSE),
                                    data[p_cols]))
                      for (p in p_cols)
                          data[[p]] <- NULL
                  } else {
                      emat <- matrix(
                          numeric(), ncol = 2, nrow = 0,
                          dimnames = list(character(), c("mz", "intensity")))
                      object@peaksData <- replicate(nrow(data), emat,
                                                    simplify = FALSE)
                  }
                  ## Get any other potential peaks columns and put them
                  ## into peaksDataFrame
                  p_cols <- peaks_cols[!peaks_cols %in% c("mz", "intensity")]
                  if (length(p_cols)) {
                      object@peaksDataFrame <- do.call(
                          mapply, c(list(FUN = cbind.data.frame,
                                         SIMPLIFY = FALSE), data[p_cols]))
                      for (p in p_cols)
                          data[[p]] <- NULL
                  } else
                      object@peaksDataFrame <- list()
              } else {
                  object@peaksData <- list()
                  object@peaksDataFrame <- list()
              }
              object@spectraData <- data
              validObject(object)
              object
          })

## #' @rdname hidden_aliases
## setMethod("backendMerge", "MsBackendDF", function(object, ...) {
##     object <- unname(c(object, ...))
##     not_empty <- lengths(object) > 0
##     if (any(not_empty))
##         res <- .combine_backend_data_frame(object[not_empty])
##     else res <- object[[1L]]
##     validObject(res)
##     res
## })

## Data accessors

#' @rdname hidden_aliases
setMethod("acquisitionNum", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "acquisitionNum")
})

#' @rdname hidden_aliases
setMethod("centroided", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "centroided")
})

#' @rdname hidden_aliases
#'
#' @aliases centroided<-,MsBackendDF-method
setReplaceMethod("centroided", "MsBackendDF", function(object, value) {
    value_len <- length(value)
    value_type <- is.logical(value)
    if (value_type && (value_len == 1L || value_len == length(object)))
        object@spectraData$centroided <- value
    else
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("collisionEnergy", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "collisionEnergy")
})

#' @rdname hidden_aliases
setReplaceMethod("collisionEnergy", "MsBackendDF",
                 function(object, value) {
                     if (!is.numeric(value) || length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$collisionEnergy <- as.numeric(value)
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("dataOrigin", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "dataOrigin")
})

#' @rdname hidden_aliases
setReplaceMethod("dataOrigin", "MsBackendDF", function(object, value) {
    if (!is.character(value) || length(value) != length(object))
        stop("'value' has to be a 'character' of length ", length(object))
    object@spectraData$dataOrigin <- as.character(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("dataStorage", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "dataStorage")
})

#' @rdname hidden_aliases
setReplaceMethod("dataStorage", "MsBackendDF", function(object, value) {
    if (!is.character(value) || length(value) != length(object))
        stop("'value' has to be a 'character' of length ", length(object))
    object@spectraData$dataStorage <- as.character(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("intensity", "MsBackendDF", function(object) {
    if (length(object)) {
        NumericList(.df_pdata_column(object@peaksData, "intensity"),
                    compress = FALSE)
    } else NumericList(compress = FALSE)
})

#' @rdname hidden_aliases
setReplaceMethod("intensity", "MsBackendDF", function(object, value) {
    if (!length(object))
        return(object)
    if (!(is.list(value) || inherits(value, "NumericList")))
        stop("'value' has to be a list or NumericList")
    if (length(value) != length(object))
        stop("length of 'value' has to match the length of 'object'")
    if (!all(lengths(value) == lengths(object)))
        stop("lengths of 'value' has to match the number of peaks ",
             "(i.e. lengths(object))")
    idx <- which(colnames(object@peaksData[[1L]]) == "intensity")
    if (length(idx)) {
        for (i in seq_along(object@peaksData)) {
            object@peaksData[[i]][, idx] <- value[[i]]
        }
    } else {
        for (i in seq_along(object@peaksData)) {
            object@peaksData[[i]] <- cbind(object@peaksData[[i]],
                                           intensity = value[[i]])
        }
    }
    validObject(object)
    object
})

#' @rdname hidden_aliases
#' @importFrom MsCoreUtils vapply1d
setMethod("ionCount", "MsBackendDF", function(object) {
    vapply1d(intensity(object), sum, na.rm = TRUE)
})

#' @rdname hidden_aliases
setMethod("isEmpty", "MsBackendDF", function(x) {
    lengths(intensity(x)) == 0
})

#' @rdname hidden_aliases
setMethod("isolationWindowLowerMz", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "isolationWindowLowerMz")
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowLowerMz", "MsBackendDF",
                 function(object, value) {
                     if (!is.numeric(value) || length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$isolationWindowLowerMz <-
                         as.numeric(value)
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("isolationWindowTargetMz", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "isolationWindowTargetMz")
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowTargetMz", "MsBackendDF",
                 function(object, value) {
                     if (!is.numeric(value) || length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$isolationWindowTargetMz <-
                         as.numeric(value)
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("isolationWindowUpperMz", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "isolationWindowUpperMz")
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowUpperMz", "MsBackendDF",
                 function(object, value) {
                     if (!is.numeric(value) || length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$isolationWindowUpperMz <-
                         as.numeric(value)
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("length", "MsBackendDF", function(x) {
    nrow(x@spectraData)
})

#' @rdname hidden_aliases
setMethod("lengths", "MsBackendDF", function(x, use.names = FALSE) {
    if (length(x))
        lengths(x@peaksData) / ncol(x@peaksData[[1L]])
    else integer()
})

#' @rdname hidden_aliases
setMethod("msLevel", "MsBackendDF", function(object, ...) {
    .get_column(object@spectraData, "msLevel")
})

#' @rdname hidden_aliases
setReplaceMethod("msLevel", "MsBackendDF", function(object, value) {
    if (!is.integer(value) && is.numeric(value))
        value <- as.integer(value)
    if (!is.integer(value) || length(value) != length(object))
        stop("'value' has to be an 'integer' of length ", length(object))
    object@spectraData$msLevel <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("mz", "MsBackendDF", function(object) {
    if (length(object)) {
        NumericList(.df_pdata_column(object@peaksData, "mz"),
                    compress = FALSE)
    } else NumericList(compress = FALSE)
})

#' @rdname hidden_aliases
setReplaceMethod("mz", "MsBackendDF", function(object, value) {
    if (!length(object))
        return(object)
    if (!(is.list(value) || inherits(value, "NumericList")))
        stop("'value' has to be a list or NumericList")
    if (length(value) != length(object))
        stop("length of 'value' has to match the length of 'object'")
    if (!all(lengths(value) == lengths(object)))
        stop("lengths of 'value' has to match the number of peaks ",
             "(i.e. lengths(object))")
    idx <- which(colnames(object@peaksData[[1L]]) == "mz")
    if (length(idx)) {
        for (i in seq_along(object@peaksData)) {
            object@peaksData[[i]][, idx] <- value[[i]]
        }
    } else {
        for (i in seq_along(object@peaksData)) {
            object@peaksData[[i]] <- cbind(object@peaksData[[i]],
                                           mz = value[[i]])
        }
    }
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("peaksData", "MsBackendDF", function(object,
                                               columns = c("mz", "intensity")) {
    if (length(object)) {
        if (!all(columns %in% peaksVariables(object)))
            stop("Some of the requested peaks variables are not ",
                 "available", call. = FALSE)
        ## quick return
        if (length(columns) == 2 &&
            all(columns == colnames(object@peaksData[[1L]])))
            return(object@peaksData)
        pcol <- intersect(columns, c("mz", "intensity"))
        pdcol <- setdiff(columns, c("mz", "intensity"))
        ## request columns only from peaksData
        if (length(pcol) & !length(pdcol))
            return(lapply(object@peaksData,
                          function(z) z[, pcol, drop = FALSE]))
        ## request columns only from peaksDataFrame
        if (length(pdcol) & !length(pcol))
            return(lapply(object@peaksDataFrame,
                          function(z) as.matrix(z[, pdcol, drop = FALSE])))
        ## request columns from both
        mapply(object@peaksData, object@peaksDataFrame, FUN = function(a, b) {
            as.matrix(cbind(a[, pcol, drop = FALSE],
                            b[, pdcol, drop = FALSE])[, columns])
        }, SIMPLIFY = FALSE)
    } else list()
})

#' @rdname hidden_aliases
setReplaceMethod("peaksData", "MsBackendDF", function(object, value) {
    if (length(object)) {
        if (!(is.list(value) || inherits(value, "SimpleList")))
            stop("'value' has to be a list-like object")
        if (length(value) != length(object))
            stop("Length of 'value' has to match length of 'object'")
        if (!is.matrix(value[[1L]]))
            stop("'value' is expected to be a 'list' of 'matrix'")
        cn <- colnames(value[[1L]])
        lcn <- length(cn)
        lapply(value, function(z) {
            cur_cn <- colnames(z)
            if (lcn != length(cur_cn) || !all(cn == cur_cn))
                stop("provided matrices don't have the same column names")
        })
        ## columns mz and intensity go into peaksData.
        if (lcn == 2 && all(cn == c("mz", "intensity")))
            object@peaksData <- value
        else {
            pcn <- intersect(c("mz", "intensity"), cn)
            if (length(pcn))
                object@peaksData <- lapply(
                    value, function(z) z[, pcn, drop = FALSE])
            pcn <- setdiff(cn, c("mz", "intensity"))
            if (length(pcn))
                object@peaksDataFrame <- lapply(
                    value, function(z) as.data.frame(z[, pcn, drop = FALSE]))
        }
        validObject(object)
    }
    object
})

#' @rdname hidden_aliases
setMethod("polarity", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "polarity")
})

#' @rdname hidden_aliases
setReplaceMethod("polarity", "MsBackendDF", function(object, value) {
    value_len <- length(value)
    value_type <- is.numeric(value)
    if (value_type && (value_len == 1L || value_len == length(object)))
        object@spectraData$polarity <- as.integer(value)
    else
        stop("'value' has to be an 'integer' of length 1 or ", length(object))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("precScanNum", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "precScanNum")
})

#' @rdname hidden_aliases
setMethod("precursorCharge", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "precursorCharge")
})

#' @rdname hidden_aliases
setMethod("precursorIntensity", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "precursorIntensity")
})

#' @rdname hidden_aliases
setMethod("precursorMz", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "precursorMz")
})

#' @rdname hidden_aliases
setMethod("rtime", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "rtime")
})

#' @rdname hidden_aliases
setReplaceMethod("rtime", "MsBackendDF", function(object, value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$rtime <- as.numeric(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("scanIndex", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "scanIndex")
})

## #' @rdname hidden_aliases
## setMethod("selectSpectraVariables", "MsBackendDF",
##           function(object, spectraVariables = spectraVariables(object)) {
##               if (!all(spectraVariables %in% spectraVariables(object)))
##                   stop("Spectra variables ",
##                        paste(spectraVariables[!(spectraVariables %in%
##                                                 spectraVariables(object))],
##                              collapse = ", "), " not available")
##               keep <- spectraVariables[spectraVariables %in%
##                                             colnames(object@spectraData)]
##               if (length(keep))
##                   object@spectraData <- object@spectraData[, keep,
##                                                            drop = FALSE]
##               msg <- .valid_spectra_data_required_columns(object@spectraData)
##               if (length(msg))
##                   stop(msg)
##               validObject(object)
##               object
## })

#' @rdname hidden_aliases
setMethod("smoothed", "MsBackendDF", function(object) {
    .get_column(object@spectraData, "smoothed")
})

#' @rdname hidden_aliases
#'
#' @aliases smoothed<-,MsBackendDF-method
setReplaceMethod("smoothed", "MsBackendDF", function(object, value) {
    value_len <- length(value)
    value_type <- is.logical(value)
    if (value_type && (value_len == 1L || value_len == length(object)))
        object@spectraData$smoothed <- value
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
#' @importMethodsFrom S4Vectors lapply
setMethod(
    "spectraData", "MsBackendDF",
    function(object, columns = spectraVariables(object)) {
        as(.df_spectra_data(object, columns), "DataFrame")
    })

#' @rdname hidden_aliases
#'
#' @note: this replaces **all** the data in the backend.
setReplaceMethod("spectraData", "MsBackendDF", function(object, value) {
    if (inherits(value, "DataFrame")) {
        if (length(object) && nrow(value) != length(object))
            stop("'value' has to be a 'DataFrame' with ",
                 length(object), " rows.")
        object <- backendInitialize(new("MsBackendDF"), value)
    } else stop("'value' is expected to be a 'DataFrame'")
    object
})

#' @rdname hidden_aliases
setMethod("spectraNames", "MsBackendDF", function(object) {
    rownames(object@spectraData)
})

#' @rdname hidden_aliases
setReplaceMethod("spectraNames", "MsBackendDF", function(object, value) {
    rownames(object@spectraData) <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("spectraVariables", "MsBackendDF", function(object) {
    unique(c(names(.SPECTRA_DATA_COLUMNS), colnames(object@spectraData),
             peaksVariables(object)))
})

#' @rdname hidden_aliases
setMethod("peaksVariables", "MsBackendDF", function(object) {
    if (length(object)) {
        pv <- colnames(object@peaksData[[1L]])
        if (length(object@peaksDataFrame))
            pv <- c(pv, colnames(object@peaksDataFrame[[1L]]))
        pv
    } else c("mz", "intensity")
})

## #' @rdname hidden_aliases
## setMethod("tic", "MsBackendDF", function(object, initial = TRUE) {
##     if (initial) {
##         if (any(colnames(object@spectraData) == "totIonCurrent"))
##             .get_column(object@spectraData, "totIonCurrent")
##         else rep(NA_real_, times = length(object))
##     } else vapply1d(intensity(object), sum, na.rm = TRUE)
## })

#' @rdname hidden_aliases
setMethod("$", "MsBackendDF", function(x, name) {
    .df_spectra_data(x, name)[, 1]
})

#' @rdname hidden_aliases
setReplaceMethod("$", "MsBackendDF", function(x, name, value) {
    if (is.list(value) || inherits(value, "SimpleList")) {
        ## Check if lengths matches those of peaksData.
        lns <- lengths(x)
        if (!length(value) == length(x) || !all(lns == lengths(value)))
            stop("length of 'value' has to match length of 'x' and the number ",
                 "of values per list-element has to match the number of peaks ",
                 "per spectrum.")
        if (name %in% c("mz", "intensity")) {
            for (i in seq_along(value))
                x@peaksData[[i]][, name] <- value[[i]]
        } else {
            if (length(x@peaksDataFrame)) {
                for (i in seq_along(value))
                    x@peaksDataFrame[[i]][[name]] <- value[[i]]
            } else {
                value <- lapply(value, function(z) {
                    df <- data.frame(z, check.names = FALSE)
                    colnames(df) <- name
                    df
                })
                x@peaksDataFrame <- value
            }
        }
    } else
        x@spectraData[[name]] <- value
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
setMethod("[", "MsBackendDF", function(x, i, j, ..., drop = FALSE) {
    .df_subset(x, i)
})

#' @rdname hidden_aliases
setMethod("split", "MsBackendDF", function(x, f, drop = FALSE, ...) {
    if (!is.factor(f))
        f <- as.factor(f)
    lapply(split(seq_along(x), f, ...), function(i) .df_subset(x, i))
})

#' @rdname hidden_aliases
setMethod("filterAcquisitionNum", "MsBackendDF",
          function(object, n = integer(), dataStorage = character(),
                   dataOrigin = character()) {
    if (!length(n) || !length(object)) return(object)
    if (!is.integer(n)) stop("'n' has to be an integer representing the ",
                             "acquisition number(s) for sub-setting")
    sel_file <- .sel_file(object, dataStorage, dataOrigin)
    sel_acq <- acquisitionNum(object) %in% n & sel_file
    object[sel_acq | !sel_file]
})
