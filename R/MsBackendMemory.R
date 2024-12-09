#' @include hidden_aliases.R
NULL

#' @title Improved in-memory MS data backend
#'
#' @name MsBackendMemory
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @noRd
#'
#' @exportClass MsBackendMemory
NULL

#' - spectraData: `data.frame` with arbitrary spectra annotations.
#' - peaksData: `list` of `numeric` matrices with m/z and intensity values.
#' - peaksDataFrame: `list` of `data.frame`s with potentially arbitrary peak
#'     annotations.
#'
#' @noRd
setClass("MsBackendMemory",
         contains = "MsBackend",
         slots = c(spectraData = "data.frame",
                   peaksData = "list",
                   peaksDataFrame = "list"),
         prototype = prototype(spectraData = data.frame(),
                               peaksData = list(),
                               peaksDataFrame = list(),
                               readonly = FALSE,
                               version = "0.1"))

setValidity("MsBackendMemory", function(object) {
    msg <- .valid_spectra_data_required_columns(object@spectraData)
    if (length(msg))
        return(msg)
    msg <-
        .valid_column_datatype(object@spectraData, .SPECTRA_DATA_COLUMNS)
    if (length(object@peaksDataFrame) && length(object@peaksData) !=
        length(object@peaksDataFrame))
        msg <- c(msg, "peaksData and peaksDataFrame have different length")
    if (length(object@peaksData) &&
        length(colnames(object@peaksData[[1L]]) %in%
               c("mz", "intensity")) == 1L)
        msg <- c(msg, paste0("If provided, both \"mz\" and \"intensity\" peak ",
                             "variables need to be defined."))
    if (is.null(msg)) TRUE
    else msg
})

#' @rdname hidden_aliases
setMethod("show", "MsBackendMemory", function(object) {
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
setMethod("backendInitialize", signature = "MsBackendMemory",
          function(object, data, peaksVariables = c("mz", "intensity"), ...) {
              if (missing(data)) data <- data.frame()
              if (is(data, "DataFrame"))
                  data <- as(data, "data.frame")
              if (!is.data.frame(data))
                  stop("'data' has to be a 'DataFrame' or 'data.frame'")
              if (nrow(data)) {
                  data$dataStorage <- "<memory>"
                  peaksVariables <- intersect(peaksVariables, colnames(data))
                  ## Check for *peaks data* columns
                  ok <- peaksVariables %in% .df_peaks_columns_data_frame(data)
                  if (any(!ok))
                      warning("'peaksVariables' ",
                              paste0(peaksVariables[!ok], collapse = ", "),
                              " don't have the correct number of elements.")
                  peaksVariables <- peaksVariables[ok]
                  ## Get m/z and intensity and put into peaksData
                  p_cols <- intersect(peaksVariables, c("mz", "intensity"))
                  if (length(p_cols))
                      object@peaksData <- do.call(
                          mapply, c(list(FUN = cbind, SIMPLIFY = FALSE),
                                    data[p_cols]))
                  else
                      object@peaksData <- .df_empty_peaks_data(nrow(data))
                  ## Get any other potential peaks columns and put them
                  ## into peaksDataFrame
                  p_cols <- peaksVariables[!peaksVariables %in%
                                           c("mz", "intensity")]
                  if (length(p_cols)) {
                      object@peaksDataFrame <- do.call(
                          mapply, c(list(FUN = cbind.data.frame,
                                         SIMPLIFY = FALSE), data[p_cols]))
                  } else
                      object@peaksDataFrame <- list()
              } else {
                  object@peaksData <- list()
                  object@peaksDataFrame <- list()
              }
              object@spectraData <- data[, !colnames(data) %in%
                                           peaksVariables]
              validObject(object)
              object
          })

#' @rdname hidden_aliases
setMethod("backendMerge", "MsBackendMemory", function(object, ...) {
    object <- unname(c(object, ...))
    not_empty <- lengths(object) > 0
    if (any(not_empty))
        res <- .df_combine(object[not_empty])
    else res <- object[[1L]]
    validObject(res)
    res
})

#' @rdname hidden_aliases
setMethod("backendRequiredSpectraVariables", "MsBackendMemory",
          function(object, ...) {
              "dataStorage"
          })

## Data accessors

#' @rdname hidden_aliases
setMethod("acquisitionNum", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "acquisitionNum")
})

#' @rdname hidden_aliases
setMethod("centroided", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "centroided")
})

#' @rdname hidden_aliases
#'
#' @aliases centroided<-,MsBackendMemory-method
setReplaceMethod("centroided", "MsBackendMemory", function(object, value) {
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
setMethod("collisionEnergy", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "collisionEnergy")
})

#' @rdname hidden_aliases
setReplaceMethod("collisionEnergy", "MsBackendMemory",
                 function(object, value) {
                     if (!is.numeric(value) || length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$collisionEnergy <- as.numeric(value)
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("dataOrigin", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "dataOrigin")
})

#' @rdname hidden_aliases
setReplaceMethod("dataOrigin", "MsBackendMemory", function(object, value) {
    if (!is.character(value) || length(value) != length(object))
        stop("'value' has to be a 'character' of length ", length(object))
    object@spectraData$dataOrigin <- as.character(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("dataStorage", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "dataStorage")
})

#' @rdname hidden_aliases
setReplaceMethod("dataStorage", "MsBackendMemory", function(object, value) {
    if (!is.character(value) || length(value) != length(object))
        stop("'value' has to be a 'character' of length ", length(object))
    object@spectraData$dataStorage <- as.character(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("extractByIndex", c("MsBackendMemory", "ANY"), function(object, i) {
    slot(object, "spectraData", check = FALSE) <-
        object@spectraData[i, , drop = FALSE]
    if (length(object@peaksData))
        slot(object, "peaksData", check = FALSE) <- object@peaksData[i]
    if (length(object@peaksDataFrame))
        slot(object, "peaksDataFrame", check = FALSE) <-
            object@peaksDataFrame[i]
    object
})

#' @rdname hidden_aliases
setMethod("intensity", "MsBackendMemory", function(object) {
    if (length(object)) {
        NumericList(.df_pdata_column(object@peaksData, "intensity"),
                    compress = FALSE)
    } else NumericList(compress = FALSE)
})

#' @rdname hidden_aliases
setReplaceMethod("intensity", "MsBackendMemory", function(object, value) {
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
setMethod("ionCount", "MsBackendMemory", function(object) {
    vapply1d(intensity(object), sum, na.rm = TRUE)
})

#' @rdname hidden_aliases
setMethod("isEmpty", "MsBackendMemory", function(x) {
    lengths(intensity(x)) == 0
})

#' @rdname hidden_aliases
setMethod("isolationWindowLowerMz", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "isolationWindowLowerMz")
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowLowerMz", "MsBackendMemory",
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
setMethod("isolationWindowTargetMz", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "isolationWindowTargetMz")
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowTargetMz", "MsBackendMemory",
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
setMethod("isolationWindowUpperMz", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "isolationWindowUpperMz")
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowUpperMz", "MsBackendMemory",
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
setMethod("length", "MsBackendMemory", function(x) {
    nrow(x@spectraData)
})

#' @rdname hidden_aliases
setMethod("lengths", "MsBackendMemory", function(x, use.names = FALSE) {
    if (length(x))
        as.integer(lengths(x@peaksData) / ncol(x@peaksData[[1L]]))
    else integer()
})

#' @rdname hidden_aliases
setMethod("msLevel", "MsBackendMemory", function(object, ...) {
    .get_column(object@spectraData, "msLevel")
})

#' @rdname hidden_aliases
setReplaceMethod("msLevel", "MsBackendMemory", function(object, value) {
    if (!is.integer(value) && is.numeric(value))
        value <- as.integer(value)
    if (!is.integer(value) || length(value) != length(object))
        stop("'value' has to be an 'integer' of length ", length(object))
    object@spectraData$msLevel <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("mz", "MsBackendMemory", function(object) {
    if (length(object)) {
        NumericList(.df_pdata_column(object@peaksData, "mz"),
                    compress = FALSE)
    } else NumericList(compress = FALSE)
})

#' @rdname hidden_aliases
setReplaceMethod("mz", "MsBackendMemory", function(object, value) {
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
            if (is.unsorted(value[[i]]))
                stop("m/z values of each spectrum are expected to be ",
                     "increasingly sorted ")
            object@peaksData[[i]][, idx] <- value[[i]]
        }
    } else {
        for (i in seq_along(object@peaksData)) {
            if (is.unsorted(value[[i]]))
                stop("m/z values of each spectrum are expected to be ",
                     "increasingly sorted ")
            object@peaksData[[i]] <- cbind(object@peaksData[[i]],
                                           mz = value[[i]])
        }
    }
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod(
    "peaksData", "MsBackendMemory", function(object,
                                             columns = c("mz", "intensity")) {
        if (length(object)) {
            cns <- colnames(object@peaksData[[1L]])
            if (length(columns) == length(cns) && all(cns == columns))
                return(object@peaksData)
            if (!all(columns %in% peaksVariables(object)))
                stop("Some of the requested peaks variables are not ",
                     "available", call. = FALSE)
            pcol <- intersect(columns, c("mz", "intensity"))
            pdcol <- setdiff(columns, c("mz", "intensity"))
            ## request columns only from peaksData
            if (length(pcol) & !length(pdcol))
                return(lapply(object@peaksData,
                              function(z) z[, pcol, drop = FALSE]))
            ## request columns only from peaksDataFrame
            if (length(pdcol) & !length(pcol))
                return(lapply(object@peaksDataFrame,
                              function(z) z[, pdcol, drop = FALSE]))
            ## request columns from both
            mapply(object@peaksData, object@peaksDataFrame,
                   FUN = function(a, b) {
                       data.frame(a[, pcol, drop = FALSE],
                                  b[, pdcol, drop = FALSE],
                                  check.names = FALSE)[, columns]
                   }, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        } else list()
    })

#' @rdname hidden_aliases
setReplaceMethod("peaksData", "MsBackendMemory", function(object, value) {
    if (length(object)) {
        .check_peaks_data_value(value, length(object))
        cn <- colnames(value[[1L]])
        lcn <- length(cn)
        ## columns mz and intensity go into peaksData.
        if (lcn == 2 && all(cn == c("mz", "intensity")))
            object@peaksData <- lapply(value, base::as.matrix)
        else {
            pcn <- intersect(c("mz", "intensity"), cn)
            if (length(pcn))
                object@peaksData <- lapply(
                    value, function(z) as.matrix(z[, pcn, drop = FALSE],
                                                 rownames.force = FALSE))
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
setMethod("polarity", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "polarity")
})

#' @rdname hidden_aliases
setReplaceMethod("polarity", "MsBackendMemory", function(object, value) {
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
setMethod("precScanNum", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "precScanNum")
})

#' @rdname hidden_aliases
setMethod("precursorCharge", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "precursorCharge")
})

#' @rdname hidden_aliases
setMethod("precursorIntensity", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "precursorIntensity")
})

#' @rdname hidden_aliases
setMethod("precursorMz", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "precursorMz")
})

#' @rdname hidden_aliases
setMethod("rtime", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "rtime")
})

#' @rdname hidden_aliases
setReplaceMethod("rtime", "MsBackendMemory", function(object, value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$rtime <- as.numeric(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("scanIndex", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "scanIndex")
})

#' @rdname hidden_aliases
setMethod("selectSpectraVariables", "MsBackendMemory",
          function(object, spectraVariables = spectraVariables(object)) {
              if (!all(spectraVariables %in% spectraVariables(object)))
                  stop("Spectra variables ",
                       paste(spectraVariables[!(spectraVariables %in%
                                                spectraVariables(object))],
                             collapse = ", "), " not available")
              ## spectraData
              keep <- intersect(spectraVariables, colnames(object@spectraData))
              object@spectraData <- object@spectraData[, keep, drop = FALSE]
              ## peaksData
              if (length(object@peaksData)) {
                  have <- colnames(object@peaksData[[1L]])
                  keep <- intersect(spectraVariables, have)
                  if (!length(keep))
                      object@peaksData <- .df_empty_peaks_data(length(object))
                  else if (length(keep) < length(have))
                      object@peaksData <- lapply(object@peaksData, function(z)
                          z[, keep, drop = FALSE])
              }
              ## peaksDataFrame
              if (length(object@peaksDataFrame)) {
                  have <- colnames(object@peaksDataFrame[[1L]])
                  keep <- intersect(spectraVariables, have)
                  if (!length(keep)) {
                      object@peaksDataFrame <- list()
                  } else {
                      if (length(keep) < length(have))
                          object@peaksDataFrame <- lapply(
                              object@peaksDataFrame, function(z)
                                  z[, keep, drop = FALSE])
                  }
              }
              msg <- .valid_spectra_data_required_columns(
                  object@spectraData, backendRequiredSpectraVariables(object))
              if (length(msg))
                  stop(msg)
              validObject(object)
              object
})

#' @rdname hidden_aliases
setMethod("smoothed", "MsBackendMemory", function(object) {
    .get_column(object@spectraData, "smoothed")
})

#' @rdname hidden_aliases
#'
#' @aliases smoothed<-,MsBackendMemory-method
setReplaceMethod("smoothed", "MsBackendMemory", function(object, value) {
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
    "spectraData", "MsBackendMemory",
    function(object, columns = spectraVariables(object)) {
        as(.df_spectra_data(object, columns), "DataFrame")
    })

#' @rdname hidden_aliases
#'
#' @note: this replaces **all** the data in the backend.
setReplaceMethod("spectraData", "MsBackendMemory", function(object, value) {
    if (inherits(value, "DataFrame")) {
        if (length(object) && nrow(value) != length(object))
            stop("'value' has to be a 'DataFrame' with ",
                 length(object), " rows.")
        object <- backendInitialize(new("MsBackendMemory"), value,
                                    peaksVariables = peaksVariables(object))
    } else stop("'value' is expected to be a 'DataFrame'")
    object
})

#' @rdname hidden_aliases
setMethod("spectraNames", "MsBackendMemory", function(object) {
    rownames(object@spectraData)
})

#' @rdname hidden_aliases
setReplaceMethod("spectraNames", "MsBackendMemory", function(object, value) {
    rownames(object@spectraData) <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("spectraVariables", "MsBackendMemory", function(object) {
    unique(c(names(.SPECTRA_DATA_COLUMNS), colnames(object@spectraData),
             peaksVariables(object)))
})

#' @rdname hidden_aliases
setMethod("peaksVariables", "MsBackendMemory", function(object) {
    if (length(object)) {
        pv <- colnames(object@peaksData[[1L]])
        if (length(object@peaksDataFrame))
            pv <- c(pv, colnames(object@peaksDataFrame[[1L]]))
        pv
    } else c("mz", "intensity")
})

#' @rdname hidden_aliases
setMethod("tic", "MsBackendMemory", function(object, initial = TRUE) {
    if (initial) {
        if (any(colnames(object@spectraData) == "totIonCurrent"))
            .get_column(object@spectraData, "totIonCurrent")
        else rep(NA_real_, times = length(object))
    } else vapply1d(intensity(object), sum, na.rm = TRUE)
})

#' @rdname hidden_aliases
setMethod("$", "MsBackendMemory", function(x, name) {
    .df_spectra_data(x, name)[, 1]
})

#' @rdname hidden_aliases
setReplaceMethod("$", "MsBackendMemory", function(x, name, value) {
    lv <- length(value)
    if (lv && name %in% peaksVariables(x)) {
        if (is.list(value) || inherits(value, "SimpleList")) {
            lns <- lengths(x)
            if (length(value) == length(x) && all(lns == lengths(value))) {
                if (name %in% c("mz", "intensity")) {
                    for (i in seq_along(value))
                        x@peaksData[[i]][, name] <- value[[i]]
                } else {
                    for (i in seq_along(value))
                        x@peaksDataFrame[[i]][[name]] <- value[[i]]
                }
                validObject(x)
                return(x)
            } else
                stop("length of 'value' has to match length of 'x' and the ",
                     "number of values per list-element has to match the ",
                     "number of peaks per spectrum.")

        }
        else stop("'", name, "' is a peaks variable and 'value' is thus ",
                  "expected to be a list of vectors with replacement values")
    }
    ## Support deleting a peaks variable (except m/z and intensity)
    if (is.null(value) && name %in% peaksVariables(x) &&
        !name %in% c("mz", "intensity")) {
        for (i in seq_along(x@peaksDataFrame))
            x@peaksDataFrame[[i]][[name]] <- value
        validObject(x)
        return(x)
    }
    ## Otherwise just add.
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
setMethod("[", "MsBackendMemory", function(x, i, j, ..., drop = FALSE) {
    .df_subset(x, i)
})

#' @importMethodsFrom methods cbind2
#'
#' @rdname hidden_aliases
setMethod("cbind2", signature = c("MsBackendMemory",
                                  "dataframeOrDataFrameOrmatrix"),
          function(x, y = data.frame(), ...) {
              if (is(y, "matrix"))
                  y <- as.data.frame(y)
              if (any(spectraVariables(x) %in% colnames(y)))
                  stop("spectra variables in 'y' are already present in 'x' ",
                       "replacing them is not allowed")

              if (nrow(y) != length(x))
                  stop("Number of row in'y' does not match the number of ",
                       "spectra in 'x'")
              x@spectraData <- cbind(x@spectraData, y)
              validObject(x)
              x
          })

#' @rdname hidden_aliases
setMethod("split", "MsBackendMemory", function(x, f, drop = FALSE, ...) {
    if (!is.factor(f))
        f <- as.factor(f)
    lapply(split(seq_along(x), f, ...), function(i) .df_subset(x, i))
})

#' @rdname hidden_aliases
setMethod("filterAcquisitionNum", "MsBackendMemory",
          function(object, n = integer(), dataStorage = character(),
                   dataOrigin = character()) {
    if (!length(n) || !length(object)) return(object)
    if (!is.integer(n)) stop("'n' has to be an integer representing the ",
                             "acquisition number(s) for sub-setting")
    sel_file <- .sel_file(object, dataStorage, dataOrigin)
    sel_acq <- acquisitionNum(object) %in% n & sel_file
    object[sel_acq | !sel_file]
})
