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
                               files = character(),
                               modCount = integer(),
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
        .valid_ms_backend_files(object@files),
        .valid_ms_backend_files_from_file(object@files,
                                          as.vector(fromFile(object))),
        .valid_ms_backend_mod_count(object@files, object@modCount),
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
#' @rdname hidden_aliases
setMethod("backendInitialize", signature = "MsBackendDataFrame",
          definition = function(object, files, spectraData, ...) {
              if (missing(files)) files <- character()
              if (missing(spectraData)) spectraData <- DataFrame()
              object@files <- files
              object@modCount <- integer(length(files))
              if (is.list(spectraData$mz))
                  spectraData$mz <- SimpleList(spectraData$mz)
              if (is.list(spectraData$intensity))
                  spectraData$intensity <- SimpleList(spectraData$intensity)
              object@spectraData <- spectraData
              validObject(object)
              object
          })

#' @rdname hidden_aliases
setMethod("backendMerge", "MsBackendDataFrame", function(object, ...) {
    object <- unname(c(object, ...))
    object <- object[lengths(object) > 0]
    if (length(object) == 1)
        return(object[[1]])
    if (!all(vapply(object, class, character(1)) == class(object[[1]])))
        stop("Can only merge backends of the same type: ", class(object[[1]]))
    mod_counts <- unlist(lapply(object, function(z) z@modCount))
    files <- unlist(lapply(object, function(z) z@files))
    from_file <- unlist(lapply(object, function(z) z@files[fromFile(z)]),
                        use.names = FALSE)
    res <- new(class(object[[1]]))
    res@files <- unique(files)
    suppressWarnings(
        res@spectraData <- do.call(
            .rbind_fill, lapply(object, function(z) z@spectraData))
    )
    res@spectraData$fromFile <- match(from_file, res@files)
    ## modCount: take the largest modCount for all with the same file
    res@modCount <- unname(vapply(
        split(mod_counts, paste0(files))[paste0(res@files)], max, integer(1)))
    validObject(res)
    res
})

## Data accessors

#' @rdname hidden_aliases
setMethod("acquisitionNum", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "acquisitionNum"))
        object@spectraData$acquisitionNum
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("centroided", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "centroided"))
        object@spectraData$centroided
    else rep(NA, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("centroided", "MsBackendDataFrame", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.logical(value) | length(value) != length(object))
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    object@spectraData$centroided <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("collisionEnergy", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "collisionEnergy"))
        object@spectraData$collisionEnergy
    else rep(NA_real_, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("collisionEnergy", "MsBackendDataFrame", function(object, value) {
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$collisionEnergy <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("fromFile", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "fromFile"))
        object@spectraData$fromFile
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("intensity", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "intensity"))
        object@spectraData$intensity
    else {
        lst <- SimpleList(numeric())
        lst[rep(1, times = length(object))]
    }
})

#' @rdname hidden_aliases
setReplaceMethod("intensity", "MsBackendDataFrame", function(object, value) {
    if (!(is.list(value) | inherits(value, "SimpleList")))
        stop("'value' has to be a list")
    if (length(value) != length(object))
        stop("length of 'value' has to match the length of 'object'")
    if (!all(lengths(value) == peaksCount(object)))
        stop("lengths of 'value' has to match the number of peaks ",
             "(i.e. peaksCount(object))")
    if (!is(value, "SimpleList"))
        value <- SimpleList(value)
    object@spectraData$intensity <- value
    object@modCount <- object@modCount + 1L
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
    if (any(colnames(object@spectraData) == "isolationWindowLowerMz"))
        object@spectraData$isolationWindowLowerMz
    else rep(NA_real_, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowLowerMz", "MsBackendDataFrame",
                 function(object, value) {
                     if (!is.numeric(value) | length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$isolationWindowLowerMz <- value
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("isolationWindowTargetMz", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "isolationWindowTargetMz"))
        object@spectraData$isolationWindowTargetMz
    else rep(NA_real_, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowTargetMz", "MsBackendDataFrame",
                 function(object, value) {
                     if (!is.numeric(value) | length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$isolationWindowTargetMz <- value
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("isolationWindowUpperMz", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "isolationWindowUpperMz"))
        object@spectraData$isolationWindowUpperMz
    else rep(NA_real_, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("isolationWindowUpperMz", "MsBackendDataFrame",
                 function(object, value) {
                     if (!is.numeric(value) | length(value) != length(object))
                         stop("'value' has to be a 'numeric' of length ",
                              length(object))
                     object@spectraData$isolationWindowUpperMz <- value
                     validObject(object)
                     object
                 })

#' @rdname hidden_aliases
setMethod("length", "MsBackendDataFrame", function(x) {
    nrow(x@spectraData)
})

#' @rdname hidden_aliases
setMethod("msLevel", "MsBackendDataFrame", function(object, ...) {
    if (any(colnames(object@spectraData) == "msLevel"))
        object@spectraData$msLevel
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("mz", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "mz"))
        object@spectraData$mz
    else {
        lst <- SimpleList(numeric())
        lst[rep(1, times = length(object))]
    }
})

#' @rdname hidden_aliases
setReplaceMethod("mz", "MsBackendDataFrame", function(object, value) {
    if (!(is.list(value) | inherits(value, "SimpleList")))
        stop("'value' has to be a list")
    if (length(value) != length(object))
        stop("length of 'value' has to match the length of 'object'")
    if (!all(lengths(value) == peaksCount(object)))
        stop("lengths of 'value' has to match the number of peaks ",
             "(i.e. peaksCount(object))")
    if (!is(value, "SimpleList"))
        value <- SimpleList(value)
    object@spectraData$mz <- value
    object@modCount <- object@modCount + 1L
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
    if (!(is.list(value) | inherits(value, "SimpleList")))
        stop("'value' has to be a list-like object")
    if (length(value) != length(object))
        stop("Length of 'value' has to match length of 'object'")
    object@modCount <- object@modCount + 1L
    vals <- lapply(value, function(z) z[, 1])
    if (!is(vals, "SimpleList"))
        vals <- SimpleList(vals)
    object@spectraData$mz <- vals
    vals <- lapply(value, function(z) z[, 2])
    if (!is(vals, "SimpleList"))
        vals <- SimpleList(vals)
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
    if (any(colnames(object@spectraData) == "polarity"))
        object@spectraData$polarity
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("polarity", "MsBackendDataFrame", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be an 'integer' of length 1 or ", length(object))
    object@spectraData$polarity <- as.integer(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("precScanNum", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "precScanNum"))
        object@spectraData$precScanNum
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("precursorCharge", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "precursorCharge"))
        object@spectraData$precursorCharge
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("precursorIntensity", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "precursorIntensity"))
        object@spectraData$precursorIntensity
    else rep(NA_real_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("precursorMz", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "precursorMz"))
        object@spectraData$precursorMz
    else rep(NA_real_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("rtime", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "rtime"))
        object@spectraData$rtime
    else rep(NA_real_, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("rtime", "MsBackendDataFrame", function(object, value) {
    if (!is.numeric(value) | length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@spectraData$rtime <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("scanIndex", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "scanIndex"))
        object@spectraData$scanIndex
    else rep(NA_integer_, times = length(object))
})

#' @rdname hidden_aliases
setMethod("selectSpectraVariables", "MsBackendDataFrame",
          function(object, spectraVariables = spectraVariables(object)) {
              if (!all(spectraVariables %in% spectraVariables(object)))
                  stop("Spectra variables ",
                       paste(spectraVariables[!(spectraVariables %in%
                                                spectraVariables(object))],
                             collapse = ", "), " not available")
              to_subset <- spectraVariables[spectraVariables %in%
                                            colnames(object@spectraData)]
              if (length(to_subset))
                  object@spectraData <- object@spectraData[, to_subset,
                                                           drop = FALSE]
              validObject(object)
              object
})

#' @rdname hidden_aliases
setMethod("smoothed", "MsBackendDataFrame", function(object) {
    if (any(colnames(object@spectraData) == "smoothed"))
        object@spectraData$smoothed
    else rep(NA, times = length(object))
})

#' @rdname hidden_aliases
setReplaceMethod("smoothed", "MsBackendDataFrame", function(object, value) {
    if (length(value) == 1)
        value <- rep(value, length(object))
    if (!is.logical(value) | length(value) != length(object))
        stop("'value' has to be a 'logical' of length 1 or ", length(object))
    object@spectraData$smoothed <- value
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
              df_columns <- columns[columns %in% colnames(object@spectraData)]
              res <- object@spectraData[, df_columns, drop = FALSE]
              other_columns <- columns[!(columns %in% colnames(object@spectraData))]
              if (length(other_columns)) {
                  other_res <- lapply(other_columns, .get_spectra_data_column,
                                      x = object)
                  names(other_res) <- other_columns
                  is_mz_int <- names(other_res) %in% c("mz", "intensity")
                  if (!all(is_mz_int))
                      res <- cbind(res, as(other_res[!is_mz_int], "DataFrame"))
                  if (any(names(other_res) == "mz"))
                      res$mz <- if (length(other_res$mz)) other_res$mz else SimpleList()
                  if (any(names(other_res) == "intensity"))
                      res$intensity <- if (length(other_res$intensity)) other_res$intensity else SimpleList()
              }
              res[, columns, drop = FALSE]
          })

#' @rdname hidden_aliases
setReplaceMethod("spectraData", "MsBackendDataFrame", function(object, value) {
    if (inherits(value, "DataFrame")) {
        if (length(object) && nrow(value) != length(object))
            stop("'value' has to be a 'DataFrame' with ", length(object), " rows.")
        if (is.list(value$mz))
            value$mz <- SimpleList(value$mz)
        if (is.list(value$intensity))
            value$intensity <- SimpleList(value$intensity)
    } else {
        if (length(value) == 1)
            value <- rep(value, length(object))
        if (length(value) != length(object))
            stop("length of 'value' has to be ", length(object))
    }
    object@spectraData <- value
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
            object@spectraData$totIonCurrent
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
        value <- SimpleList(value)
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
#' @rdname hidden_aliases
setMethod("[", "MsBackendDataFrame", function(x, i, j, ..., drop = FALSE) {
    if (!missing(j))
        stop("Subsetting by column ('j = ", j, "' is not supported")
    i <- .i_to_index(i, length(x), rownames(x@spectraData))
    x@spectraData <- x@spectraData[i, , drop = FALSE]
    orig_files <- x@files
    files_idx <- unique(fromFile(x))
    x@files <- orig_files[files_idx]
    x@modCount <- x@modCount[files_idx]
    x@spectraData$fromFile <- match(orig_files[fromFile(x)], x@files)
    validObject(x)
    x
})

#' @rdname hidden_aliases
setMethod("filterAcquisitionNum", "MsBackendDataFrame", function(object,
                                                                 n = integer(),
                                                                 file = integer()) {
    if (!length(n) | !length(object)) return(object)
    if (!length(file)) file <- unique(fromFile(object))
    if (!is.integer(n)) stop("'n' has to be an integer representing the ",
                             "acquisition number(s) for sub-setting")
    if (!is.integer(file)) stop("'file' has to be an integer with the index ",
                                "of the file(s) for subsetting")
    sel_file <- fromFile(object) %in% file
    sel_acq <- acquisitionNum(object) %in% n & sel_file
    object[sel_acq | !sel_file]
})

#' @rdname hidden_aliases
setMethod("filterEmptySpectra", "MsBackendDataFrame", function(object) {
    if (!length(object)) return(object)
    object[as.logical(peaksCount(object))]
})

#' @rdname hidden_aliases
setMethod("filterFile", "MsBackendDataFrame",
          function(object, file = integer()) {
              if (length(file)) {
                  file <- .i_to_index(file, length(fileNames(object)),
                                      fileNames(object))
                  file_names <- fileNames(object)
                  object <- object[fromFile(object) %in% file]
                  if (is.unsorted(file))
                      object[order(match(fileNames(object)[fromFile(object)],
                                         file_names[file]))]
                  else object
              } else object
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
