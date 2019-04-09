#' @include hidden_aliases.R
NULL

#' @title Mass spectrometry data backends
#'
#' @description
#'
#' [MsBackend-class] objects provide access to mass spectrometry data. Such
#' backends can be generally classidied into *in-memory* and *on-disk* backends.
#' In-memory backends, such as [MsBackendMemory()], keep all (spectra) data in
#' memory ensuring fast data access. On-disk backends like the [MsBackendMzR()]
#' keep only part of the data in memory retrieving the remaining data (mostly
#' m/z and intensity values) on-demand from disk.
#'
#' @param files For `backendInitialize`: `character` with the file names from
#'     which the data is/will be imported. Should be set to `NA_character_` if
#'     not applicable.
#'
#' @param object Object extending [MsBackend-class].
#'
#' @param value replacement value for `<-` methods. See individual method
#'     description or expected data type.
#'
#' @param x Object extending [MsBackend-class].
#'
#' @param ... Additional arguments.
#'
#' @section Implementation notes for new backend classes:
#'
#' New backend classes should extend the base `MsBackend` class and **have** to
#' implement the following methods:
#'
#' - `acquisitionNum`: get the acquisition number of each spectrum. Returns an
#'   `integer` of length equal to the number of spectra (with `NA_integer_` if
#'   not available).
#'
#' - `backendInitialize`: initialize the backend. This method is supposed to be
#'   called rights after creating an instance of the backend class and should
#'   prepare the backend (e.g. set the data for the memory backend or read
#'   the spectra header data for the [MsBackendMzR()] backend).
#'
#' - `centroided`, `centroided<-`: get or set the centroiding information of
#'   the spectra. `centroided` returns a `logical` vector of length equal to the
#'   number of spectra with `TRUE` if a spectrum is centroided, `FALSE` if it
#'   is in profile mode and `NA` if it is undefined. See also `isCentroided`
#'   for estimating from the spectrum data whether the spectrum is centroided.
#'   `value` for `centroided<-` is either a single `logical` or a `logical`
#'   of length equal to the number of spectra in `object`.
#'
#' - `collisionEnergy`, `collisionEnergy<-`: get or set the collision energy
#'   for all spectra in `object`. `collisionEnergy` returns a `numeric` with
#'   length equal to the number of spectra (`NA_real_` if not present/defined),
#'   `collisionEnergy<-` takes a `numeric` of length equal to the number of
#'   spectra in `object`.
#'
#' - `fileNames`: get the file names.
#'
#' - `fromFile`: get the file/sample assignment of each spectrum. Returns an
#'   `integer` vector of length equal to the number of spectra in `object`.
#'
#' - `intensity`: get the intensity values from the spectra. Returns a
#'   `list` of `numeric` vectors (intensity values for each spectrum). The
#'   length of the `list` is equal to the number of `spectra` in `object`.
#'
#' - `ionCount`: returns a `numeric` representing the sum of intensities for
#'   each spectrum.
#'
#' - `isCentroided`: a heuristic approach assessing if the spectra in `object`
#'   are in profile or centroided mode. The function takes the `qtl`th quantile
#'   top peaks, then calculates the difference between adjacent M/Z value and
#'   returns `TRUE` if the first quartile is greater than `k`. (See
#'   `MSnbase:::.isCentroided` for the code.)
#'
#' - `isEmpty`: whether a spectrum in `object` is empty (i.e. does not contain
#'   any peaks). Returns a `logical` vector (length equal number of spectra).
#'
#' - `length`: get the number of spectra in the object.
#'
#' - `msLevel`: get the spectra's MS level. Returns an `integer` vector (length
#'   equal to the number of spectra) with the MS level for each spectrum (or
#'   `NA_integer_` if not available).
#'
#' - `mz`: get the mass-to-charge ratios (m/z) from the spectra. Returns a
#'   `list` or length equal to the number of spectra, each element a `numeric`
#'   vector with the m/z values of one spectrum.
#'
#' - `polarity`, `polarity<-`: get or set the polarity for each spectrum.
#'   `polarity` returns an `integer` vector (length equal to the number of
#'   spectra), with `0` and `1` representing negative and positive polarity,
#'   respectively. `polarity<-` expects an integer vector of length 1 or equal
#'   to the number of spectra.
#'
#' - `peaksCount`: get the number of peaks (m/z-intensity values) per spectrum.
#'   Returns an `integer` vector (length equal to the number of spectra).
#'
#' - `precursorCharge`, `precursorIntensity`, `precursorMz`, `precScanNum`: get
#'   the charge (`integer`), intensity (`numeric`), m/z (`numeric`) and scan
#'   index (`integer`) of the precursor for MS level > 2 spectra from the
#'   object. Returns a vector of length equal to the number of spectra in
#'   `object`. `NA` are reported for MS1 spectra of if no precursor information
#'   is available.
#'
#' - `isReadOnly`: returns a `logical(1)` whether the backend is *read only* or
#'   does allow also to write/update data.
#'
#' - `rtime`, `rtime<-`: get or set the retention times for each spectrum.
#'   `rtime` returns a `numeric` vector (length equal to the number of spectra)
#'   with the retention time for each spectrum.
#'   `rtime<-` expects a numeric vector with length equal to the number of
#'   spectra.
#'
#' - `scanIndex`: get an `integer` vector with the *scan index* for each
#'   spectrum. This represents the relative index of the spectrum within each
#'   file (i.e. for each sample). Note that this can be different to the
#'   `acquisitionNum` of the spectrum which is the index of the spectrum as
#'   reported in the mzML file.
#'
#' - `smoothed`,`smoothed<-`: get or set the information whether a spectrum
#'   was *smoothed*. `smoothed` returns a `logical` vector of
#'   length equal to the number of spectra. `smoothed<-` takes a `logical`
#'   vector of length 1 or equal to the number of spectra in `object`.
#'
#' - `spectraData`, `spectraData<-`: get or set general spectrum metadata.
#'   `spectraData` returns a `DataFrame`, `spectraData<-` expects a `DataFrame`.
#'
#' - `spectraNames`: returns a `character` vector with the names of the spectra
#'   in `object`.
#'
#' - `spectraVariables`: returns a `character` vector with the available
#'   spectra variables (columns, fields or attributes) available in `object`.
#'
#' - `tic`: get the total ion current/count (sum of signal of a spectrum) for
#'   all spectra in `object`. By default, the value reported in the original
#'   raw data file is returned.
#'
#' @name MsBackend
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @md
#'
#' @exportClass MsBackend
NULL

setClass("MsBackend",
         contains = "VIRTUAL",
         slots = c(
             files = "character",  # Can also be NA_character_
             modCount = "integer", # Same length than files.
             readonly = "logical",
             version = "character"),
         prototype = prototype(spectraData = DataFrame(),
                               files = character(),
                               modCount = integer(),
                               readonly = FALSE,
                               version = "0.1"))

setValidity("MsBackend", function(object) {
    cat("MsBackend validity\n") # TODO: remove this later - just for debugging
    msg <- c(
        .valid_ms_backend_files(object@files),
        .valid_ms_backend_mod_count(object@files, object@modCount))
    if (is.null(msg)) TRUE
    else msg
})

#' @title Initialize a backend
#'
#' This generic is used to setup and initialize a backend.
#'
#' It has to be reimplemented by each backend and is supposed to initialize
#' the backend or, in case of in-memory backend, import/set the data.
#'
#' @param object An object inheriting from [MsBackend-class]
#'
#' @param files The path to the source (generally .mzML) files.
#'
#' @param spectraData A [S4Vectors::DataFrame-class]
#'
#' @param ... Other arguments passed to the methods.
#'
#' @return A fully operational [Backend-class] derivate.
#'
#' @family Backend generics
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#'
#' @noRd
setMethod("backendInitialize", signature = "MsBackend",
          definition = function(object, files, spectraData, ...) {
              if (missing(files)) files <- character()
              object@files <- files
              object@modCount <- integer(length(files))
              validObject(object)
              object
          })

#' @exportMethod acquisitionNum
#'
#' @rdname MsBackend
setMethod("acquisitionNum", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod centroided
#'
#' @rdname MsBackend
setMethod("centroided", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod centroided<-
#'
#' @rdname MsBackend
setReplaceMethod("centroided", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod collisionEnergy
#'
#' @rdname MsBackend
setMethod("collisionEnergy", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod collisionEnergy<-
#'
#' @rdname MsBackend
setReplaceMethod("collisionEnergy", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod fromFile
#'
#' @rdname MsBackend
setMethod("fromFile", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod fileNames
#'
#' @rdname MsBackend
setMethod("fileNames", "MsBackend", function(object) {
    object@files
})

#' @exportMethod intensity
#'
#' @importMethodsFrom ProtGenerics intensity
#'
#' @rdname MsBackend
setMethod("intensity", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod ionCount
#'
#' @rdname MsBackend
setMethod("ionCount", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod isCentroided
#'
#' @importMethodsFrom ProtGenerics isCentroided
#'
#' @rdname MsBackend
setMethod("isCentroided", "MsBackend", function(object, ...) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod isEmpty
#'
#' @rdname MsBackend
#'
#' @importMethodsFrom S4Vectors isEmpty
setMethod("isEmpty", "MsBackend", function(x) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod isReadOnly
#'
#' @rdname MsBackend
setMethod("isReadOnly", "MsBackend", function(object) {
    object@readonly
})

#' @exportMethod length
#'
#' @rdname MsBackend
setMethod("length", "MsBackend", function(x) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod msLevel
#'
#' @rdname MsBackend
setMethod("msLevel", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod mz
#'
#' @importMethodsFrom ProtGenerics mz
#'
#' @rdname MsBackend
setMethod("mz", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod peaks
#'
#' @importMethodsFrom ProtGenerics peaks
#'
#' @rdname MsBackend
setMethod("peaks", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod peaksCount
#'
#' @rdname MsBackend
setMethod("peaksCount", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod polarity
#'
#' @rdname MsBackend
setMethod("polarity", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod polarity<-
#'
#' @rdname MsBackend
setReplaceMethod("polarity", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod precScanNum
#'
#' @rdname MsBackend
setMethod("precScanNum", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod precursorCharge
#'
#' @rdname MsBackend
setMethod("precursorCharge", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod precursorIntensity
#'
#' @rdname MsBackend
setMethod("precursorIntensity", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod precursorMz
#'
#' @rdname MsBackend
setMethod("precursorMz", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod rtime
#'
#' @importMethodsFrom ProtGenerics rtime
#'
#' @rdname MsBackend
setMethod("rtime", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod rtime<-
#'
#' @rdname MsBackend
setReplaceMethod("rtime", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod scanIndex
#'
#' @rdname MsBackend
setMethod("scanIndex", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod smoothed
#'
#' @rdname MsBackend
setMethod("smoothed", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod smoothed<-
#'
#' @rdname MsBackend
setReplaceMethod("smoothed", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod spectraData
#'
#' @rdname MsBackend
setMethod("spectraData", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod spectraData<-
#'
#' @rdname MsBackend
setReplaceMethod("spectraData", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod spectraNames
#'
#' @rdname MsBackend
setMethod("spectraNames", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod spectraNames<-
#'
#' @rdname MsBackend
setReplaceMethod("spectraNames", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod spectraVariables
#'
#' @rdname MsBackend
setMethod("spectraVariables", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod tic
#'
#' @importMethodsFrom ProtGenerics tic
#'
#' @rdname MsBackend
setMethod("tic", "MsBackend", function(object) {
    stop("Not implemented for ", class(x), ".")
})
