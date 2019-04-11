#' @include hidden_aliases.R
NULL

#' @title Mass spectrometry data backends
#'
#' @aliases class:MsBackend MsBackend-class MsBackendDataFrame-class MsBackendMzR-class
#'
#' @description
#' 
#' Note that the classes described here are not meant to be used
#' directly by the end-users and the material in this man page is
#' aimed at package developers.
#'
#' `MsBackend` is a vitual class that defines what each different
#' backend needs to provide. `MsBackend` objects provide access to
#' mass spectrometry data. Such backends can be generally classified
#' into *in-memory* or *on-disk*, depending where the data, i.e
#' spectra (m/z and intensities) and spectra annotation (MS level,
#' charge, polarity, ...) are stored. Typically, in-memory backends
#' keep all data in memory ensuring fast data access, while on-disk
#' backends store their data on disk and access it on demand. Note
#' that some backends will have a mixed model and store parts of the
#' data in memory and parts on disk, to stike a good balance between
#' fast access and small memory footprint.
#'
#' The *Backend functions and implementation notes for new backend
#' classes* section documents the API that a backend must implement.
#'
#' Currently available backends are:
#' 
#' - `MsBackendDataFrame`: stores all data in memory using a `DataFrame`.
#' 
#' - `MsBackendMzR`: stores the m/z and intensities on-disk in raw
#'    data files (typically `mzML` or `mzXML`) and the spectra
#'    annotation information (header) in memory in a `DataFrame`. This
#'    backend requires the `mzR` package.
#'
#' See below for more details about individual backends.
#'
#' @param columns For `spectraData`: optional `character` with column
#'     names (spectra variables) that should be included in the
#'     returned `DataFrame`.
#'
#' @param files For `backendInitialize`: `character` with the file
#'     names from which the data is/will be imported. Should be set to
#'     `NA_character_` if not applicable.
#'
#' @param initial For `tic`: `logical(1)` whether the initially
#'     reported total ion current should be reported, or whether the
#'     total ion current should be (re)calculated on the actual data
#'     (`initial = FALSE`).
#'
#' @param object Object extending `MsBackend`.
#'
#' @param spectraData For `backendInitialize`: `DataFrame` with spectrum
#'     metadata/data. Format and whether the argument is required depends on
#'     the backend.
#'
#' @param value replacement value for `<-` methods. See individual method
#'     description or expected data type.
#'
#' @param x Object extending `MsBackend`.
#'
#' @param ... Additional arguments.
#'
#' @section Backend functions and implementation notes for new backend classes:
#'
#' New backend classes **must** extend the base `MsBackend` class and
#' **have** to implement the following methods:
#'
#' - `length`: returns the number of spectra in the object.
#'
#' - `acquisitionNum`: returns the acquisition number of each
#'   spectrum. Returns an `integer` of length equal to the number of
#'   spectra (with `NA_integer_` if not available).
#'
#' - `fileNames`: returns a `character` with the file names, or
#'   `NA_character_` if not relevant.
#'
#' - `fromFile`: returns an `integer` vector of length equal to the
#'   number of spectra in `object` indicating the file index from
#'   which spectra originate. If no files are available,
#'   `NA_character_` is returned for all spectra.
#'   
#' - `intensity`: gets the intensity values from the spectra. Returns
#'   a `list` of `numeric` vectors (intensity values for each
#'   spectrum). The length of the `list` is equal to the number of
#'   `spectra` in `object`.
#'
#' - `mz`: gets the mass-to-charge ratios (m/z) from the
#'   spectra. Returns a `list` or length equal to the number of
#'   spectra, each element a `numeric` vector with the m/z values of
#'   one spectrum.
#'
#' - `peaks` returns a `list` of length equal to the number of spectra
#'    in `object`. Each element of the list is a `matrix` with columns
#'    `mz` and `intensity`. For an empty spectrum, a `matrix` with 0
#'    rows and two columns (named `mz` and `intensity`) is returned.
#'
#' - `ionCount`: returns a `numeric` with the sum of intensities for
#'   each spectrum. If the spectrum is empty (see `isEmptt`),
#'   `NA_real_` is returned.
#'
#' - `isCentroided`: a heuristic approach assessing if the spectra in
#'   `object` are in profile or centroided mode. The function takes
#'   the `qtl`th quantile top peaks, then calculates the difference
#'   between adjacent m/z value and returns `TRUE` if the first
#'   quartile is greater than `k`. (See `MSnbase:::.isCentroided` for
#'   the code.)
#'
#' - `isEmpty`: checks whether a spectrum in `object` is empty
#'   (i.e. does not contain any peaks). Returns a `logical` vector of
#'   length equal number of spectra.
#'
#' - `msLevel`: gets the spectra's MS level. Returns an `integer`
#'   vector (of length equal to the number of spectra) with the MS
#'   level for each spectrum (or `NA_integer_` if not available).
#'
#' - `polarity`, `polarity<-`: gets or sets the polarity for each
#'   spectrum.  `polarity` returns an `integer` vector (length equal
#'   to the number of spectra), with `0` and `1` representing negative
#'   and positive polarities, respectively. `polarity<-` expects an
#'   integer vector of length 1 or equal to the number of spectra.
#'
#' - `peaksCount`: gets the number of peaks (m/z-intensity values) per
#'   spectrum.  Returns an `integer` vector (length equal to the
#'   number of spectra). For empty spectra, `NA_integer_` is returned.
#'
#' - `precursorCharge`, `precursorIntensity`, `precursorMz`,
#'   `precScanNum`, `precAcquisitionNum`: get the charge (`integer`),
#'   intensity (`numeric`), m/z (`numeric`), scan index (`integer`)
#'   and acquisition number (`interger`) of the precursor for MS level
#'   > 2 spectra from the object. Returns a vector of length equal to
#'   the number of spectra in `object`. `NA` are reported for MS1
#'   spectra of if no precursor information is available.
#'
#' - `isReadOnly`: returns a `logical(1)` whether the backend is *read
#'   only* or does allow also to write/update data.
#'
#' - `rtime`, `rtime<-`: gets or sets the retention times for each
#'   spectrum.  `rtime` returns a `numeric` vector (length equal to
#'   the number of spectra) with the retention time for each spectrum.
#'   `rtime<-` expects a numeric vector with length equal to the
#'   number of spectra.
#'
#' - `scanIndex`: returns an `integer` vector with the *scan index*
#'   for each spectrum. This represents the relative index of the
#'   spectrum within each file. Note that this can be different to the
#'   `acquisitionNum` of the spectrum which is the index of the
#'   spectrum as reported in the mzML file.
#'
#' - `spectraData`, `spectraData<-`: get or sets general spectrum
#'   metadata (annotation, also called header).  `spectraData` returns
#'   a `DataFrame`, `spectraData<-` expects a `DataFrame`.
#'
#' - `spectraNames`: returns a `character` vector with the names of
#'   the spectra in `object`. If names aren't set, should return
#'   `NULL`.
#'
#' - `spectraVariables`: returns a `character` vector with the
#'   available spectra variables (columns, fields or attributes)
#'   available in `object`.
#'
#' - `tic`: gets the total ion current/count (sum of signal of a
#'   spectrum) for all spectra in `object`. By default, the value
#'   reported in the original raw data file is returned. For an empty
#'   spectrum, `NA_real_` is returned.
#'
#' - `smoothed`,`smoothed<-`: geta or sets whether a spectrum is
#'   *smoothed*. `smoothed` returns a `logical` vector of length equal
#'   to the number of spectra. `smoothed<-` takes a `logical` vector
#'   of length 1 or equal to the number of spectra in `object`.
#'
#' - `centroided`, `centroided<-`: gets or sets the centroiding
#'   information of the spectra. `centroided` returns a `logical`
#'   vector of length equal to the number of spectra with `TRUE` if a
#'   spectrum is centroided, `FALSE` if it is in profile mode and `NA`
#'   if it is undefined. See also `isCentroided` for estimating from
#'   the spectrum data whether the spectrum is centroided.  `value`
#'   for `centroided<-` is either a single `logical` or a `logical` of
#'   length equal to the number of spectra in `object`.
#'
#' - `collisionEnergy`, `collisionEnergy<-`: gets or sets the
#'   collision energy for all spectra in `object`. `collisionEnergy`
#'   returns a `numeric` with length equal to the number of spectra
#'   (`NA_real_` if not present/defined), `collisionEnergy<-` takes a
#'   `numeric` of length equal to the number of spectra in `object`.
#'
#' - `backendInitialize`: initialises the backend. This method is
#'   supposed to be called rights after creating an instance of the
#'   backend class and should prepare the backend (e.g. set the data
#'   for the memory backend or read the spectra header data for the
#'   `MsBackendMzR` backend).
#'
#' @section `MsBackendDataFrame`, in-memory MS data backend:
#'
#' The `MsBackendDataFrame` objects keep all MS data in memory. New
#' objects can be created with the `MsBackendDataFrame()`
#' function. The backend can be subsequently initialized with the
#' `backendInitialize` method, taking a `DataFrame` with the MS data
#' as parameter. Suggested columns of this `DataFrame` are:
#'
#' - `"msLevel"`: `integer` with MS levels of the spectra.
#' - `"rt"`: `numeric` with retention times of the spectra.
#' - `"acquisitionNum"`: `integer` with the acquisition number of the spectrum.
#' - `"scanIndex"`: `integer` with the index of the scan/spectrum within the
#'   *mzML*/*mzXML*/*CDF* file.
#' - `"fromFile"`: `integer` indicating in which file in an experiment the
#'   spectrum was measured.
#' - `"centroided"`: `logical` whether the spectrum is centroided.
#' - `"smoothed"`: `logical` whether the spectrum was smoothed.
#' - `"polarity"`: `integer` with the polarity information of the spectra.
#' - `"precScanNum"`: `integer` specifying the index of the (MS1) spectrum
#'   containing the precursor of a (MS2) spectrum.
#' - `"precursorMz"`: `numeric` with the m/z value of the precursor.
#' - `"precursorIntensity"`: `numeric` with the intensity value of the
#'   precursor.
#' - `"precursorCharge"`: `integer` with the charge of the precursor.
#' - `"collisionEnergy"`: `numeric` with the collision energy.
#' - `"mz"`: `list` of `numeric` vectors representing the m/z values for each
#'   spectrum.
#' - `"intensity"`: `list` of `numeric` vectors representing the intensity
#'   values for each spectrum.
#'
#' Additional columns are allowed too.
#'
#' @section `MsBackendMzR`, on-disk MS data backend:
#'
#' The `MsBackendMzR` keeps only a limited amount of data in memory,
#' while the spectra data (m/z and intensity values) are fetched from
#' the raw files on-demand. This backend uses the `mzR` package for
#' data import and retrieval and hence requires that package to be
#' installed. Also, it can only be used to import and represent data
#' stored in *mzML*, *mzXML* and *CDF* files.
#'
#' New objects can be created with the `MsBackendMzR()` function which
#' can be subsequently filled with data by calling `backendInitialize`
#' passing the file names of the input data files.
#'
#' @name MsBackend
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @md
#'
#' @exportClass MsBackend MsBackendDataFrame MsBackendMzR
NULL

setClass("MsBackend",
         contains = "VIRTUAL",
         slots = c(
             files = "character",  # Can also be NA_character_
             modCount = "integer", # Same length than files.
             readonly = "logical",
             version = "character"),
         prototype = prototype(files = character(),
                               modCount = integer(),
                               readonly = FALSE,
                               version = "0.1"))

#' @importFrom methods .valueClassTest is new validObject
#'
#' @noRd
setValidity("MsBackend", function(object) {
    msg <- c(
        .valid_ms_backend_files(object@files),
        .valid_ms_backend_mod_count(object@files, object@modCount))
    if (is.null(msg)) TRUE
    else msg
})

#' @exportMethod backendInitialize
#'
#' @rdname MsBackend
setMethod("backendInitialize", signature = "MsBackend",
          definition = function(object, files, spectraData, ...) {
              if (missing(files)) files <- character()
              object@files <- files
              object@modCount <- integer(length(files))
              validObject(object)
              object
          })

## Data accessors

#' @exportMethod acquisitionNum
#'
#' @rdname MsBackend
setMethod("acquisitionNum", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod centroided
#'
#' @rdname MsBackend
setMethod("centroided", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod centroided<-
#'
#' @rdname MsBackend
setReplaceMethod("centroided", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod collisionEnergy
#'
#' @rdname MsBackend
setMethod("collisionEnergy", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod collisionEnergy<-
#'
#' @rdname MsBackend
setReplaceMethod("collisionEnergy", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod fromFile
#'
#' @rdname MsBackend
setMethod("fromFile", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
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
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod ionCount
#'
#' @rdname MsBackend
setMethod("ionCount", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod isCentroided
#'
#' @importMethodsFrom ProtGenerics isCentroided
#'
#' @rdname MsBackend
setMethod("isCentroided", "MsBackend", function(object, ...) {
    stop("Not implemented for ", class(object), ".")
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
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod mz
#'
#' @importMethodsFrom ProtGenerics mz
#'
#' @rdname MsBackend
setMethod("mz", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod peaks
#'
#' @importMethodsFrom ProtGenerics peaks
#'
#' @rdname MsBackend
setMethod("peaks", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod peaksCount
#'
#' @rdname MsBackend
setMethod("peaksCount", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod polarity
#'
#' @rdname MsBackend
setMethod("polarity", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod polarity<-
#'
#' @rdname MsBackend
setReplaceMethod("polarity", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precScanNum
#'
#' @rdname MsBackend
setMethod("precScanNum", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precursorCharge
#'
#' @rdname MsBackend
setMethod("precursorCharge", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precursorIntensity
#'
#' @rdname MsBackend
setMethod("precursorIntensity", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precursorMz
#'
#' @rdname MsBackend
setMethod("precursorMz", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod rtime
#'
#' @importMethodsFrom ProtGenerics rtime
#'
#' @rdname MsBackend
setMethod("rtime", "MsBackend", function(object) {


    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod rtime<-
#'
#' @rdname MsBackend
setReplaceMethod("rtime", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod scanIndex
#'
#' @rdname MsBackend
setMethod("scanIndex", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod smoothed
#'
#' @rdname MsBackend
setMethod("smoothed", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod smoothed<-
#'
#' @rdname MsBackend
setReplaceMethod("smoothed", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod spectraData
#'
#' @rdname MsBackend
setMethod("spectraData", "MsBackend", function(object, columns) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod spectraData<-
#'
#' @rdname MsBackend
setReplaceMethod("spectraData", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod spectraNames
#'
#' @rdname MsBackend
setMethod("spectraNames", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod spectraNames<-
#'
#' @rdname MsBackend
setReplaceMethod("spectraNames", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod spectraVariables
#'
#' @rdname MsBackend
setMethod("spectraVariables", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod tic
#'
#' @importMethodsFrom ProtGenerics tic
#'
#' @rdname MsBackend
setMethod("tic", "MsBackend", function(object, initial = TRUE) {
    stop("Not implemented for ", class(object), ".")
})
