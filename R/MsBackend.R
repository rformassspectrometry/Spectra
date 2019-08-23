#' @include hidden_aliases.R
NULL

#' @title Mass spectrometry data backends
#'
#' @aliases class:MsBackend MsBackend-class MsBackendDataFrame-class MsBackendMzR-class [,MsBackend-method
#'
#' @description
#'
#' Note that the classes described here are not meant to be used
#' directly by the end-users and the material in this man page is
#' aimed at package developers.
#'
#' `MsBackend` is a virtual class that defines what each different
#' backend needs to provide. `MsBackend` objects provide access to
#' mass spectrometry data. Such backends can be classified into
#' *in-memory* or *on-disk* backends, depending on where the data, i.e
#' spectra (m/z and intensities) and spectra annotation (MS level,
#' charge, polarity, ...) are stored.
#'
#' Typically, in-memory backends keep all data in memory ensuring fast
#' data access, while on-disk backends store (parts of) their data on
#' disk and retrieve it on demand.
#'
#' The *Backend functions and implementation notes for new backend
#' classes* section documents the API that a backend must implement.
#'
#' Currently available backends are:
#'
#' - `MsBackendDataFrame`: stores all data in memory using a `DataFrame`.
#'
#' - `MsBackendMzR`: stores the m/z and intensities on-disk in raw
#'   data files (typically `mzML` or `mzXML`) and the spectra
#'   annotation information (header) in memory in a `DataFrame`. This
#'   backend requires the `mzR` package.
#'
#' - `MsBackendHdf5Peaks`: stores the m/z and intensities on-disk in custom hdf5
#'   data files and the remaining spectra variables in memory (in a
#'   `DataFrame`). This backend requires the `rhdf5` package.
#'
#' See below for more details about individual backends.
#'
#' @param acquisitionNum for `filterPrecursorScan`: `integer` with the
#'     acquisition number of the spectra to which the object should be
#'     subsetted.
#'
#' @param columns For `spectraData` accessor: optional `character` with column
#'     names (spectra variables) that should be included in the
#'     returned `DataFrame`. By default, all columns are returned.
#'
#' @param dataOrigin For `filterDataOrigin`: `character` to define which
#'     spectra to keep.
#'     For `filterAcquisitionNum`: optionally specify if filtering should occurr
#'     only for spectra of selected `dataOrigin`.
#'
#' @param dataStorage For `filterDataStorage`: `character` to define which
#'     spectra to keep.
#'     For `filterAcquisitionNum`: optionally specify if filtering should occur
#'     only for spectra of selected `dataStorage`.
#'
#' @param drop For `[`: not considered.
#'
#' @param file For `filterFile`: index or name of the file(s) to which the data
#'     should be subsetted.
#'
#' @param initial For `tic`: `logical(1)` whether the initially
#'     reported total ion current should be reported, or whether the
#'     total ion current should be (re)calculated on the actual data
#'     (`initial = FALSE`).
#'
#' @param i For `[`: `integer`, `logical` or `character` to subset the object.
#'
#' @param j For `[`: not supported.
#'
#' @param msLevel `integer` defining the MS level of the spectra to which the
#'     function should be applied. For `filterMsLevel`: the MS level to which
#'     `object` should be subsetted.
#'
#' @param mz For `filterIsolationWindow`: `numeric(1)` with the m/z value to
#'     filter the object. For `filterPrecursorMz`: `numeric(2)` with the lower
#'     and upper m/z boundary.
#'
#' @param n for `filterAcquisitionNum`: `integer` with the acquisition numbers
#'     to filter for.
#'
#' @param name For `$` and `$<-`: the name of the spectra variable to return
#'     or set.
#'
#' @param object Object extending `MsBackend`.
#'
#' @param polarity For `filterPolarity`: `integer` specifying the polarity to
#'     to subset `object`.
#'
#' @param rt for `filterRt`: `numeric(2)` defining the retention time range to
#'     be used to subset/filter `object`.
#'
#' @param spectraData For `backendInitialize`: `DataFrame` with spectrum
#'     metadata/data. This parameter can be empty for `MsBackendMzR` backends
#'     but needs to be provided for `MsBackendDataFrame` backends.
#'
#' @param spectraVariables For `selectSpectraVariables`: `character` with the
#'     names of the spectra variables to which the backend should be subsetted.
#'
#' @param value replacement value for `<-` methods. See individual
#'     method description or expected data type.
#'
#' @param x Object extending `MsBackend`.
#'
#' @param ... Additional arguments.
#'
#' @section Backend functions:
#'
#' New backend classes **must** extend the base `MsBackend` class and
#' **have** to implement the following methods:
#'
#' - `[`: subset the backend. Only subsetting by element (*row*/`i`) is
#'   allowed
#'
#' - `$`, `$<-`: access or set/add a single spectrum variable (column) in the
#'   backend.
#'
#' - `acquisitionNum`: returns the acquisition number of each
#'   spectrum. Returns an `integer` of length equal to the number of
#'   spectra (with `NA_integer_` if not available).
#'
#' - `backendInitialize`: initialises the backend. This method is
#'   supposed to be called rights after creating an instance of the
#'   backend class and should prepare the backend (e.g. set the data
#'   for the memory backend or read the spectra header data for the
#'   `MsBackendMzR` backend). This method has to ensure to set the
#'   spectra variable `dataStorage` correctly.
#'
#' - `backendMerge`: merges (combines) `MsBackend` objects into a single
#'   instance. All objects to be merged have to be of the same type (e.g.
#'   [MsBackendDataFrame()]).
#'
#' - `dataOrigin`: gets a `character` of length equal to the number of spectra
#'   in `object` with the *data origin* of each spectrum. This could e.g. be
#'   the mzML file from which the data was read.
#'
#' - `dataStorage`: gets a `character` of length equal to the number of spectra
#'   in `object` with the data storage of each spectrum. Note that a
#'   `dataStorage` of `NA_character_` is not supported.
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
#' - `filterAcquisitionNum`: filters the object keeping only spectra matching the
#'   provided acquisition numbers (argument `n`). If `dataOrigin` or
#'   `dataStorage` is also provided, `object` is subsetted to the spectra with
#'   an acquisition number equal to `n` **in spectra with matching dataOrigin
#'   or dataStorage values** retaining all other spectra.
#'
#' - `filterDataOrigin`: filters the object retaining spectra matching the
#'   provided `dataOrigin`. Parameter `dataOrigin` has to be of type
#'   `character` and needs to match exactly the data origin value of the
#'   spectra to subset.
#'   `filterDataOrigin` should return the data ordered by the provided
#'   `dataOrigin` parameter, i.e. if `dataOrigin = c("2", "1")` was provided,
#'   the spectra in the resulting object should be ordered accordingly (first
#'   spectra from data origin `"2"` and then from `"1"`).
#'
#' - `filterDataStorage`: filters the object retaining spectra matching the
#'   provided `dataStorage`. Parameter `dataStorage` has to be of type
#'   `character` and needs to match exactly the data storage value of the
#'   spectra to subset.
#'   `filterDataStorage` should return the data ordered by the provided
#'   `dataStorage` parameter, i.e. if `dataStorage = c("2", "1")` was provided,
#'   the spectra in the resulting object should be ordered accordingly (first
#'   spectra from data storage `"2"` and then from `"1"`).
#'
#' - `filterEmptySpectra`: removes empty spectra (i.e. spectra without peaks).
#'
#' - `filterFile`: retains data of files matching the file index or file name
#'    provided with parameter `file`.
#'
#' - `filterIsolationWindow`: retains spectra that contain `mz` in their
#'   isolation window m/z range (i.e. with an `isolationWindowLowerMz` <= `mz`
#'   and `isolationWindowUpperMz` >= `mz`.
#'
#' - `filterMsLevel`: retains spectra of MS level `msLevel`.
#'
#' - `filterPolarity`: retains spectra of polarity `polarity`.
#'
#' - `filterPrecursorMz`: retains spectra with a precursor m/z within the
#'   provided m/z range.
#'
#' - `filterPrecursorScan`: retains parent (e.g. MS1) and children scans (e.g.
#'    MS2) of acquisition number `acquisitionNum`.
#'
#' - `filterRt`: retains spectra of MS level `msLevel` with retention times
#'    within (`>=`) `rt[1]` and (`<=`) `rt[2]`.
#'
#' - `intensity`: gets the intensity values from the spectra. Returns
#'   a [NumericList()] of `numeric` vectors (intensity values for each
#'   spectrum). The length of the `list` is equal to the number of
#'   `spectra` in `object`.
#'
#' - `intensity<-`: replaces the intensity values. `value` has to be a `list`
#'   (or [NumericList()]) of length equal to the number of spectra and the
#'   number of values within each list element identical to the number of
#'   peaks in each spectrum (i.e. the `peaksCount(x)`). Note that just
#'   writeable backends support this method.
#'
#' - `ionCount`: returns a `numeric` with the sum of intensities for
#'   each spectrum. If the spectrum is empty (see `isEmpty`),
#'   `NA_real_` is returned.
#'
#' - `isCentroided`: a heuristic approach assessing if the spectra in
#'   `object` are in profile or centroided mode. The function takes
#'   the `qtl`th quantile top peaks, then calculates the difference
#'   between adjacent m/z value and returns `TRUE` if the first
#'   quartile is greater than `k`. (See `Spectra:::.isCentroided` for
#'   the code.)
#'
#' - `isEmpty`: checks whether a spectrum in `object` is empty
#'   (i.e. does not contain any peaks). Returns a `logical` vector of
#'   length equal number of spectra.
#'
#' - `isolationWindowLowerMz`, `isolationWindowLowerMz<-`: gets or sets the
#'   lower m/z boundary of the isolation window.
#'
#' - `isolationWindowTargetMz`, `isolationWindowTargetMz<-`: gets or sets the
#'   target m/z of the isolation window.
#'
#' - `isolationWindowUpperMz`, `isolationWindowUpperMz<-`: gets or sets the
#'   upper m/z boundary of the isolation window.
#'
#' - `isReadOnly`: returns a `logical(1)` whether the backend is *read
#'   only* or does allow also to write/update data.
#'
#' - `length`: returns the number of spectra in the object.
#'
#' - `msLevel`: gets the spectra's MS level. Returns an `integer`
#'   vector (of length equal to the number of spectra) with the MS
#'   level for each spectrum (or `NA_integer_` if not available).
#'
#' - `mz`: gets the mass-to-charge ratios (m/z) from the
#'   spectra. Returns a [NumericList()] or length equal to the number of
#'   spectra, each element a `numeric` vector with the m/z values of
#'   one spectrum.
#'
#' - `mz<-`: replaces the m/z values. `value` has to be a `list` of length equal
#'   to the number of spectra and the number of values within each list element
#'   identical to the number of peaks in each spectrum (i.e. the
#'   `peaksCount(x)`). Note that just writeable backends support this method.
#'
#' - `peaks` returns a `list` of length equal to the number of spectra
#'   in `object`. Each element of the list is a `matrix` with columns
#'   `mz` and `intensity`. For an empty spectrum, a `matrix` with 0
#'   rows and two columns (named `mz` and `intensity`) is returned.
#'
#' - `peaks<-` replaces the peak data (m/z and intensity values) of the backend.
#'   This method expects a `list` of `matrix` objects with columns `"mz"` and
#'   `"intensity"` that has the same length than the number of spectra in the
#'   backend. Note that just writeable backends support this method.
#'
#' - `peaksCount`: gets the number of peaks (m/z-intensity values) per
#'   spectrum.  Returns an `integer` vector (length equal to the
#'   number of spectra). For empty spectra, `NA_integer_` is returned.
#'
#' - `polarity`, `polarity<-`: gets or sets the polarity for each
#'   spectrum.  `polarity` returns an `integer` vector (length equal
#'   to the number of spectra), with `0` and `1` representing negative
#'   and positive polarities, respectively. `polarity<-` expects an
#'   integer vector of length 1 or equal to the number of spectra.
#'
#' - `precursorCharge`, `precursorIntensity`, `precursorMz`,
#'   `precScanNum`, `precAcquisitionNum`: get the charge (`integer`),
#'   intensity (`numeric`), m/z (`numeric`), scan index (`integer`)
#'   and acquisition number (`interger`) of the precursor for MS level
#'   > 2 spectra from the object. Returns a vector of length equal to
#'   the number of spectra in `object`. `NA` are reported for MS1
#'   spectra of if no precursor information is available.
#'
#' - `rtime`, `rtime<-`: gets or sets the retention times for each
#'   spectrum (in seconds). `rtime` returns a `numeric` vector (length equal to
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
#' - `selectSpectraVariables`: reduces the information within the backend to
#'   the selected spectra variables.
#'
#' - `smoothed`,`smoothed<-`: gets or sets whether a spectrum is
#'   *smoothed*. `smoothed` returns a `logical` vector of length equal
#'   to the number of spectra. `smoothed<-` takes a `logical` vector
#'   of length 1 or equal to the number of spectra in `object`.
#'
#' - `spectraData`, `spectraData<-`: gets or sets general spectrum
#'   metadata (annotation, also called header).  `spectraData` returns
#'   a `DataFrame`, `spectraData<-` expects a `DataFrame` with the same number
#'   of rows as there are spectra in `object`.
#'
#' - `spectraNames`: returns a `character` vector with the names of
#'   the spectra in `object`.
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
#' @section Subsetting and merging backend classes:
#'
#' Backend classes must support (implement) the `[` method to subset the object.
#' This method should only support subsetting by spectra (rows, `i`) and has
#' to return a `MsBackend` class.
#'
#' Backends extending `MsBackend` should also implement the `backendMerge`
#' method to support combining backend instances (only backend classes of the
#' same type should be merged). Merging should follow the following rules:
#'
#' - The whole spectrum data of the various objects should be merged. The
#'   resulting merged object should contain the union of the individual objects'
#'   spectra variables (columns/fields), with eventually missing variables in
#'   one object being filled with `NA`.
#'
#' @section `MsBackendDataFrame`, in-memory MS data backend:
#'
#' The `MsBackendDataFrame` objects keep all MS data in memory.
#' To reduce memory requirement, all spectra variables with a single
#' value (e.g. if all spectra are from MS level 1) are internally represented
#' as an [Rle()] object that are converted into the original class (e.g.
#' `integer`) when the column is accessed.
#'
#' New objects can be created with the `MsBackendDataFrame()`
#' function. The backend can be subsequently initialized with the
#' `backendInitialize` method, taking a `DataFrame` with the MS data
#' as parameter. Suggested columns of this `DataFrame` are:
#'
#' - `"msLevel"`: `integer` with MS levels of the spectra.
#' - `"rt"`: `numeric` with retention times of the spectra.
#' - `"acquisitionNum"`: `integer` with the acquisition number of the spectrum.
#' - `"scanIndex"`: `integer` with the index of the scan/spectrum within the
#'   *mzML*/*mzXML*/*CDF* file.
#' - `"dataOrigin"`: `character` defining the *data origin*.
#' - `"dataStorage"`: `character` indicating grouping of spectra in different
#'   e.g. input files. Note that missing values are not supported.
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
#' - `"mz"`: [NumericList()] of `numeric` vectors representing the m/z values
#'   for each spectrum.
#' - `"intensity"`: [NumericList()] of `numeric` vectors representing the
#'   intensity values for each spectrum.
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
#' The `MsBackendMzR` backend extends the `MsBackendDataFrame` backend using
#' its `DataFrame` to keep spectra variables (except m/z and intensity) in
#' memory.
#'
#' New objects can be created with the `MsBackendMzR()` function which
#' can be subsequently filled with data by calling `backendInitialize`
#' passing the file names of the input data files with argument `files`.
#'
#' @section `MsBackendHdf5Peaks`, on-disk MS data backend:
#'
#' The `MsBackendHdf5Peaks` keeps, similar to the `MsBackendMzR`, peak data
#' (i.e. m/z and intensity values) in custom data files (in HDF5 format) on
#' disk while the remaining spectra variables are kept in memory. This backend
#' supports updating and writing of manipulated peak data to the data files.
#'
#' New objects can be created with the `MsBackendHdf5Peaks()` function which
#' can be subsequently filled with data by calling the object's
#' `backendInitialize` method passing the desired file names of the HDF5 data
#' files along with the spectra variables in form of a `DataFrame` (see
#' `MsBackendDataFrame` for the expected format). An optional parameter
#' `hdf5path` allows to specify the folder where the HDF5 data files should be
#' stored to. If provided, this is added as the path to the submitted file
#' names (parameter `files`).
#'
#' By default `backendInitialize` will store all peak data into a single HDF5
#' file which name has to be provided with the parameter `files`. To store peak
#' data across several HDF5 files `spectraData` has to contain a column
#' `"dataStorage"` that defines the grouping of spectra/peaks into files: peaks
#' for spectra with the same value in `"dataStorage"` are saved into the same
#' HDF5 file. If parameter `files` is omitted, the value in `dataStorage` is
#' used as file name (replacing any file ending with `".h5"`. To specify the
#' file names, `files`' length has to match the number of unique elements in
#' `"dataStorage"`.
#'
#' For details see examples on the [Spectra()] help page.
#'
#' @section Implementation notes:
#'
#' Backends extending `MsBackend` **must** implement all of its methods (listed
#' above). Developers of new `MsBackend`s should follow the
#' `MsBackendDataFrame` implementation.
#'
#' The `MsBackend` defines the following slots:
#'
#' - `@readonly`: `logical(1)` whether the backend supports writing/replacing
#'   of m/z or intensity values.
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
             readonly = "logical",
             version = "character"),
         prototype = prototype(readonly = FALSE,
                               version = "0.1"))

#' @importFrom methods .valueClassTest is new validObject
#'
#' @noRd
setValidity("MsBackend", function(object) {
    msg <- .valid_ms_backend_data_storage(dataStorage(object))
    if (length(dataStorage(object)) != length(object))
        msg <- c(msg, "length of object and 'dataStorage' have to match")
    if (is.null(msg)) TRUE
    else msg
})

#' @exportMethod backendInitialize
#'
#' @rdname MsBackend
setMethod("backendInitialize", signature = "MsBackend",
          definition = function(object, ...) {
              validObject(object)
              object
          })

#' @rdname MsBackend
setMethod("backendMerge", "list", function(object, ...) {
    backendMerge(object[[1]], object[-1])
})

#' @exportMethod backendMerge
#'
#' @rdname MsBackend
setMethod("backendMerge", "MsBackend", function(object, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod acquisitionNum
#'
#' @importMethodsFrom ProtGenerics acquisitionNum
#'
#' @rdname MsBackend
setMethod("acquisitionNum", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod centroided
#'
#' @importMethodsFrom ProtGenerics centroided
#'
#' @rdname MsBackend
setMethod("centroided", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod centroided<-
#'
#' @importMethodsFrom ProtGenerics centroided<-
#'
#' @rdname MsBackend
setReplaceMethod("centroided", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod collisionEnergy
#'
#' @importMethodsFrom ProtGenerics collisionEnergy
#'
#' @rdname MsBackend
setMethod("collisionEnergy", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod collisionEnergy<-
#'
#' @importMethodsFrom ProtGenerics collisionEnergy<-
#'
#' @rdname MsBackend
setReplaceMethod("collisionEnergy", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod dataOrigin
#'
#' @importMethodsFrom ProtGenerics dataOrigin
#'
#' @rdname MsBackend
setMethod("dataOrigin", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod dataOrigin<-
#'
#' @importMethodsFrom ProtGenerics dataOrigin<-
#'
#' @rdname MsBackend
setReplaceMethod("dataOrigin", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod dataStorage
#'
#' @importMethodsFrom ProtGenerics dataStorage
#'
#' @rdname MsBackend
setMethod("dataStorage", "MsBackend", function(object) {
    stop("Method 'dataStorage' is not implemented for ", class(object), ".")
})

#' @exportMethod dataStorage<-
#'
#' @importMethodsFrom ProtGenerics dataStorage<-
#'
#' @rdname MsBackend
setReplaceMethod("dataStorage", "MsBackend", function(object, value) {
    stop("Method 'dataStorage' is not implemented for ", class(object), ".")
})

#' @exportMethod filterAcquisitionNum
#'
#' @importMethodsFrom ProtGenerics filterAcquisitionNum
#'
#' @rdname MsBackend
setMethod("filterAcquisitionNum", "MsBackend", function(object, n, file, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterDataOrigin
#'
#' @importMethodsFrom ProtGenerics filterDataOrigin
#'
#' @rdname MsBackend
setMethod("filterDataOrigin", "MsBackend", function(object, dataOrigin, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterDataStorage
#'
#' @importMethodsFrom ProtGenerics filterDataStorage
#'
#' @rdname MsBackend
setMethod("filterDataStorage", "MsBackend", function(object, dataStorage, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterEmptySpectra
#'
#' @importMethodsFrom ProtGenerics filterEmptySpectra
#'
#' @rdname MsBackend
setMethod("filterEmptySpectra", "MsBackend", function(object, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterIsolationWindow
#'
#' @importMethodsFrom ProtGenerics filterIsolationWindow
#'
#' @rdname MsBackend
setMethod("filterIsolationWindow", "MsBackend", function(object, mz, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterMsLevel
#'
#' @importMethodsFrom ProtGenerics filterMsLevel
#'
#' @rdname MsBackend
setMethod("filterMsLevel", "MsBackend", function(object, msLevel) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterPolarity
#'
#' @importMethodsFrom ProtGenerics filterPolarity
#'
#' @rdname MsBackend
setMethod("filterPolarity", "MsBackend", function(object, polarity) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterPrecursorMz
#'
#' @importMethodsFrom ProtGenerics filterPrecursorMz
#'
#' @rdname MsBackend
setMethod("filterPrecursorMz", "MsBackend", function(object, mz) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterPrecursorScan
#'
#' @importMethodsFrom ProtGenerics filterPrecursorScan
#'
#' @rdname MsBackend
setMethod("filterPrecursorScan", "MsBackend", function(object,
                                                       acquisitionNum, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterRt
#'
#' @importMethodsFrom ProtGenerics filterRt
#'
#' @rdname MsBackend
setMethod("filterRt", "MsBackend", function(object, rt, msLevel, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod intensity
#'
#' @importMethodsFrom ProtGenerics intensity
#'
#' @rdname MsBackend
setMethod("intensity", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod intensity<-
#'
#' @importMethodsFrom ProtGenerics intensity<-
#'
#' @rdname MsBackend
setReplaceMethod("intensity", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod ionCount
#'
#' @importMethodsFrom ProtGenerics ionCount
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

#' @exportMethod isolationWindowLowerMz
#'
#' @importMethodsFrom ProtGenerics isolationWindowLowerMz
#'
#' @rdname MsBackend
setMethod("isolationWindowLowerMz", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod isolationWindowLowerMz<-
#'
#' @importMethodsFrom ProtGenerics isolationWindowLowerMz<-
#'
#' @rdname MsBackend
setReplaceMethod("isolationWindowLowerMz", "MsBackend",
                 function(object, value) {
                     stop("Not implemented for ", class(object), ".")
                 })

#' @exportMethod isolationWindowTargetMz
#'
#' @importMethodsFrom ProtGenerics isolationWindowTargetMz
#'
#' @rdname MsBackend
setMethod("isolationWindowTargetMz", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod isolationWindowTargetMz<-
#'
#' @importMethodsFrom ProtGenerics isolationWindowTargetMz<-
#'
#' @rdname MsBackend
setReplaceMethod("isolationWindowTargetMz", "MsBackend",
                 function(object, value) {
                     stop("Not implemented for ", class(object), ".")
                 })

#' @exportMethod isolationWindowUpperMz
#'
#' @importMethodsFrom ProtGenerics isolationWindowUpperMz
#'
#' @rdname MsBackend
setMethod("isolationWindowUpperMz", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod isolationWindowUpperMz<-
#'
#' @importMethodsFrom ProtGenerics isolationWindowUpperMz<-
#'
#' @rdname MsBackend
setReplaceMethod("isolationWindowUpperMz", "MsBackend",
                 function(object, value) {
                     stop("Not implemented for ", class(object), ".")
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
#' @importMethodsFrom ProtGenerics msLevel
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

#' @exportMethod mz<-
#'
#' @importMethodsFrom ProtGenerics mz<-
#'
#' @rdname MsBackend
setReplaceMethod("mz", "MsBackend", function(object, value) {
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

#' @exportMethod peaks<-
#'
#' @importMethodsFrom ProtGenerics peaks<-
#'
#' @rdname MsBackend
setReplaceMethod("peaks", "MsBackend", function(object, value) {
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
#' @importMethodsFrom ProtGenerics polarity
#'
#' @rdname MsBackend
setMethod("polarity", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod polarity<-
#'
#' @importMethodsFrom ProtGenerics polarity<-
#'
#' @rdname MsBackend
setReplaceMethod("polarity", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precScanNum
#'
#' @importMethodsFrom ProtGenerics precScanNum
#'
#' @rdname MsBackend
setMethod("precScanNum", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precursorCharge
#'
#' @importMethodsFrom ProtGenerics precursorCharge
#'
#' @rdname MsBackend
setMethod("precursorCharge", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precursorIntensity
#'
#' @importMethodsFrom ProtGenerics precursorIntensity
#'
#' @rdname MsBackend
setMethod("precursorIntensity", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precursorMz
#'
#' @importMethodsFrom ProtGenerics precursorMz
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
#' @importMethodsFrom ProtGenerics rtime<-
#'
#' @rdname MsBackend
setReplaceMethod("rtime", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod scanIndex
#'
#' @importMethodsFrom ProtGenerics scanIndex
#'
#' @rdname MsBackend
setMethod("scanIndex", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod selectSpectraVariables
#'
#' @rdname MsBackend
setMethod("selectSpectraVariables", "MsBackend",
          function(object, spectraVariables = spectraVariables(object)) {
              stop("Not implemented for ", class(object), ".")
})

#' @exportMethod smoothed
#'
#' @importMethodsFrom ProtGenerics smoothed
#'
#' @rdname MsBackend
setMethod("smoothed", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod smoothed<-
#'
#' @importMethodsFrom ProtGenerics smoothed<-
#'
#' @rdname MsBackend
setReplaceMethod("smoothed", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod spectraData
#'
#' @importMethodsFrom ProtGenerics spectraData
#'
#' @rdname MsBackend
setMethod("spectraData", "MsBackend",
          function(object, columns = spectraVariables(object)) {
              stop("Not implemented for ", class(object), ".")
          })

#' @exportMethod spectraData<-
#'
#' @importMethodsFrom ProtGenerics spectraData<-
#'
#' @rdname MsBackend
setReplaceMethod("spectraData", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod spectraNames
#'
#' @importMethodsFrom ProtGenerics spectraNames
#'
#' @rdname MsBackend
setMethod("spectraNames", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod spectraNames<-
#'
#' @importMethodsFrom ProtGenerics spectraNames<-
#'
#' @rdname MsBackend
setReplaceMethod("spectraNames", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod spectraVariables
#'
#' @importMethodsFrom ProtGenerics spectraVariables
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

#' @exportMethod [
#'
#' @rdname MsBackend
setMethod("[", "MsBackend", function(x, i, j, ..., drop = FALSE) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod $
#'
#' @rdname MsBackend
setMethod("$", "MsBackend", function(x, name) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod $<-
#'
#' @rdname MsBackend
setReplaceMethod("$", "MsBackend", function(x, name, value) {
    stop("Not implemented for ", class(x), ".")
})
