#' @include hidden_aliases.R
NULL

#' @title The Spectra class to manage and access MS data
#'
#' @aliases Spectra-class [,Spectra-method
#' @aliases uniqueMsLevels uniqueMsLevels,Spectra-method
#'
#' @name Spectra
#'
#' @description
#'
#' The `Spectra` class encapsules spectral mass spectrometry data and
#' related metadata.
#'
#' It supports multiple data backends, e.g. in-memory ([MsBackendMemory],
#' [MsBackendDataFrame()]), on-disk as mzML ([MsBackendMzR()]) or HDF5
#' ([MsBackendHdf5Peaks()]).
#'
#' @details
#'
#' The `Spectra` class uses by default a lazy data manipulation strategy,
#' i.e. data manipulations such as performed with `replaceIntensitiesBelow` are
#' not applied immediately to the data, but applied on-the-fly to the spectrum
#' data once it is retrieved. For some backends that allow to write data back
#' to the data storage (such as the [MsBackendMemory()], [MsBackendDataFrame()]
#' and [MsBackendHdf5Peaks()]) it is possible to apply to queue with the
#' `applyProcessing` function. See the *Data manipulation and analysis
#' methods* section below for more details.
#'
#' To apply arbitrary functions to a `Spectra` use the `spectrapply` function
#' (or directly [chunkapply()] for chunk-wise processing). See description of
#' the `spectrapply` function below for details.
#'
#' For details on plotting spectra, see [plotSpectra()].
#'
#' Clarifications regarding scan/acquisition numbers and indices:
#'
#' - A `spectrumId` (or `spectrumID`) is a vendor specific field in
#'   the mzML file that contains some information about the
#'   run/spectrum, e.g.: `controllerType=0 controllerNumber=1
#'   scan=5281 file=2`
#'
#' - `acquisitionNum` is a more a less sanitize spectrum id generated
#'   from the `spectrumId` field by `mzR` (see
#'   [here](https://github.com/sneumann/mzR/blob/master/src/pwiz/data/msdata/MSData.cpp#L552-L580)).
#'
#' - `scanIndex` is the `mzR` generated sequence number of the
#'   spectrum in the raw file (which doesn't have to be the same as
#'   the `acquisitionNum`)
#'
#' See also [this issue](https://github.com/lgatto/MSnbase/issues/525).
#'
#' @section Creation of objects, conversion, changing the backend and export:
#'
#' `Spectra` classes can be created with the `Spectra` constructor function
#' which supports the following formats:
#'
#' - parameter `object` is a `data.frame` or `DataFrame` containing the
#'   spectrum data. The provided `backend` (by default a
#'   [MsBackendMemory-class]) will be initialized with that data.
#'
#' - parameter `object` is a [MsBackend-class] (assumed to be already
#'   initialized).
#'
#' - parameter `object` is missing, in which case it is supposed that the data
#'   is provided by the [MsBackend-class] class passed along with the `backend`
#'   argument.
#'
#' - parameter `object` is of type `character` and is expected to be the file
#'   names(s) from which spectra should be imported. Parameter `source` allows
#'   to define a [MsBackend-class] that is able to import the data from the
#'   provided source files. The default value for `source` is [MsBackendMzR()]
#'   which allows to import spectra data from mzML, mzXML or CDF files.
#'
#' With `...` additional arguments can be passed to the backend's
#' [backendInitialize()] method. Parameter `backend` allows to specify which
#' [MsBackend-class] should be used for data storage.
#'
#' The backend of a `Spectra` object can be changed with the `setBackend`
#' method that takes an instance of the new backend as second parameter
#' `backend`. A call to `setBackend(sps, backend = MsBackendDataFrame())`
#' would for example change the backend of `sps` to the *in-memory*
#' `MsBackendDataFrame`. Changing to a backend is only supported if that
#' backend has a `data` parameter in its `backendInitialize` method and if
#' `supportsSetBackend` returns `TRUE` for that backend. `setBackend` will
#' transfer the full spectra data from the originating backend as a
#' `DataFrame` to the new backend.
#' Most *read-only* backends do not support `setBackend`. It is for example
#' not possible to change the backend to a *read-only* backend (such as
#' the [MsBackendMzR()] backend).
#'
#' The definition of the function is:
#' `setBackend(object, backend, ..., f = dataStorage(object),
#'     BPPARAM = bpparam())` and its parameters are:
#'
#' - parameter `object`: the `Spectra` object.
#'
#' - parameter `backend`: an instance of the new backend, e.g.
#'   `[MsBackendMemory()]`.
#'
#' - parameter `f`: factor allowing to parallelize the change of the backends.
#'   By default the process of copying the spectra data from the original to the
#'   new backend is performed separately (and in parallel) for each file. Users
#'   are advised to use the default setting.
#'
#' - parameter `...`: optional additional arguments passed to the
#'   [backendInitialize()] method of the new `backend`.
#'
#' - parameter `BPPARAM`: setup for the parallel processing. See [bpparam()] for
#'   details.
#'
#' Data from a `Spectra` object can be **exported** to a file with the `export`
#' function. The actual export of the data has to be performed by the `export`
#' method of the [MsBackend] class defined with the mandatory parameter
#' `backend`. Note however that not all backend classes support export of data.
#' From the `MsBackend` classes in the `Spectra` package currently only the
#' `MsBackendMzR` backend supports data export (to mzML/mzXML file(s));
#' see the help page of the [MsBackend-class] for information on its arguments
#' or the examples below or the vignette for examples.
#'
#' The definition of the function is
#' `export(object, backend,  ...)` and its
#' parameters are:
#'
#' - `object`: the `Spectra` object to be exported.
#'
#' - `backend`: instance of a class extending [MsBackend] which supports export
#'   of the data (i.e. which has a defined `export` method).
#'
#' - `...`: additional parameters specific for the `MsBackend` passed with
#'   parameter `backend`.
#'
#'
#' @section Accessing spectra data:
#'
#' - `$`, `$<-`: gets (or sets) a spectra variable for all spectra in `object`.
#'   See examples for details. Note that replacing values of a peaks variable
#'   is not supported with a non-empty processing queue, i.e. if any filtering
#'   or data manipulations on the peaks data was performed. In these cases
#'   [applyProcessing()] needs to be called first to apply all cached data
#'   operations.
#'
#' - `[[`, `[[<-`: access or set/add a single spectrum variable (column) in the
#'   backend.
#'
#' - `acquisitionNum`: returns the acquisition number of each
#'   spectrum. Returns an `integer` of length equal to the number of
#'   spectra (with `NA_integer_` if not available).
#'
#' - `centroided`, `centroided<-`: gets or sets the centroiding
#'   information of the spectra. `centroided` returns a `logical`
#'   vector of length equal to the number of spectra with `TRUE` if a
#'   spectrum is centroided, `FALSE` if it is in profile mode and `NA`
#'   if it is undefined. See also `isCentroided` for estimating from
#'   the spectrum data whether the spectrum is centroided. `value`
#'   for `centroided<-` is either a single `logical` or a `logical` of
#'   length equal to the number of spectra in `object`.
#'
#' - `collisionEnergy`, `collisionEnergy<-`: gets or sets the
#'   collision energy for all spectra in `object`. `collisionEnergy`
#'   returns a `numeric` with length equal to the number of spectra
#'   (`NA_real_` if not present/defined), `collisionEnergy<-` takes a
#'   `numeric` of length equal to the number of spectra in `object`.
#'
#' - `coreSpectraVariables`: returns the *core* spectra variables along with
#'   their expected data type.
#'
#' - `dataOrigin`, `dataOrigin<-`: gets or sets the *data origin* for each
#'   spectrum. `dataOrigin` returns a `character` vector (same length than
#'   `object`) with the origin of the spectra. `dataOrigin<-` expects a
#'   `character` vector (same length than `object`) with the replacement
#'   values for the data origin of each spectrum.
#'
#' - `dataStorage`: returns a `character` vector (same length than `object`)
#'   with the data storage location of each spectrum.
#'
#' - `intensity`: gets the intensity values from the spectra. Returns
#'   a [NumericList()] of `numeric` vectors (intensity values for each
#'   spectrum). The length of the list is equal to the number of
#'   `spectra` in `object`.
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
#' - `containsMz`: checks for each of the spectra whether they contain mass
#'   peaks with an m/z equal to `mz` (given acceptable difference as defined by
#'   parameters `tolerance` and `ppm` - see [common()] for details). Parameter
#'   `which` allows to define whether any (`which = "any"`, the default) or
#'   all (`which = "all"`) of the `mz` have to match. The function returns
#'   `NA` if `mz` is of length 0 or is `NA`.
#'
#' - `containsNeutralLoss`: checks for each spectrum in `object` if it has a
#'   peak with an m/z value equal to its precursor m/z - `neutralLoss` (given
#'   acceptable difference as defined by parameters `tolerance` and `ppm`).
#'   Returns `NA` for MS1 spectra (or spectra without a precursor m/z).
#'
#' - `length`: gets the number of spectra in the object.
#'
#' - `lengths`: gets the number of peaks (m/z-intensity values) per
#'   spectrum. Returns an `integer` vector (length equal to the
#'   number of spectra). For empty spectra, `0` is returned.
#'
#' - `msLevel`: gets the spectra's MS level. Returns an integer vector (names
#'   being spectrum names, length equal to the number of spectra) with the MS
#'   level for each spectrum.
#'
#' - `mz`: gets the mass-to-charge ratios (m/z) from the
#'   spectra. Returns a [NumericList()] or length equal to the number of
#'   spectra, each element a `numeric` vector with the m/z values of
#'   one spectrum.
#'
#' - `peaksData`: gets the *peaks* data for all spectra in `object`. Peaks data
#'   consist of the m/z and intensity values as well as possible additional
#'   annotations (variables) of all peaks of each spectrum. The function
#'   returns a [SimpleList()] of two dimensional arrays (either `matrix` or
#'   `data.frame`), with each array providing the values for the requested
#'   *peak variables* (by default `"mz"` and `"intensity"`). Optional parameter
#'   `columns` is passed to the backend's `peaksData` function to allow
#'   selection of specific (or additional) peaks variables (columns) that
#'   should be extracted (if available). Importantly,
#'   it is **not** guaranteed that each backend supports this parameter (while
#'   each backend must support extraction of `"mz"` and `"intensity"` columns).
#'   Parameter `columns` defaults to `c("mz", "intensity")` but any value
#'   returned from `peaksVariables` is supported.
#'   Note also that it is possible to extract the peak data with
#'   `as(x, "list")` and `as(x, "SimpleList")` as a `list` and `SimpleList`,
#'   respectively. Note however that, in contrast to `peaksData`, `as` does not
#'   support the parameter `columns`.
#'
#' - `peaksVariables`: lists the available variables for mass peaks provided by
#'   the backend. Default peak variables are `"mz"` and `"intensity"` (which
#'   all backends need to support and provide), but some backends might provide
#'   additional variables.
#'   These variables correspond to the column names of the peak data array
#'   returned by `peaksData`.
#'
#' - `polarity`, `polarity<-`: gets or sets the polarity for each
#'   spectrum.  `polarity` returns an `integer` vector (length equal
#'   to the number of spectra), with `0` and `1` representing negative
#'   and positive polarities, respectively. `polarity<-` expects an
#'   `integer` vector of length 1 or equal to the number of spectra.
#'
#' - `precursorCharge`, `precursorIntensity`, `precursorMz`,
#'   `precScanNum`, `precAcquisitionNum`: gets the charge (`integer`),
#'   intensity (`numeric`), m/z (`numeric`), scan index (`integer`)
#'   and acquisition number (`interger`) of the precursor for MS level >
#'   2 spectra from the object. Returns a vector of length equal to
#'   the number of spectra in `object`. `NA` are reported for MS1
#'   spectra of if no precursor information is available.
#'
#' - `rtime`, `rtime<-`: gets or sets the retention times (in seconds)
#'   for each spectrum.  `rtime` returns a `numeric` vector (length
#'   equal to the number of spectra) with the retention time for each
#'   spectrum.  `rtime<-` expects a numeric vector with length equal
#'   to the number of spectra.
#'
#' - `scanIndex`: returns an `integer` vector with the *scan index*
#'   for each spectrum. This represents the relative index of the
#'   spectrum within each file. Note that this can be different to the
#'   `acquisitionNum` of the spectrum which represents the index of the
#'   spectrum during acquisition/measurement (as reported in the mzML file).
#'
#' - `smoothed`,`smoothed<-`: gets or sets whether a spectrum is
#'   *smoothed*. `smoothed` returns a `logical` vector of length equal
#'   to the number of spectra. `smoothed<-` takes a `logical` vector
#'   of length 1 or equal to the number of spectra in `object`.
#'
#' - `spectraData`: gets general spectrum metadata (annotation, also called
#'   header). `spectraData` returns a `DataFrame`. Note that this
#'   method does by default **not** return m/z or intensity values.
#'
#' - `spectraData<-`: **replaces** the full spectra data of the `Spectra`
#'   object with the one provided with `value`. The `spectraData<-` function
#'   expects a `DataFrame` to be passed as value with the same number of rows
#'   as there a spectra in `object`. Note that replacing values of
#'   peaks variables is not supported with a non-empty processing queue, i.e.
#'   if any filtering or data manipulations on the peaks data was performed.
#'   In these cases [applyProcessing()] needs to be called first to apply all
#'   cached data operations and empty the processing queue.
#'
#' - `spectraNames`, `spectraNames<-`: gets or sets the spectra names.
#'
#' - `spectraVariables`: returns a `character` vector with the
#'   available spectra variables (columns, fields or attributes of each
#'   spectrum) available in `object`. Note that `spectraVariables` does not
#'   list the *peak variables* (`"mz"`, `"intensity"` and eventual additional
#'   annotations for each MS peak). Peak variables are returned by
#'   `peaksVariables`.
#'
#' - `tic`: gets the total ion current/count (sum of signal of a
#'   spectrum) for all spectra in `object`. By default, the value
#'   reported in the original raw data file is returned. For an empty
#'   spectrum, `0` is returned.
#'
#' - `uniqueMsLevels`: get the unique MS levels available in `object`. This
#'   function is supposed to be more efficient than `unique(msLevel(object))`.
#'
#' @section Data subsetting, filtering and merging:
#'
#' Subsetting and filtering of `Spectra` objects can be performed with the below
#' listed methods.
#'
#' - `[`: subsets the spectra keeping only selected elements (`i`). The method
#'   **always** returns a `Spectra` object.
#'
#' - `dropNaSpectraVariables`: removes spectra variables (i.e. columns in the
#'   object's `spectraData` that contain only missing values (`NA`). Note that
#'   while columns with only `NA`s are removed, a `spectraData` call after
#'   `dropNaSpectraVariables` might still show columns containing `NA` values
#'   for *core* spectra variables.
#'
#' - `filterAcquisitionNum`: filters the object keeping only spectra matching
#'   the provided acquisition numbers (argument `n`). If `dataOrigin` or
#'   `dataStorage` is also provided, `object` is subsetted to the spectra with
#'   an acquisition number equal to `n` **in spectra with matching dataOrigin
#'   or dataStorage values** retaining all other spectra.
#'   Returns the filtered `Spectra`.
#'
#' - `filterDataOrigin`: filters the object retaining spectra matching the
#'   provided `dataOrigin`. Parameter `dataOrigin` has to be of type
#'   `character` and needs to match exactly the data origin value of the
#'   spectra to subset.
#'   Returns the filtered `Spectra` object (with spectra ordered according to
#'   the provided `dataOrigin` parameter).
#'
#' - `filterDataStorage`: filters the object retaining spectra stored in the
#'   specified `dataStorage`. Parameter `dataStorage` has to be of type
#'   `character` and needs to match exactly the data storage value of the
#'   spectra to subset.
#'   Returns the filtered `Spectra` object (with spectra ordered according to
#'   the provided `dataStorage` parameter).
#'
#' - `filterEmptySpectra`: removes empty spectra (i.e. spectra without peaks).
#'   Returns the filtered `Spectra` object (with spectra in their
#'   original order).
#'
#' - `filterFourierTransformArtefacts`: remove (Orbitrap) fast fourier
#'   artefact peaks from spectra (see examples below). The function iterates
#'   through all intensity ordered peaks in a spectrum and removes all peaks
#'   with an m/z within +/- `halfWindowSize` of the current peak if their
#'   intensity is lower than `threshold` times the current peak's intensity.
#'   Additional parameters `keepIsotopes`, `maxCharge` and `isotopeTolerance`
#'   allow to avoid removing of potential `[13]C` isotope peaks (`maxCharge`
#'   being the maximum charge that should be considered and `isotopeTolerance`
#'   the absolute acceptable tolerance for matching their m/z).
#'   See [filterFourierTransformArtefacts()] for details and background.
#'
#' - `filterIsolationWindow`: retains spectra that contain `mz` in their
#'   isolation window m/z range (i.e. with an `isolationWindowLowerMz` <= `mz`
#'   and `isolationWindowUpperMz` >= `mz`. Returns the filtered `Spectra`
#'   object (with spectra in their original order).
#'
#' - `filterMsLevel`: filters object by MS level keeping only spectra matching
#'   the MS level specified with argument `msLevel`. Returns the filtered
#'   `Spectra` (with spectra in their original order).
#'
#' - `filterMzRange`: filters the object keeping or removing peaks in each
#'   spectrum that are within the provided m/z range. Whether peaks are
#'   retained or removed can be configured with parameter `keep` (default
#'   `keep = TRUE`).
#'
#' - `filterMzValues`: filters the object keeping **all** peaks in each
#'   spectrum that match the provided m/z value(s) (for `keep = TRUE`, the
#'   default) or removing **all** of them (for `keep = FALSE`). The m/z
#'   matching considers also the absolute `tolerance` and m/z-relative
#'   `ppm` values. `tolerance` and `ppm` have to be of length 1.
#'
#' - `filterPolarity`: filters the object keeping only spectra matching the
#'   provided polarity. Returns the filtered `Spectra` (with spectra in their
#'   original order).
#'
#' - `filterPrecursorMzRange` (previously `filterPrecursorMz` which is now
#'   deprecated): retains spectra with a precursor m/z within the
#'   provided m/z range. See examples for details on selecting spectra with
#'   a precursor m/z for a target m/z accepting a small difference in *ppm*.
#'
#' - `filterPrecursorMzValues`: retains spectra with precursor m/z matching any
#'   of the provided m/z values (given `ppm` and `tolerance`). Spectra with
#'   missing precursor m/z value (e.g. MS1 spectra) are dropped.
#'
#' - `filterPrecursorCharge`: retains spectra with the defined precursor
#'   charge(s).
#'
#' - `filterPrecursorScan`: retains parent (e.g. MS1) and children scans (e.g.
#'   MS2) of acquisition number `acquisitionNum`. Returns the filtered
#'   `Spectra` (with spectra in their original order). Parameter `f` allows to
#'   define which spectra belong to the same sample or original data file (
#'   defaults to `f = dataOrigin(object)`).
#'
#' - `filterRt`: retains spectra of MS level `msLevel` with retention
#'   times (in seconds) within (`>=`) `rt[1]` and (`<=`)
#'   `rt[2]`. Returns the filtered `Spectra` (with spectra in their
#'   original order).
#'
#' - `reset`: restores the data to its original state (as much as possible):
#'   removes any processing steps from the lazy processing queue and calls
#'   `reset` on the backend which, depending on the backend, can also undo e.g.
#'   data filtering operations. Note that a `reset` call after `applyProcessing`
#'   will not have any effect. See examples below for more information.
#'
#' - `selectSpectraVariables`: reduces the information within the object to
#'   the selected spectra variables: all data for variables not specified will
#'   be dropped. For mandatory columns (i.e., those listed by
#'   [coreSpectraVariables()], such as *msLevel*, *rtime* ...) only
#'   the values will be dropped but not the variable itself. Additional (or
#'   user defined) spectra variables will be completely removed.
#'   Returns the filtered `Spectra`.
#'
#' - `split`: splits the `Spectra` object based on parameter `f` into a `list`
#'   of `Spectra` objects.
#'
#' - `joinSpectraData`: Individual spectra variables can be directly
#'    added with the `$<-` or `[[<-` syntax. The `joinSpectraData()`
#'    function allows to merge a `DataFrame` to the existing spectra
#'    data. This function diverges from the [merge()] method in two
#'    main ways:
#'     - The `by.x` and `by.y` column names must be of length 1.
#'     - If variable names are shared in `x` and `y`, the spectra
#'       variables of `x` are not modified. It's only the `y`
#'       variables that are appended the suffix defined in
#'       `suffix.y`. This is to avoid modifying any core spectra
#'       variables that would lead to an invalid object.
#'
#'    - Duplicated Spectra keys (i.e. `x[[by.x]]`) are not
#'      allowed. Duplicated keys in the `DataFrame` (i.e `y[[by.y]]`)
#'      throw a warning and only the last occurrence is kept. These
#'      should be explored and ideally be removed using for
#'      `QFeatures::reduceDataFrame()`, `PMS::reducePSMs()` or similar
#'      functions.
#'
#' Several `Spectra` objects can be concatenated into a single object with the
#' `c` or the `concatenateSpectra` function. Concatenation will fail if the
#' processing queue of any of the `Spectra` objects is not empty or if
#' different backends are used in the `Spectra` objects. The spectra variables
#' of the resulting `Spectra` object is the union of the spectra variables of
#' the individual `Spectra` objects.
#'
#'
#' @section Data manipulation and analysis methods:
#'
#' Many data manipulation operations, such as those listed in this section, are
#' not applied immediately to the spectra, but added to a
#' *lazy processing/manipulation queue*. Operations stored in this queue are
#' applied on-the-fly to spectra data each time it is accessed. This lazy
#' execution guarantees the same functionality for `Spectra` objects with
#' any backend, i.e. backends supporting to save changes to spectrum data
#' ([MsBackendMemory], [MsBackendDataFrame()] or [MsBackendHdf5Peaks()]) as
#' well as read-only backends (such as the [MsBackendMzR()]).
#' Note that for the former it is possible to apply the processing queue and
#' write the modified peak data back to the data storage with the
#' `applyProcessing` function.
#'
#' - `addProcessing`: adds an arbitrary function that should be applied to the
#'   peaks matrix of every spectrum in `object`. The function (can be passed
#'   with parameter `FUN`) is expected to take a peaks matrix as input and to
#'   return a peaks matrix. A peaks matrix is a numeric matrix with two columns,
#'   the first containing the m/z values of the peaks and the second the
#'   corresponding intensities. The function has to have `...` in its
#'   definition. Additional arguments can be passed with `...`. With parameter
#'   `spectraVariables` it is possible to define additional spectra variables
#'   from `object` that should be passed to the function `FUN`. These will be
#'   passed by their name (e.g. specifying `spectraVariables = "precursorMz"`
#'   will pass the spectra's precursor m/z as a parameter named `precursorMz`
#'   to the function. The only exception is the spectra's MS level, these will
#'   be passed to the function as a parameter called `spectrumMsLevel` (i.e.
#'   with `spectraVariables = "msLevel"` the MS levels of each spectrum will be
#'   submitted to the function as a parameter called `spectrumMsLevel`).
#'   Examples are provided in the package vignette.
#'
#' - `applyProcessing`: for `Spectra` objects that use a **writeable** backend
#'   only: apply all steps from the lazy processing queue to the peak data and
#'   write it back to the data storage. Parameter `f` allows to specify how
#'   `object` should be split for parallel processing. This should either be
#'   equal to the `dataStorage`, or `f = rep(1, length(object))` to disable
#'   parallel processing alltogether. Other partitionings might result in
#'   errors (especially if a `MsBackendHdf5Peaks` backend is used).
#'
#' - `bin`: aggregates individual spectra into discrete (m/z) bins. Binning is
#'   performed only on spectra of the specified MS level(s) (parameter
#'   `msLevel`, by default all MS levels of `x`). The bins can be defined with
#'   parameter `breaks` which by default are equally sized bins, with size
#'   being defined by parameter `binSize`, from the minimal to the maximal m/z
#'   of all spectra (of MS level `msLevel`) within `x`. The same bins are used
#'   for all spectra in `x`. All intensity values for peaks falling into the
#'   same bin are aggregated using the function provided with parameter `FUN`
#'   (defaults to `FUN = sum`, i.e. all intensities are summed up). Note that
#'   the binning operation is applied to the peak data on-the-fly upon data
#'   access and it is possible to *revert* the operation with the `reset`
#'   function (see description of `reset` above).
#'
#' - `combineSpectra`: combine sets of spectra into a single spectrum per set.
#'   For each spectrum group (set), spectra variables from the first spectrum
#'   are used and the peak matrices are combined using the function specified
#'   with `FUN`, which defaults to [combinePeaks()]. Please refer to the
#'   [combinePeaks()] help page for details and options of the actual
#'   combination of peaks across the sets of spectra and to the package vignette
#'   for examples and alternative ways to aggregate spectra.
#'   The sets of spectra can be specified with parameter `f`.
#'   In addition it is possible to define, with parameter `p` if and how to
#'   split the input data for parallel processing.
#'   This defaults to `p = x$dataStorage` and hence a per-file parallel
#'   processing is applied for `Spectra` with file-based backends (such as the
#'   [MsBackendMzR()]).
#'   Prior combination of the spectra all processings queued in the lazy
#'   evaluation queue are applied. Be aware that calling `combineSpectra` on a
#'   `Spectra` object with certain backends that allow modifications might
#'   **overwrite** the original data. This does not happen with a
#'   `MsBackendMemory` or `MsBackendDataFrame` backend, but with a
#'   `MsBackendHdf5Peaks` backend the m/z and intensity values in the original
#'   hdf5 file(s) will be overwritten.
#'   The function returns a `Spectra` of length equal to the unique levels
#'   of `f`.
#'
#' - `compareSpectra`: compare each spectrum in `x` with each spectrum in `y`
#'   using the function provided with `FUN` (defaults to [ndotproduct()]). If
#'   `y` is missing, each spectrum in `x` is compared with each other spectrum
#'   in `x`.
#'   The matching/mapping of peaks between the compared spectra is done with the
#'   `MAPFUN` function. The default [joinPeaks()] matches peaks of both spectra
#'   and allows to keep all peaks from the first spectrum (`type = "left"`),
#'   from the second (`type = "right"`), from both (`type = "outer"`) and to
#'   keep only matching peaks (`type = "inner"`); see [joinPeaks()] for more
#'   information and examples). The `MAPFUN` function should have parameters
#'   `x`, `y`, `xPrecursorMz` and `yPrecursorMz` as these values are passed to
#'   the function. In addition to `joinPeaks()` also [joinPeaksGnps()] is
#'   supported for GNPS-like similarity score calculations. Note that
#'   `joinPeaksGnps` should only be used in combination with
#'   `FUN = MsCoreUtils::gnps` (see [joinPeaksGnps()] for more information and
#'   details).
#'   `FUN` is supposed to be a function to compare intensities of (matched)
#'   peaks of the two spectra that are compared. The function needs to take two
#'   matrices with columns `"mz"` and `"intensity"` as input and is supposed
#'   to return a single numeric as result. In addition to the two peak matrices
#'   the spectra's precursor m/z values are passed to the function as parameters
#'   `xPrecursorMz` (precursor m/z of the `x` peak matrix) and `yPrecursorMz`
#'   (precursor m/z of the `y` peak matrix). Additional parameters to functions
#'   `FUN` and `MAPFUN` can be passed with `...`.
#'   The function returns a `matrix` with the results of `FUN` for each
#'   comparison, number of rows equal to `length(x)` and number of columns
#'   equal `length(y)` (i.e. element in row 2 and column 3 is the result from
#'   the comparison of `x[2]` with `y[3]`). If `SIMPLIFY = TRUE` the `matrix`
#'   is *simplified* to a `numeric` if length of `x` or `y` is one.
#'
#' - `deisotopeSpectra`: *deisotope* each spectrum keeping only the monoisotopic
#'   peak for groups of isotopologues. Isotopologues are estimated using the
#'   [isotopologues()] function from the *MetaboCoreUtils* package. Note that
#'   the default parameters for isotope prediction/detection have been
#'   determined using data from the Human Metabolome Database (HMDB) and
#'   isotopes for elements other than CHNOPS might not be detected. See
#'   parameter `substDefinition` in the documentation of [isotopologues()] for
#'   more information. The approach and code to define the parameters for
#'   isotope prediction is described
#'   [here](https://github.com/EuracBiomedicalResearch/isotopologues).
#'
#' - `estimatePrecursorIntensity`: define the precursor intensities for MS2
#'   spectra using the intensity of the matching MS1 peak from the
#'   closest MS1 spectrum (i.e. the last MS1 spectrum measured before the
#'   respective MS2 spectrum). With `method = "interpolation"` it is also
#'   possible to calculate the precursor intensity based on an interpolation of
#'   intensity values (and retention times) of the matching MS1 peaks from the
#'   previous and next MS1 spectrum. See below for an example.
#'
#' - `filterIntensity`: filters each spectrum keeping only peaks with
#'   intensities that are within the provided range or match the criteria of
#'   the provided function. For the former, parameter `intensity` has to be a
#'   `numeric` defining the intensity range, for the latter a `function` that
#'   takes the intensity values of the spectrum and returns a `logical` whether
#'   the peak should be retained or not (see examples below for details) -
#'   additional parameters to the function can be passed with `...`. To
#'   remove only peaks with intensities below a certain threshold, say 100, use
#'   `intensity = c(100, Inf)`. Note: also a single value can be passed with
#'   the `intensity` parameter in which case an upper limit of `Inf` is used.
#'   Note that this function removes also peaks with missing intensities
#'   (i.e. an intensity of `NA`). Parameter `msLevel.` allows to restrict the
#'   filtering to spectra of the specified MS level(s).
#'
#' - `neutralLoss`: calculate neutral loss spectra for fragment spectra. See
#'   [neutralLoss()] for detailed documentation.
#'
#' - `processingLog`: returns a `character` vector with the processing log
#'   messages.
#'
#' - `reduceSpectra`: for groups of peaks within highly similar m/z values
#'   within each spectrum (given `ppm` and `tolerance`), this function keeps
#'   only the peak with the highest intensity removing all other peaks hence
#'   *reducing* each spectrum to the highest intensity peaks per *peak group*.
#'   Peak groups are defined using the [group()] function from the
#'   *MsCoreUtils* package.
#'
#' - `spectrapply`: apply a given function to each individual spectrum or sets
#'   of a `Spectra` object. By default, the `Spectra` is split into individual
#'   spectra (i.e. `Spectra` of length 1) and the function `FUN` is applied to
#'   each of them. An alternative splitting can be defined with parameter `f`.
#'   Parameters for `FUN` can be passed using `...`.
#'   The returned result and its order depend on the function `FUN` and how
#'   `object` is split (hence on `f`, if provided). Parallel processing is
#'   supported and can be configured with parameter `BPPARAM`, is however only
#'   suggested for computational intense `FUN`.
#'   As an alternative to the (eventual parallel) processing of the full
#'   `Spectra`, `spectrapply` supports also a chunk-wise processing. For this,
#'   parameter `chunkSize` needs to be specified. `object` is then split into
#'   chunks of size `chunkSize` which are then (stepwise) processed by `FUN`.
#'   This guarantees a lower memory demand (especially for on-disk backends)
#'   since only the data for one chunk needs to be loaded into memory in each
#'   iteration. Note that by specifying `chunkSize`, parameters `f` and
#'   `BPPARAM` will be ignored.
#'   See also [chunkapply()] or examples below for details on chunk-wise
#'   processing.
#'
#' - `smooth`: smooth individual spectra using a moving window-based approach
#'   (window size = `2 * halfWindowSize`). Currently, the
#'   Moving-Average- (`method = "MovingAverage"`),
#'   Weighted-Moving-Average- (`method = "WeightedMovingAverage")`,
#'   weights depending on the distance of the center and calculated
#'   `1/2^(-halfWindowSize:halfWindowSize)`) and
#'   Savitzky-Golay-Smoothing (`method = "SavitzkyGolay"`) are supported.
#'   For details how to choose the correct `halfWindowSize` please see
#'   [`MsCoreUtils::smooth()`].
#'
#' - `pickPeaks`: picks peaks on individual spectra using a moving window-based
#'   approach (window size = `2 * halfWindowSize`). For noisy spectra there
#'   are currently two different noise estimators available,
#'   the *M*edian *A*bsolute *D*eviation (`method = "MAD"`) and
#'   Friedman's Super Smoother (`method = "SuperSmoother"`),
#'   as implemented in the [`MsCoreUtils::noise()`].
#'   The method supports also to optionally *refine* the m/z value of
#'   the identified centroids by considering data points that belong (most
#'   likely) to the same mass peak. Therefore the m/z value is calculated as an
#'   intensity weighted average of the m/z values within the peak region.
#'   The peak region is defined as the m/z values (and their respective
#'   intensities) of the `2 * k` closest signals to the centroid or the closest
#'   valleys (`descending = TRUE`) in the `2 * k` region. For the latter the `k`
#'   has to be chosen general larger. See [`MsCoreUtils::refineCentroids()`] for
#'   details.
#'   If the ratio of the signal to the highest intensity of the peak is below
#'   `threshold` it will be ignored for the weighted average.
#'
#' - `replaceIntensitiesBelow`: replaces intensities below a specified
#'   threshold with the provided `value`. Parameter `threshold` can be either
#'   a single numeric value or a function which is applied to all non-`NA`
#'   intensities of each spectrum to determine a threshold value for each
#'   spectrum. The default is `threshold = min` which replaces all values
#'   which are <= the minimum intensity in a spectrum with `value` (the
#'   default for `value` is `0`). Note that the function specified with
#'   `threshold` is expected to have a parameter `na.rm` since `na.rm = TRUE`
#'   will be passed to the function. If the spectrum is in profile mode,
#'   ranges of successive non-0 peaks <= `threshold` are set to 0.
#'   Parameter `msLevel.` allows to apply this to only spectra of certain MS
#'   level(s).
#'
#' @section Parallel processing:
#'
#' Some `Spectra` functions have build-in parallel processing that can be
#' configured by passing the parallel processing setup with the `BPPARAM`
#' function argument (which defaults to `BPPARAM = bpparam()`, thus uses
#' the default set up). Most functions have an additional parameter `f` that
#' allows to define how `Spectra` will be split to perform parallel processing.
#' This parameter `f` defaults to `f = dataStorage(object)` and hence
#' parallel processing is performed *by file* (if a file-based, on-disk
#' backend such as `MsBackendMzR` is used). Some `MsBackend` classes might
#' however not support parallel processing. The `backendBpparam` function
#' allows to evaluate wheter a `Spectra` (respectively its `MsBackend`)
#' supports a certain parallel processing setup. Calling
#' `backendBpparam(sps, BPPARAM = MulticoreParam(3))` on a `Spectra` object
#' `sps` would return `SerialParam()` in case the backend of the `Spectra`
#' object does not support parallel processing. All functions listed below
#' use this same function to eventually disable parallel processing to
#' avoid failure of a function call.
#'
#' The functions with build-in parallel processing capabilities are:
#'
#' - `applyProcessing`.
#' - `combineSpectra`.
#' - `containsMz` (does not provide a parameter `f`, but performs parallel
#'   processing separate for `dataStorage`).
#' - `containsNeutralLoss` (same as `containsMz`).
#' - `estimatePrecursorIntensity`.
#' - `setBackend`.
#' - `Spectra` (that passes the `BPPARAM` to the `backendInitialize` of the
#'   used `MsBackend`).
#' - `spectrapply`.
#'
#'
#' @return See individual method description for the return value.
#'
#' @param acquisitionNum for `filterPrecursorScan`: `integer` with the
#'     acquisition number of the spectra to which the object should be
#'     subsetted.
#'
#' @param backend For `Spectra`: [MsBackend-class] to be used as backend. See
#'     section on creation of `Spectra` objects for details. For `setBackend`:
#'     instance of [MsBackend-class] that supports `setBackend` (i.e. for
#'     which `supportsSetBackend` returns `TRUE`). Such backends have a
#'     parameter `data` in their `backendInitialize` function that support
#'     passing the full spectra data to the initialize method. See section on
#'     creation of `Spectra` objects for details.
#'     For `export`: [MsBackend-class] to be used to export the data.
#'
#' @param binSize For `bin`: `numeric(1)` defining the size for the m/z bins.
#'     Defaults to `binSize = 1`.
#'
#' @param BPPARAM Parallel setup configuration. See [bpparam()] for more
#'     information. This is passed directly to the [backendInitialize()] method
#'     of the [MsBackend-class].
#'
#' @param breaks For `bin`: `numeric` defining the m/z breakpoints between bins.
#'
#' @param by.x A `character(1)` specifying the spectra variable used
#'     for merging. Default is `"spectrumId"`.
#'
#' @param by.y A `character(1)` specifying the column used for
#'     merging. Set to `by.x` if missing.
#'
#' @param charge For `deisotopeSpectra`: expected charge of the ionized
#'     compounds. See [isotopologues()] for details.
#'
#' @param chunkSize For `spectrapply`: size of the chunks into which `Spectra`
#'     should be split. This parameter overrides parameters `f` and `BPPARAM`.
#'
#' @param columns For `spectraData` accessor: optional `character` with column
#'     names (spectra variables) that should be included in the
#'     returned `DataFrame`. By default, all columns are returned.
#'     For `peaksData` accessor: optional `character` with requested columns in
#'     the individual `matrix` of the returned `list`. Defaults to
#'     `c("mz", "value")` but any values returned by `peaksVariables(object)`
#'     with `object` being the `Spectra` object are supported.
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
#' @param descending For `pickPeaks`: `logical`, if `TRUE` just values between
#'     the nearest valleys around the peak centroids are used.
#
#' @param drop For `[`, `split`: not considered.
#'
#' @param f For `split`: factor defining how to split `x`. See [base::split()]
#'     for details. For `setBackend`: factor defining how to split the data for
#'     parallelized copying of the spectra data to the new backend. For some
#'     backends changing this parameter can lead to errors.
#'     For `combineSpectra`: `factor` defining the grouping of the spectra that
#'     should be combined. For `spectrapply`: `factor` how `object` should be
#'     splitted. For `estimatePrecursorIntensity` and `filterPrecursorScan`:
#'     defining which spectra belong to the same original data file (sample).
#'     Defaults to `f = dataOrigin(x)`.
#'
#' @param FUN For `addProcessing`: function to be applied to the peak matrix
#'     of each spectrum in `object`. For `compareSpectra`: function to compare
#'     intensities of peaks between two spectra with each other.
#'     For `combineSpectra`: function to combine the (peak matrices) of the
#'     spectra. See section *Data manipulations* and examples below for more
#'     details.
#'     For `bin`: function to aggregate intensity values of peaks falling into
#'     the same bin. Defaults to `FUN = sum` thus summing up intensities.
#'     For `spectrapply` and `chunkapply`: function to be applied to `Spectra`.
#'
#' @param halfWindowSize
#' - For `pickPeaks`: `integer(1)`, used in the
#'   identification of the mass peaks: a local maximum has to be the maximum
#'   in the window from `(i - halfWindowSize):(i + halfWindowSize)`.
#' - For `smooth`: `integer(1)`, used in the smoothing algorithm, the window
#'   reaches from `(i - halfWindowSize):(i + halfWindowSize)`.
#' - For `filterFourierTransformArtefacts`: `numeric(1)` defining the m/z window
#'   left and right of a peak where to remove fourier transform artefacts.
#'
#' @param i For `[`: `integer`, `logical` or `character` to subset the object.
#'
#' @param j For `[`: not supported.
#'
#' @param initial For `tic`: `logical(1)` whether the initially
#'     reported total ion current should be reported, or whether the
#'     total ion current should be (re)calculated on the actual data
#'     (`initial = FALSE`, same as `ionCount`).
#'
#' @param intensity For `filterIntensity`: `numeric` of length 1 or 2 defining
#'     either the lower or the lower and upper intensity limit for the
#'     filtering, or a `function` that takes the intensities as input and
#'     returns a `logical` (same length then peaks in the spectrum) whether the
#'     peak should be retained or not. Defaults to `intensity = c(0, Inf)` thus
#'     only peaks with `NA` intensity are removed.
#'
#' @param isotopeTolerance For `filterFourierTransformArtefacts`: the m/z
#'     `tolerance` to be used to define whether peaks might be isotopes of
#'     the current tested peak.
#'
#' @param k For `pickPeaks`: `integer(1)`, number of values left and right of
#'  the peak that should be considered in the weighted mean calculation.
#'
#' @param keep For `filterMzValues` and `filterMzRange`: `logical(1)` whether
#'     the matching peaks should be retained (`keep = TRUE`, the default`)
#'     or dropped (`keep = FALSE`).
#'
#' @param keepIsotopes For `filterFourierTransformArtefacts`: whether isotope
#'     peaks should not be removed as fourier artefacts.
#'
#' @param maxCharge For `filterFourierTransformArtefacts`: the maximum charge
#'     to be considered for isotopes.
#'
#' @param MAPFUN For `compareSpectra`: function to map/match peaks between the
#'     two compared spectra. See [joinPeaks()] for more information and possible
#'     functions.
#'
#' @param method
#' - For `pickPeaks`: `character(1)`, the noise estimators that
#'   should be used, currently the the *M*edian *A*bsolute *D*eviation
#'   (`method = "MAD"`) and Friedman's Super Smoother
#'   (`method = "SuperSmoother"`) are supported.
#' - For `smooth`: `character(1)`, the smoothing function that should be used,
#'   currently, the Moving-Average- (`method = "MovingAverage"`),
#'   Weighted-Moving-Average- (`method = "WeightedMovingAverage")`,
#'   Savitzky-Golay-Smoothing (`method = "SavitzkyGolay"`) are supported.
#' - For `estimatePrecursorIntensity`: `character(1)` defining whether the
#'   precursor intensity should be estimated on the previous MS1 spectrum
#'   (`method = "previous"`, the default) or based on an interpolation on the
#'   previous and next MS1 spectrum (`method = "interpolation"`).
#'
#' @param metadata For `Spectra`: optional `list` with metadata information.
#'
#' @param msLevel. `integer` defining the MS level(s) of the spectra to which
#'     the function should be applied (defaults to all MS levels of `object`.
#'     For `filterMsLevel`: the MS level to which `object` should be subsetted.
#'
#' @param mz For `filterIsolationWindow`: `numeric(1)` with the m/z value to
#'     filter the object. For `filterPrecursorMz` and `filterMzRange`:
#'     `numeric(2)` defining the lower and upper m/z boundary.
#'     For `filterMzValues` and `filterPrecursorMzValues`: `numeric` with the
#'     m/z values to match peaks or precursor m/z against.
#'
#' @param n for `filterAcquisitionNum`: `integer` with the acquisition numbers
#'     to filter for.
#'
#' @param name For `$` and `$<-`: the name of the spectra variable to return
#'     or set.
#'
#' @param neutralLoss for `containsNeutralLoss`: `numeric(1)` defining the value
#'     which should be subtracted from the spectrum's precursor m/z.
#'
#' @param object For `Spectra`: either a `DataFrame` or `missing`. See section
#'     on creation of `Spectra` objects for details. For all other methods a
#'     `Spectra` object.
#'
#' @param p For `combineSpectra`: `factor` defining how to split the input
#'     `Spectra` for parallel processing. Defaults to `x$dataStorage`, i.e.,
#'     depending on the used backend, per-file parallel processing will be
#'     performed.
#'
#' @param polarity for `filterPolarity`: `integer` specifying the polarity to
#'     to subset `object`.
#'
#' @param ppm For `compareSpectra`, `containsMz`, `deisotopeSpectra`,
#'     `filterMzValues` and `reduceSpectra`: `numeric(1)`
#'     defining a relative, m/z-dependent, maximal accepted difference between
#'     m/z values for peaks to be matched (or grouped).
#'
#' @param processingQueue For `Spectra`: optional `list` of
#'     [ProcessingStep-class] objects.
#'
#' @param rt for `filterRt`: `numeric(2)` defining the retention time range to
#'     be used to subset/filter `object`.
#'
#' @param SIMPLIFY For `compareSpectra` whether the result matrix should be
#'     *simplified* to a `numeric` if possible (i.e. if either `x` or `y` is
#'     of length 1).
#'
#' @param snr For `pickPeaks`: `double(1)` defining the
#'     *S*ignal-to-*N*oise-*R*atio. The intensity of a local maximum has to be
#'     higher than `snr * noise` to be considered as peak.
#'
#' @param source For `Spectra`: instance of [MsBackend-class] that can be used
#'     to import spectrum data from the provided files. See section *Creation
#'     of objects, conversion and changing the backend* for more details.
#'
#' @param spectraVariables For `selectSpectraVariables`: `character` with the
#'     names of the spectra variables to which the backend should be subsetted.
#'     For `addProcessing`: `character` with additional spectra variables that
#'     should be passed along to the function defined with `FUN`. See function
#'     description for details.
#'
#' @param substDefinition For `deisotopeSpectra`: `matrix` or `data.frame`
#'     with definitions of isotopic substitutions. Uses by default isotopic
#'     substitutions defined from all compounds in the Human Metabolome
#'     Database (HMDB). See [isotopologues()] or [isotopicSubstitutionMatrix()]
#'     for details.
#'
#' @param suffix.y A `character(1)` specifying the suffix to be used
#'     for making the names of columns in the merged spectra variables
#'     unique. This suffix will be used to amend `names(y)`, while
#'     `spectraVariables(x)` will remain unchanged.
#'
#' @param tolerance For `compareSpectra`, `containsMz`, `deisotopeSpectra`,
#'     `filterMzValues` and `reduceSpectra`: `numeric(1)` allowing to define
#'     a constant maximal accepted difference between m/z values for peaks
#'     to be matched (or grouped). For `containsMz` it can also be of length
#'     equal `mz` to specify a different tolerance for each m/z value.
#'
#' @param threshold
#' - For `pickPeaks`: a `double(1)` defining the proportion of the maximal peak
#'      intensity. Just values above are used for the weighted mean calculation.
#' - For `replaceIntensitiesBelow`: a `numeric(1)` defining the threshold or
#'   a `function` to calculate the threshold for each spectrum on its intensity
#'   values. Defaults to `threshold = min`.
#' - For `filterFourierTransformArtefacts`: the relative intensity (to a peak)
#'   below which peaks are considered fourier artefacts. Defaults to
#'   `threshold = 0.2` hence removing peaks that have an intensity below 0.2
#'   times the intensity of the tested peak (within the selected
#'   `halfWindowSize`).
#'
#' @param use.names For `lengths`: ignored.
#'
#' @param value replacement value for `<-` methods. See individual
#'     method description or expected data type.
#'
#' @param which for `containsMz`: either `"any"` or `"all"` defining whether any
#'     (the default) or all provided `mz` have to be present in the spectrum.
#'
#' @param x A `Spectra` object.
#'
#' @param y A `Spectra` object. A `DataFrame` for `joinSpectraData()`.
#'
#' @param z For `filterPrecursorCharge`: `integer()` with the precursor charges
#'     to be used as filter.
#'
#' @param ... Additional arguments.
#'
#' @author Sebastian Gibb, Johannes Rainer, Laurent Gatto
#'
#' @md
#'
#' @exportClass Spectra
#'
#' @exportMethod Spectra
#'
#' @examples
#'
#' ## Create a Spectra providing a `DataFrame` containing the spectrum data.
#'
#' spd <- DataFrame(msLevel = c(1L, 2L), rtime = c(1.1, 1.2))
#' spd$mz <- list(c(100, 103.2, 104.3, 106.5), c(45.6, 120.4, 190.2))
#' spd$intensity <- list(c(200, 400, 34.2, 17), c(12.3, 15.2, 6.8))
#'
#' data <- Spectra(spd)
#' data
#'
#' ## Get the number of spectra
#' length(data)
#'
#' ## Get the number of peaks per spectrum
#' lengths(data)
#'
#' ## Create a Spectra from mzML files and use the `MsBackendMzR` on-disk
#' ## backend.
#' sciex_file <- dir(system.file("sciex", package = "msdata"),
#'     full.names = TRUE)
#' sciex <- Spectra(sciex_file, backend = MsBackendMzR())
#' sciex
#'
#' ## The MS data is on disk and will be read into memory on-demand. We can
#' ## however change the backend to a MsBackendMemory backend which will
#' ## keep all of the data in memory.
#' sciex_im <- setBackend(sciex, MsBackendMemory())
#' sciex_im
#'
#' ## The `MsBackendMemory()` supports the `setBackend` method:
#' supportsSetBackend(MsBackendMemory())
#'
#' ## Thus, it is possible to change to that backend with `setBackend`. Most
#' ## read-only backends however don't support that, such as the
#' ## `MsBackendMzR` and `setBackend` would fail to change to that backend.
#' supportsSetBackend(MsBackendMzR())
#'
#' ## The on-disk object `sciex` is light-weight, because it does not keep the
#' ## MS peak data in memory. The `sciex_im` object in contrast keeps all the
#' ## data in memory and its size is thus much larger.
#' object.size(sciex)
#' object.size(sciex_im)
#'
#' ## The spectra variable `dataStorage` returns for each spectrum the location
#' ## where the data is stored. For in-memory objects:
#' head(dataStorage(sciex_im))
#'
#' ## While objects that use an on-disk backend will list the files where the
#' ## data is stored.
#' head(dataStorage(sciex))
#'
#' ## The spectra variable `dataOrigin` returns for each spectrum the *origin*
#' ## of the data. If the data is read from e.g. mzML files, this will be the
#' ## original mzML file name:
#' head(dataOrigin(sciex))
#' head(dataOrigin(sciex_im))
#'
#' ## ---- ACCESSING AND ADDING DATA ----
#'
#' ## Get the MS level for each spectrum.
#' msLevel(data)
#'
#' ## Alternatively, we could also use $ to access a specific spectra variable.
#' ## This could also be used to add additional spectra variables to the
#' ## object (see further below).
#' data$msLevel
#'
#' ## Get the intensity and m/z values.
#' intensity(data)
#' mz(data)
#'
#' ## Determine whether one of the spectra has a specific m/z value
#' containsMz(data, mz = 120.4)
#'
#' ## Accessing spectra variables works for all backends:
#' intensity(sciex)
#' intensity(sciex_im)
#'
#' ## Get the m/z for the first spectrum.
#' mz(data)[[1]]
#'
#' ## Get the peak data (m/z and intensity values).
#' pks <- peaksData(data)
#' pks
#' pks[[1]]
#' pks[[2]]
#'
#' ## Note that we could get the same resulb by coercing the `Spectra` to
#' ## a `list` or `SimpleList`:
#' as(data, "list")
#' as(data, "SimpleList")
#'
#' ## List all available spectra variables (i.e. spectrum data and metadata).
#' spectraVariables(data)
#'
#' ## For all *core* spectrum variables accessor functions are available. These
#' ## return NA if the variable was not set.
#' centroided(data)
#' dataStorage(data)
#' rtime(data)
#' precursorMz(data)
#'
#' ## The core spectra variables are:
#' coreSpectraVariables()
#'
#' ## Add an additional metadata column.
#' data$spectrum_id <- c("sp_1", "sp_2")
#'
#' ## List spectra variables, "spectrum_id" is now also listed
#' spectraVariables(data)
#'
#' ## Get the values for the new spectra variable
#' data$spectrum_id
#'
#' ## Extract specific spectra variables.
#' spectraData(data, columns = c("spectrum_id", "msLevel"))
#'
#' ## Drop spectra variable data and/or columns.
#' res <- selectSpectraVariables(data, c("mz", "intensity"))
#'
#' ## This removed the additional columns "spectrum_id" and deleted all values
#' ## for all spectra variables, except "mz" and "intensity".
#' spectraData(res)
#'
#' ## Compared to the data before selectSpectraVariables.
#' spectraData(data)
#'
#'
#' ## ---- SUBSETTING, FILTERING AND COMBINING
#'
#' ## Subset to all MS2 spectra.
#' data[msLevel(data) == 2]
#'
#' ## Same with the filterMsLevel function
#' filterMsLevel(data, 2)
#'
#' ## Below we combine the `data` and `sciex_im` objects into a single one.
#' data_comb <- c(data, sciex_im)
#'
#' ## The combined Spectra contains a union of all spectra variables:
#' head(data_comb$spectrum_id)
#' head(data_comb$rtime)
#' head(data_comb$dataStorage)
#' head(data_comb$dataOrigin)
#'
#' ## Filter a Spectra for a target precursor m/z with a tolerance of 10ppm
#' spd$precursorMz <- c(323.4, 543.2302)
#' data_filt <- Spectra(spd)
#' filterPrecursorMzRange(data_filt, mz = 543.23 + ppm(c(-543.23, 543.23), 10))
#'
#' ## Filter a Spectra keeping only peaks matching certain m/z values
#' sps_sub <- filterMzValues(data, mz = c(103, 104), tolerance = 0.3)
#' mz(sps_sub)
#'
#' ## This function can also be used to remove specific peaks from a spectrum
#' ## by setting `keep = FALSE`.
#' sps_sub <- filterMzValues(data, mz = c(103, 104),
#'     tolerance = 0.3, keep = FALSE)
#' mz(sps_sub)
#'
#' ## Note that `filterMzValues` keeps or removes all peaks with a matching
#' ## m/z given the provided `ppm` and `tolerance` parameters.
#'
#' ## Filter a Spectra keeping only peaks within a m/z range
#' sps_sub <- filterMzRange(data, mz = c(100, 300))
#' mz(sps_sub)
#'
#' ## Remove empty spectra variables
#' sciex_noNA <- dropNaSpectraVariables(sciex)
#'
#' ## Available spectra variables before and after dropNaSpectraVariables
#' spectraVariables(sciex)
#' spectraVariables(sciex_noNA)
#'
#'
#' ## Adding new spectra variables
#' sciex1 <- filterDataOrigin(sciex, dataOrigin(sciex)[1])
#' spv <- DataFrame(spectrumId = sciex1$spectrumId[3:12], ## used for merging
#'                  var1 = rnorm(10),
#'                  var2 = sample(letters, 10))
#' spv
#'
#' sciex2 <- joinSpectraData(sciex1, spv, by.y = "spectrumId")
#'
#' spectraVariables(sciex2)
#' spectraData(sciex2)[1:13, c("spectrumId", "var1", "var2")]
#'
#' ## Removing fourier transform artefacts seen in Orbitra data.
#'
#' ## Loading an Orbitrap spectrum with artefacts.
#' data(fft_spectrum)
#' plotSpectra(fft_spectrum, xlim = c(264.5, 265.5))
#' plotSpectra(fft_spectrum, xlim = c(264.5, 265.5), ylim = c(0, 5e6))
#'
#' fft_spectrum <- filterFourierTransformArtefacts(fft_spectrum)
#' fft_spectrum
#' plotSpectra(fft_spectrum, xlim = c(264.5, 265.5), ylim = c(0, 5e6))
#'
#' ## Using a few examples peaks in your data you can optimize the parameters
#' fft_spectrum_filtered <- filterFourierTransformArtefacts(fft_spectrum,
#'                                                halfWindowSize = 0.2,
#'                                                threshold = 0.005,
#'                                                keepIsotopes = TRUE,
#'                                                maxCharge = 5,
#'                                                isotopeTolerance = 0.005
#'                                                )
#'
#' fft_spectrum_filtered
#' length(mz(fft_spectrum_filtered)[[1]])
#' plotSpectra(fft_spectrum_filtered, xlim = c(264.5, 265.5), ylim = c(0, 5e6))
#'
#' ## ---- DATA MANIPULATIONS AND OTHER OPERATIONS ----
#'
#' ## Set the data to be centroided
#' centroided(data) <- TRUE
#'
#' ## Replace peak intensities below 40 with 3.
#' res <- replaceIntensitiesBelow(data, threshold = 40, value = 3)
#' res
#'
#' ## Get the intensities of the first and second spectrum.
#' intensity(res)[[1]]
#' intensity(res)[[2]]
#'
#' ## Remove all peaks with an intensity below 40.
#' res <- filterIntensity(res, intensity = c(40, Inf))
#'
#' ## Get the intensities of the first and second spectrum.
#' intensity(res)[[1]]
#' intensity(res)[[2]]
#'
#' ## Lengths of spectra is now different
#' lengths(mz(res))
#' lengths(mz(data))
#'
#' ## In addition it is possible to pass a function to `filterIntensity`: in
#' ## the example below we want to keep only peaks that have an intensity which
#' ## is larger than one third of the maximal peak intensity in that spectrum.
#' keep_peaks <- function(x, prop = 3) {
#'     x > max(x, na.rm = TRUE) / prop
#' }
#' res2 <- filterIntensity(data, intensity = keep_peaks)
#' intensity(res2)[[1L]]
#' intensity(data)[[1L]]
#'
#' ## We can also change the proportion by simply passing the `prop` parameter
#' ## to the function. To keep only peaks that have an intensity which is
#' ## larger than half of the maximum intensity:
#' res2 <- filterIntensity(data, intensity = keep_peaks, prop = 2)
#' intensity(res2)[[1L]]
#' intensity(data)[[1L]]
#'
#' ## Since data manipulation operations are by default not directly applied to
#' ## the data but only added to the internal lazy evaluation queue, it is also
#' ## possible to remove these data manipulations with the `reset` function:
#' res_rest <- reset(res)
#' res_rest
#' lengths(mz(res_rest))
#' lengths(mz(res))
#' lengths(mz(data))
#'
#' ## `reset` after a `applyProcessing` can not restore the data, because the
#' ## data in the backend was changed. Similarly, `reset` after any filter
#' ## operations can not restore data for a `Spectra` with a
#' ## `MsBackendMemory` or `MsBackendDataFrame`.
#' res_2 <- applyProcessing(res)
#' res_rest <- reset(res_2)
#' lengths(mz(res))
#' lengths(mz(res_rest))
#'
#'
#' ## Compare spectra: comparing spectra 2 and 3 against spectra 10:20 using
#' ## the normalized dotproduct method.
#' res <- compareSpectra(sciex_im[2:3], sciex_im[10:20])
#' ## first row contains comparisons of spectrum 2 with spectra 10 to 20 and
#' ## the second row comparisons of spectrum 3 with spectra 10 to 20
#' res
#'
#' ## To use a simple Pearson correlation instead we can define a function
#' ## that takes the two peak matrices and calculates the correlation for
#' ## their second columns (containing the intensity values).
#' correlateSpectra <- function(x, y, use = "pairwise.complete.obs", ...) {
#'     cor(x[, 2], y[, 2], use = use)
#' }
#' res <- compareSpectra(sciex_im[2:3], sciex_im[10:20],
#'     FUN = correlateSpectra)
#' res
#'
#' ## Use compareSpectra to determine the number of common (matching) peaks
#' ## with a ppm of 10:
#' ## type = "inner" uses a *inner join* to match peaks, i.e. keeps only
#' ## peaks that can be mapped betwen both spectra. The provided FUN returns
#' ## simply the number of matching peaks.
#' compareSpectra(sciex_im[2:3], sciex_im[10:20], ppm = 10, type = "inner",
#'     FUN = function(x, y, ...) nrow(x))
#'
#' ## Apply an arbitrary function to each spectrum in a Spectra.
#' ## In the example below we calculate the mean intensity for each spectrum
#' ## in a subset of the sciex_im data. Note that we can access all variables
#' ## of each individual spectrum either with the `$` operator or the
#' ## corresponding method.
#' res <- spectrapply(sciex_im[1:20], FUN = function(x) mean(x$intensity[[1]]))
#' head(res)
#'
#' ## It is however important to note that dedicated methods to access the
#' ## data (such as `intensity`) are much more efficient than using `lapply`:
#' res <- lapply(intensity(sciex_im[1:20]), mean)
#' head(res)
#'
#' ## As an alternative, applying a function `FUN` to a `Spectra` can be
#' ## performed *chunk-wise*. The advantage of this is, that only the data for
#' ## one chunk at a time needs to be loaded into memory reducing the memory
#' ## demand. This type of processing can be performed by specifying the size
#' ## of the chunks (i.e. number of spectra per chunk) with the `chunkSize`
#' ## parameter
#' spectrapply(sciex_im[1:20], lengths, chunkSize = 5L)
#'
#' ## Calculating the precursor intensity for MS2 spectra:
#' ##
#' ## Some MS instrument manufacturer don't report the precursor intensities
#' ## for MS2 spectra. The `estimatePrecursorIntensity` function can be used
#' ## in these cases to calculate the precursor intensity on MS1 data. Below
#' ## we load an mzML file from a vendor providing precursor intensities and
#' ## compare the estimated and reported precursor intensities.
#' tmt <- Spectra(msdata::proteomics(full.names = TRUE)[5],
#'     backend = MsBackendMzR())
#' pmi <- estimatePrecursorIntensity(tmt)
#' plot(pmi, precursorIntensity(tmt))
#'
#' ## We can also replace the original precursor intensity values with the
#' ## newly calculated ones
#' tmt$precursorIntensity <- pmi
#'
#' ## ---- DATA EXPORT ----
#'
#' ## Some `MsBackend` classes provide an `export` method to export the data to
#' ## the file format supported by the backend. The `MsBackendMzR` for example
#' ## allows to export MS data to mzML or mzXML file(s), the `MsBackendMgf`
#' ## (defined in the MsBackendMgf R package) would allow to export the data
#' ## in mgf file format. Below we export the MS data in `data`. We
#' ## call the `export` method on this object, specify the backend that should
#' ## be used to export the data (and which also defines the output format) and
#' ## provide a file name.
#' fl <- tempfile()
#' export(data, MsBackendMzR(), file = fl)
#'
#' ## This exported our data in mzML format. Below we read the first 6 lines
#' ## from that file.
#' readLines(fl, n = 6)
#'
#' ## If only a single file name is provided, all spectra are exported to that
#' ## file. To export data with the `MsBackendMzR` backend to different files, a
#' ## file name for each individual spectrum has to be provided.
#' ## Below we export each spectrum to its own file.
#' fls <- c(tempfile(), tempfile())
#' export(data, MsBackendMzR(), file = fls)
#'
#' ## Reading the data from the first file
#' res <- Spectra(backendInitialize(MsBackendMzR(), fls[1]))
#'
#' mz(res)
#' mz(data)
NULL

#' The Spectra class
#'
#' The [Spectra-class] encapsulates data and meta-data for mass
#' spectrometry experiments.
#'
#'
#' @slot backend A derivate of [MsBackend-class] holding/controlling the spectra
#' data.
#' @slot processingQueue `list` of `ProcessingStep` objects.
#' @slot processingQueueVariables `character` of spectraVariables that should
#'     be passed to the processing step function.
#' @slot processing A `character` storing logging information.
#' @slot metadata A `list` storing experiment metadata.
#' @slot version A `characher(1)` containing the class version.
#'
#' @name Spectra-class
#' @docType class
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#'
#' @importClassesFrom S4Vectors DataFrame
#'
#' @importMethodsFrom S4Vectors lapply
#'
#' @importFrom S4Vectors DataFrame
#'
#' @noRd
setClass(
    "Spectra",
    slots = c(
        backend = "MsBackend",
        processingQueue = "list",
        processingQueueVariables = "character",
        ## logging
        processing = "character",
        ## metadata
        metadata = "list",
        version = "character"
    ),
    prototype = prototype(version = "0.2")
)

setValidity("Spectra", function(object) {
    msg <- .valid_processing_queue(object@processingQueue)
    if (length(msg)) msg
    else TRUE
})

#' @rdname hidden_aliases
#'
#' @importMethodsFrom methods show
#'
#' @importFrom utils capture.output
#'
#' @exportMethod show
setMethod("show", "Spectra",
    function(object) {
        cat("MSn data (", class(object)[1L], ") with ",
            length(object@backend), " spectra in a ", class(object@backend),
            " backend:\n", sep = "")
        if (length(object@backend)) {
            txt <- capture.output(show(object@backend))
            cat(txt[-1], sep = "\n")
        }
        if (length(object@processingQueue))
            cat("Lazy evaluation queue:", length(object@processingQueue),
                "processing step(s)\n")
        lp <- length(object@processing)
        if (lp) {
            lps <- object@processing
            if (lp > 3) {
                lps <- lps[1:3]
                lps <- c(lps, paste0("...", lp - 3, " more processings. ",
                                     "Use 'processingLog' to list all."))
            }
            cat("Processing:\n", paste(lps, collapse="\n "), "\n")
        }
    })

#' @rdname Spectra
setMethod("Spectra", "missing", function(object, processingQueue = list(),
                                         metadata = list(), ...,
                                         backend = MsBackendMemory(),
                                         BPPARAM = bpparam()) {
    new("Spectra", metadata = metadata, processingQueue = processingQueue,
        backend = backend)
})

#' @rdname Spectra
setMethod("Spectra", "MsBackend", function(object, processingQueue = list(),
                                           metadata = list(), ...,
                                           BPPARAM = bpparam()) {
    new("Spectra", metadata = metadata, processingQueue = processingQueue,
        backend = object)
})

#' @rdname Spectra
#'
#' @importFrom methods callNextMethod
setMethod("Spectra", "character", function(object, processingQueue = list(),
                                           metadata = list(),
                                           source = MsBackendMzR(),
                                           backend = source,
                                           ..., BPPARAM = bpparam()) {
    if (!length(object))
        Spectra(backend, metadata = metadata,
                processingQueue = processingQueue)
    else
        callNextMethod(object = object, processingQueue = processingQueue,
                       metadata = metadata, source = source, backend = backend,
                       ..., BPPARAM = BPPARAM)
})

#' @rdname Spectra
setMethod("Spectra", "ANY", function(object, processingQueue = list(),
                                     metadata = list(),
                                     source = MsBackendMemory(),
                                     backend = source,
                                     ..., BPPARAM = bpparam()) {
    sp <- new("Spectra", metadata = metadata, processingQueue = processingQueue,
              backend = backendInitialize(
                  source, object, ...,
                  BPPARAM = backendBpparam(source, BPPARAM)))
    if (class(source)[1L] != class(backend)[1L])
        setBackend(sp, backend, ..., BPPARAM = backendBpparam(backend, BPPARAM))
    else sp
})

#' @rdname Spectra
#'
#' @exportMethod setBackend
setMethod(
    "setBackend", c("Spectra", "MsBackend"),
    function(object, backend, f = dataStorage(object), ...,
             BPPARAM = bpparam()) {
        backend_class <- class(object@backend)
        BPPARAM <- backendBpparam(object@backend, BPPARAM)
        if (!supportsSetBackend(backend))
            stop(class(backend), " does not support 'setBackend'")
        if (!length(object)) {
            bknds <- backendInitialize(
                backend, data = spectraData(object@backend), ...)
        } else {
            f <- force(factor(f, levels = unique(f)))
            if (length(f) != length(object))
                stop("length of 'f' has to match the length of 'object'")
            if (length(levels(f)) == 1L) {
                bknds <- backendInitialize(
                    backend, data = spectraData(object@backend), ...)
            } else {
                bknds <- bplapply(
                    split(object@backend, f = f),
                    function(z, ...) {
                        backendInitialize(backend,
                                          data = spectraData(z), ...,
                                          BPPARAM = SerialParam())
                    }, ..., BPPARAM = BPPARAM)
                bknds <- backendMerge(bknds)
                ## That below ensures the backend is returned in its original
                ## order - unsplit does unfortunately not work.
                if (is.unsorted(f))
                    bknds <- bknds[order(unlist(split(seq_along(bknds), f),
                                                use.names = FALSE))]
            }
        }
        object@backend <- bknds
        object@processing <- .logging(object@processing,
                                      "Switch backend from ",
                                      backend_class, " to ",
                                      class(object@backend))
        object
    })

#' @rdname Spectra
#'
#' @importFrom MsCoreUtils vapply1c
#'
#' @exportMethod c
setMethod("c", "Spectra", function(x, ...) {
    .concatenate_spectra(unname(list(unname(x), ...)))
})

#' @rdname Spectra
setMethod("split", "Spectra", function(x, f, drop = FALSE, ...) {
    bcknds <- split(x@backend, f, ...)
    lapply(bcknds, function(b) {
        slot(x, "backend", check = FALSE) <- b
        x
    })
})

#' @rdname Spectra
#'
#' @export
setMethod("export", "Spectra",
          function(object, backend, ...) {
              if (missing(backend))
                  stop("Parameter 'backend' is required.")
              export(backend, object, ...)
          })

#### ---------------------------------------------------------------------------
##
##                          ACCESSOR METHODS
##
#### ---------------------------------------------------------------------------

#' @rdname Spectra
setMethod("acquisitionNum", "Spectra", function(object)
    acquisitionNum(object@backend))

#' @rdname Spectra
setMethod(
    "peaksData", "Spectra",
    function(object, columns = c("mz", "intensity"), ..., BPPARAM = bpparam()) {
        BPPARAM <- backendBpparam(object, BPPARAM)
        if (is(BPPARAM, "SerialParam"))
            f <- NULL
        else f <- factor(dataStorage(object))
        SimpleList(.peaksapply(object, f = f, columns = columns,
                               BPPARAM = BPPARAM, ...))
})

#' @rdname Spectra
setMethod("peaksVariables", "Spectra", function(object)
    peaksVariables(object@backend))

#' @importFrom methods setAs
setAs("Spectra", "list", function(from, to) {
    BPPARAM <- backendBpparam(from)
    if (is(BPPARAM, "SerialParam"))
        f <- NULL
    else f <- factor(dataStorage(from))
    .peaksapply(from, f = f, BPPARAM = BPPARAM)
})

setAs("Spectra", "SimpleList", function(from, to) {
    peaksData(from)
})

#' @rdname Spectra
setMethod("centroided", "Spectra", function(object) {
    centroided(object@backend)
})

#' @rdname Spectra
setReplaceMethod("centroided", "Spectra", function(object, value) {
    centroided(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("collisionEnergy", "Spectra", function(object) {
    collisionEnergy(object@backend)
})

#' @rdname Spectra
setReplaceMethod("collisionEnergy", "Spectra", function(object, value) {
    collisionEnergy(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("dataOrigin", "Spectra", function(object) dataOrigin(object@backend))

#' @rdname Spectra
setReplaceMethod("dataOrigin", "Spectra", function(object, value) {
    dataOrigin(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("dataStorage", "Spectra",
          function(object) dataStorage(object@backend))

#' @rdname Spectra
setMethod("dropNaSpectraVariables", "Spectra", function(object) {
    object@backend <- dropNaSpectraVariables(object@backend)
    object
})

#' @rdname Spectra
setMethod("intensity", "Spectra", function(object, ...) {
    if (length(object@processingQueue))
        NumericList(.peaksapply(object, FUN = function(z, ...) z[, 2], ...),
                    compress = FALSE)
    else intensity(object@backend) # Disables also parallel proc. (issue #44)
})

#' @rdname Spectra
setMethod("ionCount", "Spectra", function(object) {
    BPPARAM <- backendBpparam(object)
    if (is(BPPARAM, "SerialParam"))
        f <- NULL
    else f <- factor(dataStorage(object))
    if (length(object))
        unlist(.peaksapply(object, f = f, BPPARAM = BPPARAM,
                           FUN = function(pks, ...)
                               sum(pks[, 2], na.rm = TRUE)),
               use.names = FALSE)
    else numeric()
})

#' @rdname Spectra
setMethod("isCentroided", "Spectra", function(object, ...) {
    BPPARAM <- backendBpparam(object)
    if (is(BPPARAM, "SerialParam"))
        f <- NULL
    else f <- factor(dataStorage(object))
    if (length(object))
        unlist(.peaksapply(object, f = f, BPPARAM = BPPARAM,
                           FUN = .peaks_is_centroided),
               use.names = FALSE)
    else logical()
})

#' @rdname Spectra
setMethod("isEmpty", "Spectra", function(x) {
    BPPARAM <- backendBpparam(x)
    if (is(BPPARAM, "SerialParam"))
        f <- NULL
    else f <- factor(dataStorage(x))
    if (length(x))
        unlist(.peaksapply(x, f = f, BPPARAM = BPPARAM,
                           FUN = function(pks, ...) nrow(pks) == 0),
               use.names = FALSE)
    else logical()
})

#' @rdname Spectra
setMethod("isolationWindowLowerMz", "Spectra", function(object) {
    isolationWindowLowerMz(object@backend)
})

#' @rdname Spectra
setReplaceMethod("isolationWindowLowerMz", "Spectra", function(object, value) {
    isolationWindowLowerMz(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("isolationWindowTargetMz", "Spectra", function(object) {
    isolationWindowTargetMz(object@backend)
})

#' @rdname Spectra
setReplaceMethod("isolationWindowTargetMz", "Spectra", function(object, value) {
    isolationWindowTargetMz(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("isolationWindowUpperMz", "Spectra", function(object) {
    isolationWindowUpperMz(object@backend)
})

#' @rdname Spectra
setReplaceMethod("isolationWindowUpperMz", "Spectra", function(object, value) {
    isolationWindowUpperMz(object@backend) <- value
    object
})

#' @rdname Spectra
#'
#' @exportMethod containsMz
setMethod("containsMz", "Spectra", function(object, mz = numeric(),
                                            tolerance = 0,
                                            ppm = 20, which = c("any", "all"),
                                            BPPARAM = bpparam()) {
    cond_fun <- match.fun(match.arg(which))
    if (all(is.na(mz)))
        return(rep(NA, length(object)))
    mz <- unique(sort(mz))
    BPPARAM <- backendBpparam(object@backend, BPPARAM)
    if (is(BPPARAM, "SerialParam"))
        .has_mz(object, mz, tolerance = tolerance, ppm = ppm,
                condFun = cond_fun, parallel = BPPARAM)
    else {
        sp <- SerialParam()
        f <- as.factor(dataStorage(object))
        res <- .lapply(object, FUN = .has_mz, mz = mz, tolerance = tolerance,
                       condFun = cond_fun, parallel = sp, f = f,
                       BPPARAM = BPPARAM)
        unsplit(res, f = f)
    }
})

#' @rdname Spectra
#'
#' @exportMethod containsNeutralLoss
setMethod("containsNeutralLoss", "Spectra", function(object, neutralLoss = 0,
                                                     tolerance = 0, ppm = 20,
                                                     BPPARAM = bpparam()) {
    BPPARAM <- backendBpparam(object@backend, BPPARAM)
    if (is(BPPARAM, "SerialParam")) {
        .has_mz_each(object, precursorMz(object) - neutralLoss,
                     tolerance = tolerance, ppm = ppm, parallel = BPPARAM)
    } else {
        sp <- SerialParam()
        f <- as.factor(dataStorage(object))
        res <- .lapply(object, FUN = function(obj, n, tol, ppm, par) {
            .has_mz_each(obj, precursorMz(obj) - n, tolerance = tol,
                         ppm = ppm, parallel = sp)
        }, n = neutralLoss, tol = tolerance, ppm = ppm, par = sp, f = f,
                       BPPARAM = BPPARAM)
        unsplit(res, f = f)
    }
})

#' @rdname Spectra
#'
#' @importMethodsFrom ProtGenerics spectrapply
#'
#' @exportMethod spectrapply
setMethod("spectrapply", "Spectra", function(object, FUN, ...,
                                             chunkSize = integer(),
                                             f = factor(),
                                             BPPARAM = SerialParam()) {
    if (missing(FUN))
        FUN <- identity
    if (length(chunkSize))
        return(chunkapply(object, FUN, ..., chunkSize = chunkSize))
    if (!length(f))
        f <- as.factor(seq_along(object))
    .lapply(object, FUN = FUN, f = f, ...,
            BPPARAM = backendBpparam(object@backend, BPPARAM))
})

#' @rdname Spectra
#'
#' @exportMethod length
setMethod("length", "Spectra", function(x) length(x@backend))

#' @rdname Spectra
setMethod("msLevel", "Spectra", function(object) msLevel(object@backend))

#' @rdname Spectra
setMethod("mz", "Spectra", function(object, ...) {
    if (length(object@processingQueue))
        NumericList(.peaksapply(object, FUN = function(z, ...) z[, 1], ...),
                    compress = FALSE)
    else mz(object@backend) # Disables also parallel processing (issue #44)
})

#' @rdname Spectra
#'
#' @exportMethod lengths
setMethod("lengths", "Spectra", function(x, use.names = FALSE) {
    if (length(x)) {
        if (length(x@processingQueue))
            unlist(.peaksapply(x, FUN = function(pks, ...) nrow(pks)),
                   use.names = use.names)
        else lengths(x@backend, use.names = use.names)
    } else integer()
})

#' @rdname Spectra
setMethod("polarity", "Spectra", function(object) {
    polarity(object@backend)
})

#' @rdname Spectra
setReplaceMethod("polarity", "Spectra", function(object, value) {
    polarity(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("precScanNum", "Spectra", function(object) {
    precScanNum(object@backend)
})

#' @rdname Spectra
setMethod("precursorCharge", "Spectra", function(object) {
    precursorCharge(object@backend)
})

#' @rdname Spectra
setMethod("precursorIntensity", "Spectra", function(object) {
    precursorIntensity(object@backend)
})

#' @rdname Spectra
setMethod("precursorMz", "Spectra", function(object) {
    precursorMz(object@backend)
})

#' @rdname Spectra
setMethod("rtime", "Spectra", function(object) {
    rtime(object@backend)
})

#' @rdname Spectra
setReplaceMethod("rtime", "Spectra", function(object, value) {
    rtime(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("scanIndex", "Spectra", function(object) {
    scanIndex(object@backend)
})

#' @rdname Spectra
setMethod(
    "selectSpectraVariables", "Spectra",
    function(object, spectraVariables = union(spectraVariables(object),
                                              peaksVariables(object))) {
        spectraVariables <- union(spectraVariables, "dataStorage")
        object@backend <- selectSpectraVariables(
            object@backend, spectraVariables = spectraVariables)
        object
    })

#' @rdname Spectra
setMethod("smoothed", "Spectra", function(object) {
    smoothed(object@backend)
})

#' @rdname Spectra
setReplaceMethod("smoothed", "Spectra", function(object, value) {
    smoothed(object@backend) <- value
    object
})

#' @rdname Spectra
#'
#' @importMethodsFrom ProtGenerics spectraData
#'
#' @exportMethod spectraData
setMethod(
    "spectraData", "Spectra",
    function(object, columns = spectraVariables(object)) {
        if (length(object@processingQueue) &&
            length(pcns <- intersect(columns, peaksVariables(object)))) {
            scns <- setdiff(columns, pcns)
            if (length(scns))
                spd <- spectraData(object@backend, columns = scns)
            else
                spd <- make_zero_col_DFrame(nrow = length(object))
            pkd <- peaksData(object, columns = pcns)
            for (pcn in pcns) {
                vals <- lapply(pkd, `[`, , pcn)
                if (pcn %in% c("mz", "intensity"))
                    vals <- NumericList(vals, compress = FALSE)
                spd <- do.call(`[[<-`, list(spd, i = pcn, value = vals))
            }
            spd
        } else
            spectraData(object@backend, columns = columns)
    })

#' @rdname Spectra
#'
#' @importMethodsFrom ProtGenerics spectraData<-
#'
#' @exportMethod spectraData<-
setReplaceMethod("spectraData", "Spectra", function(object, value) {
    if (!inherits(value, "DataFrame"))
        stop("'spectraData<-' expects a 'DataFrame' as input.", call. = FALSE)
    pvs <- peaksVariables(object)
    if (length(object@processingQueue) &&
        any(colnames(value) %in% pvs))
        stop("Can not replace peaks variables with a non-empty processing ",
             "queue. Please use 'object <- applyProcessing(object)' to apply ",
             "and clear the processing queue. Note that 'applyProcessing' ",
             "requires a *writeable* backend. Use e.g. 'object <- ",
             "setBackend(object, MsBackendMemory())' if needed.")
    pvs <- setdiff(pvs, colnames(value))
    if (length(pvs)) {
        sd <- spectraData(object, pvs)
        for (pv in pvs) {
            value <- do.call("$<-", list(value, name = pv, sd[, pv]))
        }
        object@processingQueue <- list()
    }
    spectraData(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("spectraNames", "Spectra", function(object) {
    spectraNames(object@backend)
})

#' @rdname Spectra
setReplaceMethod("spectraNames", "Spectra", function(object, value) {
    spectraNames(object@backend) <- value
    object
})

#' @rdname Spectra
setMethod("spectraVariables", "Spectra", function(object) {
    setdiff(spectraVariables(object@backend), peaksVariables(object@backend))
})

#' @rdname Spectra
setMethod("tic", "Spectra", function(object, initial = TRUE) {
    if (!length(object))
        return(numeric())
    if (initial)
        tic(object@backend, initial = initial)
    else ionCount(object)
})

#' @rdname Spectra
#'
#' @importMethodsFrom S4Vectors $
#'
#' @export
setMethod("$", "Spectra", function(x, name) {
    if (!(name %in% c(spectraVariables(x@backend), peaksVariables(x@backend))))
        stop("No spectra variable '", name, "' available")
    if (name == "mz")
        mz(x)
    else if (name == "intensity")
        intensity(x)
    else {
        if (length(x@processingQueue) && name %in% peaksVariables(x))
            .peaksapply(x, FUN = function(z, ...) z[, name],
                        columns = c("mz", "intensity", name))
        else
            do.call("$", list(x@backend, name))
    }
})

#' @rdname Spectra
#'
#' @export
setReplaceMethod("$", "Spectra", function(x, name, value) {
    if (length(x@processingQueue) &&
        any(name %in% peaksVariables(x)))
        stop("Can not replace peaks variables with a non-empty processing ",
             "queue. Please use 'object <- applyProcessing(object)' to apply ",
             "and clear the processing queue. Note that 'applyProcessing' ",
             "requires a *writeable* backend. Use e.g. 'object <- ",
             "setBackend(object, MsBackendMemory())' if needed.")
    x@backend <- do.call("$<-", list(x@backend, name, value))
    x
})

#' @rdname Spectra
#'
#' @export
setMethod("[[", "Spectra", function(x, i, j, ...) {
    if (!is.character(i))
        stop("'i' is supposed to be a character defining the spectra ",
             "variable to access.")
    if (!missing(j))
        stop("'j' is not supported.")
    if (!(i %in% c(spectraVariables(x), "mz", "intensity")))
        stop("No spectra variable '", i, "' available")
    if (i == "mz")
        mz(x)
    else if (i == "intensity")
        intensity(x)
    else
        do.call("[[", list(x@backend, i))
})

#' @rdname Spectra
#'
#' @export
setReplaceMethod("[[", "Spectra", function(x, i, j, ..., value) {
    if (!is.character(i))
        stop("'i' is supposed to be a character defining the spectra ",
             "variable to replace or create.")
    if (!missing(j))
        stop("'j' is not supported.")
    x@backend <- do.call("[[<-", list(x@backend, i = i, value = value))
    x
})

#### ---------------------------------------------------------------------------
##
##                      FILTERING AND SUBSETTING
##
#### ---------------------------------------------------------------------------

#' @rdname Spectra
setMethod("[", "Spectra", function(x, i, j, ..., drop = FALSE) {
    if (!missing(j))
        stop("Subsetting 'Spectra' by columns is not (yet) supported")
    if (missing(i))
        return(x)
    slot(x, "backend", check = FALSE) <- x@backend[i = i]
    x
})

#' @rdname Spectra
setMethod("filterAcquisitionNum", "Spectra", function(object, n = integer(),
                                                      dataStorage = character(),
                                                      dataOrigin = character()) {
    if (length(dataStorage) && !is.character(dataStorage))
        stop("'dataStorage' is expected to be of type character")
    if (length(dataOrigin) && !is.character(dataOrigin))
        stop("'dataOrigin' is expected to be of type character")
    object@backend <- filterAcquisitionNum(object@backend, n,
                                           dataStorage, dataOrigin)
    object@processing <- .logging(object@processing,
                                  "Filter: select by: ", length(n),
                                  " acquisition number(s) in ",
                                  max(length(dataStorage), length(dataOrigin)),
                                  " file(s)")
    object
})

#' @rdname Spectra
setMethod("filterEmptySpectra", "Spectra", function(object) {
    object@backend <- object@backend[as.logical(lengths(object))]
    object@processing <- .logging(object@processing,
                                  "Filter: removed empty spectra.")
    object
})

#' @rdname Spectra
setMethod("filterDataOrigin", "Spectra", function(object,
                                                  dataOrigin = character()) {
    if (length(dataOrigin) && !is.character(dataOrigin))
        stop("'dataOrigin' is expected to be of type character")
    object@backend <- filterDataOrigin(object@backend, dataOrigin = dataOrigin)
    object@processing <- .logging(object@processing,
                                  "Filter: select data origin(s) ",
                                  paste0(dataOrigin, collapse = ", "))
    object
})

#' @rdname Spectra
setMethod("filterDataStorage", "Spectra", function(object,
                                                   dataStorage = character()) {
    if (length(dataStorage) && !is.character(dataStorage))
        stop("'dataStorage' is expected to be of type character")
    object@backend <- filterDataStorage(object@backend, dataStorage)
    object@processing <- .logging(object@processing,
                                  "Filter: select data storage(s) ",
                                  paste0(dataStorage, collapse = ", "))
    object
})

#' @rdname Spectra
#'
#' @exportMethod filterFourierTransformArtefacts
setMethod("filterFourierTransformArtefacts", "Spectra",
          function(object, halfWindowSize = 0.05, threshold = 0.2,
                   keepIsotopes = TRUE, maxCharge = 5,
                   isotopeTolerance = 0.005) {
              object <- addProcessing(object, .peaks_remove_fft_artifact,
                                      halfWindowSize = halfWindowSize,
                                      threshold = threshold,
                                      keepIsotopes = keepIsotopes,
                                      maxCharge = maxCharge,
                                      isotopeTolerance = isotopeTolerance)
              object@processing <- .logging(
                  object@processing, "Remove fast fourier artefacts.")
              object
          })

#' @rdname Spectra
#'
#' @importMethodsFrom ProtGenerics filterIntensity
#'
#' @exportMethod filterIntensity
setMethod("filterIntensity", "Spectra",
          function(object, intensity = c(0, Inf),
                   msLevel. = uniqueMsLevels(object), ...) {
              if (!.check_ms_level(object, msLevel.))
                  return(object)
              if (is.numeric(intensity)) {
                  if (length(intensity) == 1)
                      intensity <- c(intensity, Inf)
                  if (length(intensity) != 2)
                      stop("'intensity' should be of length specifying a ",
                           "lower intensity limit or of length two defining ",
                           "a lower and upper limit.")
                  object <- addProcessing(object, .peaks_filter_intensity,
                                          intensity = intensity,
                                          msLevel = msLevel.,
                                          spectraVariables = "msLevel")
                  object@processing <- .logging(
                      object@processing, "Remove peaks with intensities ",
                      "outside [", intensity[1], ", ", intensity[2],
                      "] in spectra of MS level(s) ",
                      paste0(msLevel., collapse = ", "), ".")
              } else {
                  if (is.function(intensity)) {
                      object <- addProcessing(
                          object, .peaks_filter_intensity_function,
                          intfun = intensity, msLevel = msLevel.,
                          args = list(...), spectraVariables = "msLevel")
                      object@processing <- .logging(
                          object@processing, "Remove peaks based on their ",
                          "intensities and a user-provided function ",
                          "in spectra of MS level(s) ",
                          paste0(msLevel., collapse = ", "), ".")
                  }
                  else stop("'intensity' has to be numeric or a function")
              }
              object
          })


#' @rdname Spectra
setMethod("filterIsolationWindow", "Spectra", function(object, mz = numeric()) {
    object@backend <- filterIsolationWindow(object@backend, mz = mz)
    object@processing <- .logging(object@processing,
                                  "Filter: select spectra containing m/z ",
                                  mz, " in their isolation window")
    object
})

#' @rdname Spectra
setMethod("filterMsLevel", "Spectra", function(object, msLevel. = integer()) {
    object@backend <- filterMsLevel(object@backend, msLevel = msLevel.)
    object@processing <- .logging(object@processing,
                                  "Filter: select MS level(s) ",
                                  paste0(unique(msLevel.), collapse = " "))
    object
})

#' @rdname Spectra
#'
#' @export
setMethod("filterMzRange", "Spectra",
          function(object, mz = numeric(), msLevel. = uniqueMsLevels(object),
                   keep = TRUE) {
              if (!.check_ms_level(object, msLevel.))
                  return(object)
              if (!length(mz)) mz <- c(-Inf, Inf)
              else mz <- range(mz)
              object <- addProcessing(object, .peaks_filter_mz_range, mz = mz,
                                      msLevel = msLevel., keep = keep,
                                      spectraVariables = "msLevel")
              if (keep) keep_or_remove <- "select"
              else keep_or_remove <- "remove"
              object@processing <- .logging(
                  object@processing, "Filter: ", keep_or_remove,
                  " peaks with an m/z within [", mz[1L], ", ", mz[2L], "]")
              object
          })

#' @rdname Spectra
#'
#' @export
setMethod("filterMzValues", "Spectra",
          function(object, mz = numeric(), tolerance = 0, ppm = 20,
                   msLevel. = uniqueMsLevels(object), keep = TRUE) {
              if (!.check_ms_level(object, msLevel.))
                  return(object)
              l <- length(mz)
              if (length(tolerance) != 1)
                  stop("'tolerance' should be of length 1")
              if (length(ppm) != 1)
                  stop("'ppm' should be of length 1")
              if (is.unsorted(mz)) {
                  idx <- order(mz)
                  mz <- mz[idx]
                  if (length(tolerance) == l)
                      tolerance <- tolerance[idx]
                  if (length(ppm) == l)
                      ppm <- ppm[idx]
              }
              object <- addProcessing(object, .peaks_filter_mz_value,
                                      mz = mz, tolerance = tolerance,
                                      ppm = ppm, msLevel = msLevel.,
                                      keep = keep, spectraVariables = "msLevel")
              if (length(mz) <= 3)
                  what <- paste0(format(mz, digits = 4), collapse = ", ")
              else what <- ""
              if (keep)
                  keep_or_remove <- "select"
              else keep_or_remove <- "remove"
              object@processing <- .logging(
                  object@processing, "Filter: ", keep_or_remove,
                  " peaks matching provided m/z values ", what)
              object
          })

#' @rdname Spectra
setMethod("filterPolarity", "Spectra", function(object, polarity = integer()) {
    object@backend <- filterPolarity(object@backend, polarity = polarity)
    object@processing <- .logging(object@processing,
                                  "Filter: select spectra with polarity ",
                                  paste0(polarity, collapse = " "))
    object
})

#' @rdname Spectra
#'
#' @export
setMethod("filterPrecursorMz", "Spectra",
          function(object, mz = numeric()) {
              .Deprecated(
                  msg = paste0("'filterPrecursorMz' is deprecated. Please use",
                               " 'filterPrecursorMzRange' instead."))
              object@backend <- filterPrecursorMzRange(object@backend, mz)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select spectra with a precursor m/z within [",
                  paste0(mz, collapse = ", "), "]")
              object
          })

#' @rdname Spectra
setMethod("filterPrecursorMzRange", "Spectra",
          function(object, mz = numeric()) {
              object@backend <- filterPrecursorMzRange(object@backend, mz)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select spectra with a precursor m/z within [",
                  paste0(mz, collapse = ", "), "]")
              object
          })

#' @rdname Spectra
setMethod("filterPrecursorMzValues", "Spectra",
          function(object, mz = numeric(), ppm = 20, tolerance = 0) {
              object@backend <- filterPrecursorMzValues(
                  object@backend, sort(mz), ppm = ppm, tolerance = tolerance)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select spectra with precursor m/z matching ",
                  paste0(mz, collapse = ", "), "")
              object
          })

#' @rdname Spectra
setMethod("filterPrecursorCharge", "Spectra",
          function(object, z = integer()) {
              z <- unique(z)
              object@backend <- filterPrecursorCharge(object@backend, z)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select spectra with a precursor charge ",
                  paste0(z, collapse = ", "))
              object
          })

#' @rdname Spectra
setMethod("filterPrecursorScan", "Spectra",
          function(object, acquisitionNum = integer(), f = dataOrigin(object)) {
              if (!all(f %in% unique(dataOrigin(object))))
                  stop("'f' must be in dataOrigin().")
              object@backend <- filterPrecursorScan(object@backend,
                                                    acquisitionNum,
                                                    f = dataOrigin(object))
              object@backend <- filterDataOrigin(object@backend, dataOrigin = f)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select parent/children scans for ",
                  paste0(acquisitionNum, collapse = " "))
              object
          })

#' @rdname Spectra
setMethod("filterRt", "Spectra",
          function(object, rt = numeric(), msLevel. = uniqueMsLevels(object)) {
              if (!is.numeric(msLevel.))
                  stop("Please provide a numeric MS level.")
              if (length(rt) != 2L || !is.numeric(rt) || rt[1] >= rt[2])
                  stop("Please provide a lower and upper numeric retention",
                       " time range.")
              if (length(rt))
                  rt <- range(rt)
              else rt <- c(-Inf, Inf)
              object@backend <- filterRt(object@backend, rt, msLevel.)
              object@processing <- .logging(
                  object@processing,
                  "Filter: select retention time [", rt[1], "..", rt[2],
                  "] on MS level(s) ", paste0(msLevel., collapse = " "))
              object
          })

#' @rdname Spectra
setMethod("reset", "Spectra", function(object, ...) {
    object@backend <- reset(object@backend)
    object@processingQueue <- list()
    if (!.hasSlot(object, "processingQueueVariables"))
        object <- updateObject(object, check = FALSE)
    object@processingQueueVariables <- character()
    object@processing <- .logging(object@processing, "Reset object.")
    object
})

#### ---------------------------------------------------------------------------
##
##                      DATA MANIPULATION METHODS
##
#### ---------------------------------------------------------------------------

#' @rdname Spectra
#'
#' @importMethodsFrom ProtGenerics bin
#'
#' @exportMethod bin
setMethod("bin", "Spectra", function(x, binSize = 1L, breaks = NULL,
                                     msLevel. = uniqueMsLevels(x),
                                     FUN = sum) {
    if (!.check_ms_level(x, msLevel.))
        return(x)
    if (!length(breaks)) {
        mzr <- range(.peaksapply(filterMsLevel(x, msLevel.),
                                 function(z, ...) z[c(1L, nrow(z))]
                                 ), na.rm = TRUE)
        breaks <- seq(floor(mzr[1]), ceiling(mzr[2]), by = binSize)
        breaks <- MsCoreUtils:::.fix_breaks(breaks, mzr)
    }
    mids <- (breaks[-length(breaks)] + breaks[-1L]) / 2
    x <- addProcessing(x, .peaks_bin, breaks = breaks, mids = mids,
                       agg_fun = FUN, msLevel = msLevel.,
                       spectraVariables = "msLevel")
    x@processing <- .logging(x@processing,
                             "Spectra of MS level(s) ",
                             paste0(msLevel., collapse = ", "),
                             " binned.")
    x
})

#' @rdname Spectra
#'
#' @exportMethod compareSpectra
#'
#' @importFrom MsCoreUtils ndotproduct
#'
#' @importMethodsFrom ProtGenerics compareSpectra
#'
#' @exportMethod compareSpectra
setMethod("compareSpectra", signature(x = "Spectra", y = "Spectra"),
          function(x, y, MAPFUN = joinPeaks, tolerance = 0, ppm = 20,
                   FUN = ndotproduct, ..., SIMPLIFY = TRUE) {
              mat <- .compare_spectra_chunk(x, y, MAPFUN = MAPFUN,
                                            tolerance = tolerance,
                                            ppm = ppm, FUN = FUN, ...)
              if (SIMPLIFY && (length(x) == 1 || length(y) == 1))
                  mat <- as.vector(mat)
              mat
          })
#' @rdname Spectra
setMethod("compareSpectra", signature(x = "Spectra", y = "missing"),
          function(x, y = NULL, MAPFUN = joinPeaks, tolerance = 0, ppm = 20,
                   FUN = ndotproduct, ..., SIMPLIFY = TRUE) {
              if (length(x) == 1)
                  return(compareSpectra(x, x, MAPFUN = MAPFUN,
                                        tolerance = tolerance,
                                        ppm = ppm, FUN = FUN, ...,
                                        SIMPLIFY = SIMPLIFY))
              mat <- .compare_spectra_self(x, MAPFUN = MAPFUN, FUN = FUN,
                                           tolerance = tolerance, ppm = ppm,
                                           ...)
              if (SIMPLIFY && length(x) == 1)
                  mat <- as.vector(mat)
              mat
          })

## estimateMzResolution

## estimateNoise

## normalize

#' @rdname Spectra
#'
#' @exportMethod pickPeaks
setMethod("pickPeaks", "Spectra",
          function(object, halfWindowSize = 2L,
                   method = c("MAD", "SuperSmoother"), snr = 0, k = 0L,
                   descending = FALSE, threshold = 0,
                   msLevel. = uniqueMsLevels(object), ...) {
    if (!.check_ms_level(object, msLevel.))
        return(object)
    if (!is.integer(halfWindowSize) || length(halfWindowSize) != 1L ||
        halfWindowSize <= 0L)
        stop("Argument 'halfWindowSize' has to be an integer of length 1 ",
             "and > 0.")
    if (!is.numeric(snr) || length(snr) != 1L || snr < 0L)
        stop("Argument 'snr' has to be a numeric of length 1 that is >= 0.")
    if (!is.integer(k) || length(k) != 1L || k < 0L)
        stop("Argument 'k' has to be an integer of length 1 that is >= 0.")
    if (!is.logical(descending) || length(descending) != 1L ||
        is.na(descending))
        stop("Argument 'descending' has to be just TRUE or FALSE")
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        threshold < 0L || threshold > 1L)
        stop("Argument 'threshold' has to be a numeric of length 1 ",
             "that is >= 0 and <= 1.")

    method <- match.arg(method)

    object <- addProcessing(object, .peaks_pick,
                            halfWindowSize = halfWindowSize, method = method,
                            snr = snr, k = k, descending = descending,
                            threshold = threshold, msLevel = msLevel., ...,
                            spectraVariables = c("msLevel", "centroided"))
    object$centroided[msLevel(object) %in% msLevel.] <- TRUE
    object@processing <- .logging(object@processing,
                                  "Peak picking with ", method,
                                  " noise estimation, hws = ", halfWindowSize,
                                  ", snr = ", snr,
                                  if (k > 0) " and centroid refinement")
    object
})

## quantify

## removeReporters

#' @rdname Spectra
#'
#' @exportMethod replaceIntensitiesBelow
setMethod("replaceIntensitiesBelow", "Spectra",
          function(object, threshold = min, value = 0,
                   msLevel. = uniqueMsLevels(object)) {
              if (!is.numeric(threshold) && !is.function(threshold))
                  stop("Argument 'threshold' has to be either numeric or ",
                       "a function.")
              if (!.check_ms_level(object, msLevel.))
                  return(object)
              object <- addProcessing(
                  object, .peaks_replace_intensity, threshold = threshold,
                  value = value, msLevel = msLevel.,
                  spectraVariables = c("msLevel", "centroided"))
              msg <- ifelse(
                  is.function(threshold),
                  yes = "a threshold defined by a provided function",
                  no = threshold)
              object@processing <- .logging(object@processing,
                                            "Signal <= ", msg,
                                            " in MS level(s) ",
                                            paste0(msLevel., collapse = ", "),
                                            " set to 0")
              object
          })


#' @rdname Spectra
#'
#' @importFrom ProtGenerics smooth
#' @importFrom MsCoreUtils coefMA coefWMA coefSG
#' @exportMethod smooth
setMethod("smooth", "Spectra",
          function(x, halfWindowSize = 2L,
                   method = c("MovingAverage", "WeightedMovingAverage",
                              "SavitzkyGolay"),
                   msLevel. = uniqueMsLevels(x), ...) {
    if (!.check_ms_level(x, msLevel.))
        return(x)
    if (!is.integer(halfWindowSize) || length(halfWindowSize) != 1L ||
        halfWindowSize <= 0L)
        stop("Argument 'halfWindowSize' has to be an integer of length 1 ",
             "and > 0.")

    method <- match.arg(method)
    coef <- switch(method,
                   MovingAverage = coefMA(halfWindowSize),
                   WeightedMovingAverage = coefWMA(halfWindowSize),
                   SavitzkyGolay = coefSG(halfWindowSize, ...))

    x <- addProcessing(x, .peaks_smooth, halfWindowSize = halfWindowSize,
                       coef = coef, msLevel = msLevel., ...,
                       spectraVariables = "msLevel")
    x@processing <- .logging(x@processing, "Spectra smoothing with ", method,
                                           ", hws = ", halfWindowSize)
    x
})

#' @exportMethod addProcessing
#'
#' @importFrom ProtGenerics ProcessingStep
#'
#' @importMethodsFrom ProtGenerics addProcessing
#'
#' @importClassesFrom ProtGenerics ProcessingStep
#'
#' @importFrom methods .hasSlot
#'
#' @importFrom BiocGenerics updateObject
#'
#' @rdname Spectra
setMethod("addProcessing", "Spectra", function(object, FUN, ...,
                                               spectraVariables = character()) {
    if (missing(FUN))
        return(object)
    object@processingQueue <- c(object@processingQueue,
                                list(ProcessingStep(FUN, ARGS = list(...))))
    if (!.hasSlot(object, "processingQueueVariables"))
        object <- updateObject(object)
    object@processingQueueVariables <- union(object@processingQueueVariables,
                                             spectraVariables)
    validObject(object)
    object
})

#' @rdname Spectra
#'
#' @export
coreSpectraVariables <- function() .SPECTRA_DATA_COLUMNS

#' @rdname Spectra
setMethod("uniqueMsLevels", "Spectra", function(object, ...) {
    uniqueMsLevels(object@backend, ...)
})

#' @rdname Spectra
setMethod("backendBpparam", "Spectra", function(object, BPPARAM = bpparam()) {
    backendBpparam(object@backend, BPPARAM)
})
