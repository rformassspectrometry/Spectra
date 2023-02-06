#' @include hidden_aliases.R

#' @title Mass spectrometry data backends
#'
#' @aliases class:MsBackend MsBackend-class MsBackendDataFrame-class
#' @aliases MsBackendMzR-class [,MsBackend-method
#' @aliases uniqueMsLevels,MsBackend-method
#' @aliases MsBackendMemory-class
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
#' - `MsBackendMemory` and `MsBackendDataFrame`: store all data in memory. The
#'   `MsBackendMemory` is optimized for accessing and processing the peak data
#'   (i.e. the numerical matrices with the m/z and intensity values) while the
#'   `MsBackendDataFrame` keeps all data in a `DataFrame`.
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
#'     For `peaksData` accessor: optional `character` with requested columns in
#'     the individual `matrix` of the returned `list`. Defaults to
#'     `peaksVariables(object)` and depends on what *peaks variables* the
#'     backend provides.
#'
#' @param data For `backendInitialize`: `DataFrame` with spectrum
#'     metadata/data. This parameter can be empty for `MsBackendMzR` backends
#'     but needs to be provided for `MsBackendDataFrame` backends.
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
#' @param f `factor` defining the grouping to split `x`. See [split()]. For
#'     `filterPrecursorScan`: factor defining from which original data files
#'     the spectra derive to avoid selecting spectra from different
#'     samples/files. Defaults to `f = dataOrigin(object)`.
#'
#' @param file For `filterFile`: index or name of the file(s) to which the data
#'     should be subsetted. For `export`: `character` of length 1 or equal to
#'     the number of spectra.
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
#' @param msLevel. same as `msLevel` above.
#'
#' @param mz For `filterIsolationWindow`: `numeric(1)` with the m/z value to
#'     filter the object. For `filterPrecursorMzRange`: `numeric(2)` with the
#'     lower and upper m/z boundary. For `filterPrecursorMzValues`: `numeric`
#'     with the m/z value(s) to filter the object.
#'
#' @param peaksVariables For `backendInitialize` for `MsBackendMemory`:
#'     `character` specifying which of the columns of the provided `data`
#'     contain *peaks variables* (i.e. information for individual mass
#'     peaks). Defaults to `peaksVariables = c("mz", "intensity")`. `"mz"`
#'     and `"intensity"` should **always** be specified.
#'
#' @param ppm For `filterPrecursorMzValues`: `numeric(1)` with the m/z-relative
#'     maximal acceptable difference for a m/z to be considered matching. See
#'     [closest()] for details.
#'
#' @param z For `filterPrecursorCharge`: `integer()` with the precursor charges
#'     to be used as filter.
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
#' @param spectraVariables For `selectSpectraVariables`: `character` with the
#'     names of the spectra variables to which the backend should be subsetted.
#'
#' @param tolerance For `filterPrecursorMzValues`: `numeric(1)` with the
#'     maximal absolute acceptable difference for a m/z value to be considered
#'     matching. See [closest()] for details.
#'
#' @param use.names For `lengths`: whether spectrum names should be used.
#'
#' @param value replacement value for `<-` methods. See individual
#'     method description or expected data type.
#'
#' @param x Object extending `MsBackend`.
#'
#' @param ... Additional arguments.
#'
#'
#' @section Implementation notes:
#'
#' Backends extending `MsBackend` **must** implement all of its methods (listed
#' above). Developers of new `MsBackend`s should follow the
#' `MsBackendDataFrame` implementation. To ensure a new implementation being
#' conform with the `MsBackend` definition, developers should included test
#' suites provided by this package in their unit test setup. For that a variable
#' `be` should be created in the package's `"testthat.R"` file that represents
#' a (initialized) instance of the developed backend. Then the path to the
#' test suites should be defined with
#' `test_suite <- system.file("test_backends", "test_MsBackend",
#' package = "Spectra")` followed by `test_dir(test_suite)` to run all test
#' files in that directory. Individual unit test files could be run with
#' `test_file(file.path(test_suite, "test_spectra_variables.R"),
#' stop_on_failure = TRUE)` (note that without `stop_on_failure = TRUE` tests
#' would fail silently) . Adding this code to the packages `"testthat.R"` file
#' ensures that all tests checking the validity of an `MsBackend` instance
#' defined in the `Spectra` package are also run on the newly develped backend
#' class.
#'
#' The `MsBackend` defines the following slots:
#'
#' - `@readonly`: `logical(1)` whether the backend supports writing/replacing
#'   of m/z or intensity values.
#'
#' @section Backend functions:
#'
#' New backend classes **must** extend the base `MsBackend` class and
#' **have** to implement the following methods:
#'
#' - `[`: subset the backend. Only subsetting by element (*row*/`i`) is
#'   allowed. Parameter `i` should support `integer` indices and `logical`
#'   and should throw an error if `i` is out of bounds. The
#'   `MsCoreUtils::i2index` could be used to check the input `i`.
#'   For `i = integer()` an empty backend should be returned.
#'
#' - `$`, `$<-`: access or set/add a single spectrum variable (column) in the
#'   backend. Using a `value` of `NULL` should allow deleting the specified
#'   spectra variable. An error should be thrown if the spectra variable is not
#'   available.
#'
#' - `[[`, `[[<-`: access or set/add a single spectrum variable (column) in the
#'   backend. The default implementation uses `$`, thus these methods don't have
#'   to be implemented for new classes extending `MsBackend`.
#'
#' - `acquisitionNum`: returns the acquisition number of each
#'   spectrum. Returns an `integer` of length equal to the number of
#'   spectra (with `NA_integer_` if not available).
#'
#' - `backendInitialize`: initialises the backend. This method is
#'   supposed to be called rights after creating an instance of the
#'   backend class and should prepare the backend (e.g. set the data
#'   for the memory backend or read the spectra header data for the
#'   `MsBackendMzR` backend). Parameters can be defined freely for each
#'   backend, depending on what is needed to initialize the backend. It
#'   is however suggested to also support a parameter `data` that can be
#'   used to submit the full spectra data as a `DataFrame` to the
#'   backend. This would allow the backend to be also usable for the
#'   [setBackend()] function from `Spectra`.
#'   The `backendInitialize` method has also to ensure to correctly set
#'   spectra variable `dataStorage`.
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
#' - `dropNaSpectraVariables`: removes spectra variables (i.e. columns in the
#'   object's `spectraData` that contain only missing values (`NA`). Note that
#'   while columns with only `NA`s are removed, a `spectraData` call after
#'   `dropNaSpectraVariables` might still show columns containing `NA` values
#'   for *core* spectra variables.
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
#' - `export`: exports data from a `Spectra` class to a file. This method is
#'   called by the `export,Spectra` method that passes itself as a second
#'   argument to the function. The `export,MsBackend` implementation is thus
#'   expected to take a `Spectra` class as second argument from which all data
#'   is exported. Taking data from a `Spectra` class ensures that also all
#'   eventual data manipulations (cached in the `Spectra`'s lazy evaluation
#'   queue) are applied prior to export - this would not be possible with only a
#'   [MsBackend] class. An example implementation is the `export` method
#'   for the `MsBackendMzR` backend that supports export of the data in
#'   *mzML* or *mzXML* format. See the documentation for the `MsBackendMzR`
#'   class below for more information.
#'
#' - `filterAcquisitionNum`: filters the object keeping only spectra matching
#'   the provided acquisition numbers (argument `n`). If `dataOrigin` or
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
#'   Implementation of this method is optional since a default implementation
#'   for `MsBackend` is available.
#'
#' - `filterDataStorage`: filters the object retaining spectra matching the
#'   provided `dataStorage`. Parameter `dataStorage` has to be of type
#'   `character` and needs to match exactly the data storage value of the
#'   spectra to subset.
#'   `filterDataStorage` should return the data ordered by the provided
#'   `dataStorage` parameter, i.e. if `dataStorage = c("2", "1")` was provided,
#'   the spectra in the resulting object should be ordered accordingly (first
#'   spectra from data storage `"2"` and then from `"1"`).
#'   Implementation of this method is optional since a default implementation
#'   for `MsBackend` is available.
#'
#' - `filterEmptySpectra`: removes empty spectra (i.e. spectra without peaks).
#'   Implementation of this method is optional since a default implementation
#'   for `MsBackend` is available.
#'
#' - `filterFile`: retains data of files matching the file index or file name
#'    provided with parameter `file`.
#'
#' - `filterIsolationWindow`: retains spectra that contain `mz` in their
#'   isolation window m/z range (i.e. with an `isolationWindowLowerMz` `<=` `mz`
#'   and `isolationWindowUpperMz` `>=` `mz`.
#'   Implementation of this method is optional since a default implementation
#'   for `MsBackend` is available.
#'
#' - `filterMsLevel`: retains spectra of MS level `msLevel`.
#'   Implementation of this method is optional since a default implementation
#'   for `MsBackend` is available.
#'
#' - `filterPolarity`: retains spectra of polarity `polarity`.
#'   Implementation of this method is optional since a default implementation
#'   for `MsBackend` is available.
#'
#' - `filterPrecursorMzRange` (previously `filterPrecursorMz`): retains spectra
#'   with a precursor m/z within the provided m/z range.
#'   Implementation of this method is optional since a default implementation
#'   for `MsBackend` is available.
#'
#' - `filterPrecursorMzValues`: retains spectra with a precursor m/z matching
#'   any of the provided m/z values (given `ppm` and `tolerance`).
#'   Implementation of this method is optional since a default implementation
#'   for `MsBackend` is available.
#'
#' - `filterPrecursorCharge`: retains spectra with the defined precursor
#'   charge(s).
#'   Implementation of this method is optional since a default implementation
#'   for `MsBackend` is available.
#'
#' - `filterPrecursorScan`: retains parent (e.g. MS1) and children scans (e.g.
#'    MS2) of acquisition number `acquisitionNum`. Parameter `f` is supposed to
#'   define the origin of the spectra (i.e. the original data file) to ensure
#'   related spectra from the same file/sample are selected and retained.
#'   Implementation of this method is optional since a default implementation
#'   for `MsBackend` is available.
#'
#' - `filterRt`: retains spectra of MS level `msLevel` with retention times
#'    within (`>=`) `rt[1]` and (`<=`) `rt[2]`.
#'   Implementation of this method is optional since a default implementation
#'   for `MsBackend` is available.
#'
#' - `intensity`: gets the intensity values from the spectra. Returns
#'   a [NumericList()] of `numeric` vectors (intensity values for each
#'   spectrum). The length of the `list` is equal to the number of
#'   `spectra` in `object`.
#'
#' - `intensity<-`: replaces the intensity values. `value` has to be a `list`
#'   (or [NumericList()]) of length equal to the number of spectra and the
#'   number of values within each list element identical to the number of
#'   peaks in each spectrum (i.e. the `lengths(x)`). Note that just
#'   writeable backends support this method.
#'
#' - `ionCount`: returns a `numeric` with the sum of intensities for
#'   each spectrum. If the spectrum is empty (see `isEmpty`),
#'   `NA_real_` is returned.
#'
#' - `isCentroided`: a heuristic approach assessing if the spectra in
#'   `object` are in profile or centroided mode. The function takes
#'   the `qtl` th quantile top peaks, then calculates the difference
#'   between adjacent m/z value and returns `TRUE` if the first
#'   quartile is greater than `k`. (See `Spectra:::.peaks_is_centroided` for
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
#' - `lengths`: gets the number of peaks (m/z-intensity values) per
#'   spectrum.  Returns an `integer` vector (length equal to the
#'   number of spectra). For empty spectra, `0` is returned.
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
#'   `lengths(x)`). Note that just writeable backends support this method.
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
#'   2 and above spectra from the object. Returns a vector of length equal to
#'   the number of spectra in `object`. `NA` are reported for MS1
#'   spectra of if no precursor information is available.
#'
#' - `peaksData` returns a `list` with the spectras' peak data, i.e. numeric
#'   `matrix` with peak values. The length of the list is equal to the number
#'   of spectra in `object`. Each element of the list is a `numeric` `matrix`
#'   with columns depending on the provided `columns` parameter (by default
#'   `"mz"` and `"intensity"`, but depends on the backend's available
#'   `peaksVariables`). For an empty spectrum, a `matrix` with 0 rows and
#'   columns according to `columns` is returned. The optional parameter
#'   `columns`, if supported by the backend, allows to define which peak
#'   variables should be returned in the `numeric` peak `matrix`. As a default
#'   `c("mz", "intensity")` should be used.
#'
#' - `peaksData<-` replaces the peak data (m/z and intensity values) of the
#'   backend. This method expects a `list` of `matrix` objects with columns
#'   `"mz"` and `"intensity"` that has the same length as the number of
#'   spectra in the backend. Note that just writeable backends support this
#'   method.
#'
#' - `peaksVariables`: lists the available variables for mass peaks. Default
#'   peak variables are `"mz"` and `"intensity"` (which all backends need to
#'   support and provide), but some backends might provide additional variables.
#'   These variables correspond to the column names of the `numeric` `matrix`
#'   representing the peak data (returned by `peaksData`).
#'
#' - `reset` a backend (if supported). This method will be called on the backend
#'   by the `reset,Spectra` method that is supposed to restore the data to its
#'   original state (see `reset,Spectra` for more details). The function
#'   returns the *reset* backend. The default implementation for `MsBackend`
#'   returns the backend as-is.
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
#'   the selected spectra variables. It is suggested to **not** remove values
#'   for the `"dataStorage"` variable, since this might be required for some
#'   backends to work properly (such as the `MsBackendMzR`).
#'
#' - `smoothed`,`smoothed<-`: gets or sets whether a spectrum is
#'   *smoothed*. `smoothed` returns a `logical` vector of length equal
#'   to the number of spectra. `smoothed<-` takes a `logical` vector
#'   of length 1 or equal to the number of spectra in `object`.
#'
#' - `spectraData`, `spectraData<-`: gets or sets general spectrum
#'   metadata (annotation, also called header).  `spectraData` returns
#'   a `DataFrame`, `spectraData<-` expects a `DataFrame` with the same number
#'   of rows as there are spectra in `object`. Note that `spectraData` has to
#'   return the full data, i.e. also the m/z and intensity values (as a `list`
#'   or `SimpleList` in columns `"mz"` and `"intensity"`.
#'
#' - `spectraNames`: returns a `character` vector with the names of
#'   the spectra in `object` or `NULL` if not set. `spectraNames<-` allows to
#'   set spectra names (if the object is not read-only).
#'
#' - `spectraVariables`: returns a `character` vector with the
#'   available spectra variables (columns, fields or attributes)
#'   available in `object`. This should return **all** spectra variables which
#'   are present in `object`, also `"mz"` and `"intensity"` (which are by
#'   default not returned by the `spectraVariables,Spectra` method).
#'
#' - `split`: splits the backend into a `list` of backends (depending on
#'   parameter `f`). The default method for `MsBackend` uses [split.default()],
#'   thus backends extending `MsBackend` don't necessarily need to implement
#'   this method.
#'
#' - `tic`: gets the total ion current/count (sum of signal of a
#'   spectrum) for all spectra in `object`. By default, the value
#'   reported in the original raw data file is returned. For an empty
#'   spectrum, `NA_real_` is returned.
#'
#' - `uniqueMsLevels`: gets the unique MS levels of all spectra in `object`.
#'   The default implementation calls `unique(msLevel(object))` but more
#'   efficient implementations could be defined for specific backends.
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
#' @section In-memory data backends: `MsBackendMemory` and `MsBackendDataFrame`:
#'
#' The `MsBackendMemory` and `MsBackendDataFrame` objects keep all MS data in
#' memory are thus ideal for fast data processing. Due to their large memory
#' footprint they are however not suited for large scale experiments. The two
#' backends store the data different. The `MsBackendDataFrame` stores
#' all data in a `DataFrame` and thus supports also S4-classes as
#' spectra variables. Also, sepratate access to m/z or intensity values (i.e.
#' using the `mz` and `intensity` methods) is faster for the
#' `MsBackendDataFrame`. The `MsBackendMemory` on the other hand, due to the
#' way the data is organized internally, provides much faster access to the
#' full peak data (i.e. the numerical matrices of m/z and intensity values).
#' Also subsetting and access to any spectra variable (except `"mz"` and
#' `"intensity"` is fastest for the `MsBackendMemory`. Finally, the
#' `MsBackendMemory` supports also arbitrary peak annotations while the
#' `MsBackendDataFrame` does not have support for such additional peak
#' variables.
#'
#' Thus, for most use cases, the `MsBackendMemory` provides a higher
#' performance and flexibility than the `MsBackendDataFrame` and should thus be
#' preferred. See also issue
#' [246](https://github.com/rformassspectrometry/Spectra/issues/246) for a
#' performance comparison.
#'
#' New objects can be created with the `MsBackendMemory()` and
#' `MsBackendDataFrame()` function, respectively. The backend can be
#' subsequently initialized with the `backendInitialize` method, taking a
#' `DataFrame` (or `data.frame`) with the MS data as first parameter `data`.
#' `backendInitialize` for `MsBackendMemory` has a second parameter
#' `peaksVariables` (default `peaksVariables = c("mz", "intensity")` that
#' allows to specify which of the columns in the provided data frame should
#' be considered as a *peaks variable* (i.e. information of an individual
#' mass peak) rather than a *spectra variable* (i.e. information of an
#' individual spectrum). Note that it is important to also include `"mz"` and
#' `"intensity"` in `peaksVariables` as these would otherwise be considered
#' to be spectra variables! Also, while it is possible to change the values of
#' existing peaks variables using the `$<-` method, this method does **not**
#' allow to add new peaks variables to an existing `MsBackendMemory`. New
#' peaks variables should be added using the `backendInitialize` method.
#'
#' Suggested columns of this `DataFrame` are:
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
#' For the `MsBackendMemory`, any column in the provided `data.frame` which
#' contains a `list` of vectors each with length equal to the number of peaks
#' for a spectrum will be used as additional *peak variable* (see examples
#' below for details).
#'
#' The `MsBackendDataFrame` ignores parameter `columns` of the `peaksData`
#' function and returns **always** m/z and intensity values.
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
#' This backend provides an `export` method to export data from a `Spectra` in
#' *mzML* or *mzXML* format. The definition of the function is:
#'
#' `export(object, x, file = tempfile(), format = c("mzML", "mzXML"),
#'         copy = FALSE)`
#'
#' The parameters are:
#' - `object`: an instance of the `MsBackendMzR` class.
#' - `x`: the [Spectra-class] object to be exported.
#' - `file`: `character` with the (full) output file name(s). Should be
#'   of length 1 or equal `length(x)`. If a single file is specified, all
#'   spectra are exported to that file. Alternatively it is possible to specify
#'   for each spectrum in `x` the name of the file to which it should be
#'   exported (and hence `file` has to be of length equal `length(x)`).
#' - `format`: `character(1)`, either `"mzML"` or `"mzXML"` defining the output
#'   file format.
#' - `copy`: `logical(1)` whether general file information should be copied from
#'   the original MS data files. This only works if `x` uses a `MsBackendMzR`
#'   backend and if `dataOrigin(x)` contains the original MS data file names.
#' - `BPPARAM`: parallel processing settings.
#'
#' See examples in [Spectra-class] or the vignette for more details and
#' examples.
#'
#' The `MsBackendMzR` ignores parameter `columns` of the `peaksData`
#' function and returns **always** m/z and intensity values.
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
#' data across several HDF5 files `data` has to contain a column
#' `"dataStorage"` that defines the grouping of spectra/peaks into files: peaks
#' for spectra with the same value in `"dataStorage"` are saved into the same
#' HDF5 file. If parameter `files` is omitted, the value in `dataStorage` is
#' used as file name (replacing any file ending with `".h5"`. To specify the
#' file names, `files`' length has to match the number of unique elements in
#' `"dataStorage"`.
#'
#' For details see examples on the [Spectra()] help page.
#'
#' The `MsBackendHdf5Peaks` ignores parameter `columns` of the `peaksData`
#' function and returns **always** m/z and intensity values.
#'
#' @section Implementation notes:
#'
#' Backends extending `MsBackend` **must** implement all of its methods (listed
#' above). Developers of new `MsBackend`s should follow the
#' `MsBackendDataFrame` implementation.
#'
#' The [MsBackendCached()] backend provides a caching mechanism to allow
#' *read only* backends to add or change spectra variables. This
#' backend shouldn't be used on its own, but is meant to be extended. See
#' [MsBackendCached()] for details.
#'
#' The `MsBackend` defines the following slots:
#'
#' - `@readonly`: `logical(1)` whether the backend supports writing/replacing
#'   of m/z or intensity values.
#'
#' @name MsBackend
#'
#' @return See documentation of respective function.
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @md
#'
#' @exportClass MsBackend MsBackendDataFrame MsBackendMzR
#'
#' @examples
#'
#' ## The MsBackend class is a virtual class and can not be instantiated
#' ## directly. Below we define a new backend class extending this virtual
#' ## class
#' MsBackendDummy <- setClass("MsBackendDummy", contains = "MsBackend")
#' MsBackendDummy()
#'
#' ## This class inherits now all methods from `MsBackend`, all of which
#' ## however throw an error. These methods would have to be implemented
#' ## for the new backend class.
#' try(mz(MsBackendDummy()))
#'
#' ## See `MsBackendDataFrame` as a reference implementation for a backend
#' ## class (in the *R/MsBackendDataFrame.R* file).
#'
#' ## MsBackendDataFrame
#' ##
#' ## The `MsBackendDataFrame` uses a `S4Vectors::DataFrame` to store all MS
#' ## data. Below we create such a backend by passing a `DataFrame` with all
#' ## data to it.
#' data <- DataFrame(msLevel = c(1L, 2L, 1L), scanIndex = 1:3)
#' data$mz <- list(c(1.1, 1.2, 1.3), c(1.4, 54.2, 56.4, 122.1), c(15.3, 23.2))
#' data$intensity <- list(c(3, 2, 3), c(45, 100, 12.2, 1), c(123, 12324.2))
#'
#' ## Backends are supposed to be created with their specific constructor
#' ## function
#' be <- MsBackendDataFrame()
#'
#' be
#'
#' ## The `backendInitialize` method initializes the backend filling it with
#' ## data. This method can take any parameters needed for the backend to
#' ## get loaded with the data (e.g. a file name from which to load the data,
#' ## a database connection or, in this case, a data frame containing the data).
#' be <- backendInitialize(be, data)
#'
#' be
#'
#' ## Data can be accessed with the accessor methods
#' msLevel(be)
#'
#' mz(be)
#'
#' ## Even if no data was provided for all spectra variables, its accessor
#' ## methods are supposed to return a value.
#' precursorMz(be)
#'
#' ## The `peaksData` method is supposed to return the peaks of the spectra as
#' ## a `list`.
#' peaksData(be)
#'
#' ## List available peaks variables
#' peaksVariables(be)
#'
#' ## Use columns to extract specific peaks variables. Below we extract m/z and
#' ## intensity values, but in reversed order to the default.
#' peaksData(be, columns = c("intensity", "mz"))
#'
#' ## List available spectra variables (i.e. spectrum metadata)
#' spectraVariables(be)
#'
#' ## Extract precursor m/z, rtime, MS level spectra variables
#' spectraData(be, c("precursorMz", "rtime", "msLevel"))
#'
#' ## MsBackendMemory
#' ##
#' ## The `MsBackendMemory` uses a more efficient internal data organization
#' ## and allows also adding arbitrary additional peaks variables (annotations)
#' ## Below we thus add a column "peak_ann" with arbitrary names/ids for each
#' ## peak and add the name of this columnt to the `peaksVariables` parameter
#' ## of the `backendInitialize` method (in addition to `"mz"` and
#' ## `"intensity"` that should **always** be specified.
#' data$peak_ann <- list(c("a", "", "d"), c("", "d", "e", "f"), c("h", "i"))
#' be <- backendInitialize(MsBackendMemory(), data,
#'     peaksVariables = c("mz", "intensity", "peak_ann"))
#' be
#'
#' spectraVariables(be)
#'
#' ## peak_ann is also listed as a peaks variable
#' peaksVariables(be)
#'
#' ## The additional peaks variable can be accessed using the peaksData
#' ## function
#' peaksData(be, "peak_ann")
#'
#' ## The $<- method can be used to replace values of an existing peaks
#' ## variable. It is important that the number of elements matches the
#' ## number of peaks per spectrum.
#' be$peak_ann <- list(1:3, 1:4, 1:2)
#'
#' ## A peaks variable can again be removed by setting it to NULL
#' be$peak_ann <- NULL
#'
#' peaksVariables(be)
NULL

setClass(
    "MsBackend",
    contains = "VIRTUAL",
    slots = c(
        readonly = "logical",
        version = "character"),
    prototype = prototype(readonly = FALSE, version = "0.1"))

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
setMethod("backendInitialize", signature = "MsBackend", function(object, ...) {
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

#' @rdname MsBackend
setMethod("export", "MsBackend", function(object, ...) {
    stop(class(object), " does not support export of data; please provide a ",
         "backend that supports data export with parameter 'backend'.")
})

#' @exportMethod acquisitionNum
#'
#' @importMethodsFrom ProtGenerics acquisitionNum
#'
#' @rdname MsBackend
setMethod("acquisitionNum", "MsBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod peaksData
#'
#' @rdname MsBackend
setMethod("peaksData", "MsBackend", function(object,
                                             columns = c("mz", "intensity")) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod peaksVariables
#'
#' @rdname MsBackend
setMethod("peaksVariables", "MsBackend", function(object) {
    c("mz", "intensity")
})

#' @exportMethod centroided
#'
#' @aliases centroided<-,MsBackend-method
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

#' @exportMethod dropNaSpectraVariables
#'
#' @rdname MsBackend
setMethod("dropNaSpectraVariables", "MsBackend", function(object) {
    svs <- spectraVariables(object)
    svs <- svs[!(svs %in% c("mz", "intensity"))]
    spd <- spectraData(object, columns = svs)
    keep <- !vapply1l(spd, function(z) {
        allna <- all(is.na(z))
        if (length(allna) > 1)
            FALSE
        else allna
    })
    selectSpectraVariables(object, c(svs[keep], "mz", "intensity"))
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
setMethod("filterDataOrigin", "MsBackend",
          function(object, dataOrigin = character()) {
              if (length(dataOrigin)) {
                  object <- object[dataOrigin(object) %in% dataOrigin]
                  if (is.unsorted(dataOrigin))
                      object[order(match(dataOrigin(object), dataOrigin))]
                  else object
              } else object
          })

#' @exportMethod filterDataStorage
#'
#' @importMethodsFrom ProtGenerics filterDataStorage
#'
#' @rdname MsBackend
setMethod("filterDataStorage", "MsBackend",
          function(object, dataStorage = character()) {
              if (length(dataStorage)) {
                  object <- object[dataStorage(object) %in% dataStorage]
                  if (is.unsorted(dataStorage))
                      object[order(match(dataStorage(object), dataStorage))]
                  else object
              } else object
          })

#' @exportMethod filterEmptySpectra
#'
#' @importMethodsFrom ProtGenerics filterEmptySpectra
#'
#' @rdname MsBackend
setMethod("filterEmptySpectra", "MsBackend", function(object, ...) {
    if (!length(object)) return(object)
    object[as.logical(lengths(object))]
})

#' @exportMethod filterIsolationWindow
#'
#' @importMethodsFrom ProtGenerics filterIsolationWindow
#'
#' @rdname MsBackend
setMethod("filterIsolationWindow", "MsBackend",
          function(object, mz = numeric(), ...) {
              if (length(mz)) {
                  if (length(mz) > 1)
                      stop("'mz' is expected to be a single m/z value")
                  keep <- which(isolationWindowLowerMz(object) <= mz &
                                isolationWindowUpperMz(object) >= mz)
                  object[keep]
              } else object
          })

#' @exportMethod filterMsLevel
#'
#' @importMethodsFrom ProtGenerics filterMsLevel
#'
#' @rdname MsBackend
setMethod("filterMsLevel", "MsBackend",
          function(object, msLevel = integer()) {
              if (length(msLevel)) {
                  object[msLevel(object) %in% msLevel]
              } else object
          })

#' @exportMethod filterPolarity
#'
#' @importMethodsFrom ProtGenerics filterPolarity
#'
#' @rdname MsBackend
setMethod("filterPolarity", "MsBackend",
          function(object, polarity = integer()) {
              if (length(polarity))
                  object[polarity(object) %in% polarity]
              else object
          })

#' @exportMethod filterPrecursorMzRange
#'
#' @rdname MsBackend
setMethod("filterPrecursorMzRange", "MsBackend",
          function(object, mz = numeric()) {
              if (length(mz)) {
                  mz <- range(mz)
                  keep <- which(between(precursorMz(object), mz))
                  object[keep]
              } else object
          })

#' @importMethodsFrom ProtGenerics filterPrecursorMz
#'
#' @rdname MsBackend
setMethod("filterPrecursorMz", "MsBackend",
          function(object, mz = numeric()) {
              filterPrecursorMzRange(object, mz)
          })

#' @exportMethod filterPrecursorMzValues
#'
#' @rdname MsBackend
setMethod("filterPrecursorMzValues", "MsBackend",
          function(object, mz = numeric(), ppm = 20, tolerance = 0) {
              if (length(mz)) {
                  object[.values_match_mz(precursorMz(object), mz = mz,
                                          ppm = ppm, tolerance = tolerance)]
              } else object
          })

#' @exportMethod filterPrecursorCharge
#'
#' @importMethodsFrom ProtGenerics filterPrecursorCharge
#'
#' @rdname MsBackend
setMethod("filterPrecursorCharge", "MsBackend",
          function(object, z = integer()) {
              if (length(z)) {
                  keep <- which(precursorCharge(object) %in% z)
                  object[keep]
              } else object
          })

#' @exportMethod filterPrecursorScan
#'
#' @importMethodsFrom ProtGenerics filterPrecursorScan
#'
#' @rdname MsBackend
setMethod("filterPrecursorScan", "MsBackend",
          function(object, acquisitionNum = integer(), f = dataOrigin(object)) {
              if (length(acquisitionNum) && length(f)) {
                  if (!is.factor(f))
                      f <- factor(f, exclude = character())
                  keep <- unsplit(lapply(split(object, f = f), function(z, an) {
                      .filterSpectraHierarchy(acquisitionNum(z),
                                              precScanNum(z),
                                              an)
                  }, an = acquisitionNum), f = f)
                  object[keep]
              } else object
          })

#' @exportMethod filterRt
#'
#' @importMethodsFrom ProtGenerics filterRt
#'
#' @rdname MsBackend
setMethod("filterRt", "MsBackend",
          function(object, rt = numeric(), msLevel. = uniqueMsLevels(object)) {
              if (length(rt)) {
                  rt <- range(rt)
                  sel_ms <- msLevel(object) %in% msLevel.
                  sel_rt <- between(rtime(object), rt) & sel_ms
                  object[sel_rt | !sel_ms]
              } else object
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
#' @importFrom MsCoreUtils vapply1d
#'
#' @rdname MsBackend
setMethod("ionCount", "MsBackend", function(object) {
    vapply1d(intensity(object), sum, na.rm = TRUE)
})

#' @exportMethod isCentroided
#'
#' @importMethodsFrom ProtGenerics isCentroided
#' @importFrom MsCoreUtils vapply1l
#'
#' @rdname MsBackend
setMethod("isCentroided", "MsBackend", function(object, ...) {
    vapply1l(peaksData(object), .peaks_is_centroided)
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
setReplaceMethod("isolationWindowLowerMz", "MsBackend", function(object,
                                                                 value) {
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
setReplaceMethod("isolationWindowTargetMz", "MsBackend", function(object,
                                                                  value) {
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
setReplaceMethod("isolationWindowUpperMz", "MsBackend", function(object,
                                                                 value) {
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

#' @rdname MsBackend
setMethod("lengths", "MsBackend", function(x, use.names = FALSE) {
    stop("Not implemented for ", class(x), ".")
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

#' @exportMethod peaksData<-
#'
#' @rdname MsBackend
setReplaceMethod("peaksData", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod reset
#'
#' @rdname MsBackend
setMethod("reset", "MsBackend", function(object) {
    object
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
setMethod(
    "selectSpectraVariables", "MsBackend",
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
#' @aliases smoothed<-,MsBackend-method
#'
#' @importMethodsFrom ProtGenerics smoothed<-
#'
#' @rdname MsBackend
setReplaceMethod("smoothed", "MsBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod spectraData
#'
#' @rdname MsBackend
setMethod(
    "spectraData", "MsBackend",
    function(object, columns = spectraVariables(object)) {
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

#' @exportMethod split
#'
#' @importMethodsFrom S4Vectors split
#'
#' @rdname MsBackend
setMethod("split", "MsBackend", function(x, f, drop = FALSE, ...) {
    split.default(x, f, drop = drop, ...)
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

#' @exportMethod [[
#'
#' @rdname MsBackend
setMethod("[[", "MsBackend", function(x, i, j, ...) {
    if (!is.character(i))
        stop("'i' is supposed to be a character defining the spectra ",
             "variable to access.")
    if (!missing(j))
        stop("'j' is not supported.")
    do.call("$", list(x, i))
})

#' @exportMethod [[<-
#'
#' @rdname MsBackend
setReplaceMethod("[[", "MsBackend", function(x, i, j, ..., value) {
    if (!is.character(i))
        stop("'i' is supposed to be a character defining the spectra ",
             "variable to replace or create.")
    if (!missing(j))
        stop("'j' is not supported.")
    do.call("$<-", list(x, i, value))
})

#' @exportMethod uniqueMsLevels
#'
#' @rdname MsBackend
setMethod("uniqueMsLevels", "MsBackend", function(object, ...) {
    unique(msLevel(object))
})
