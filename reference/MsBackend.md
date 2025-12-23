# Mass spectrometry data backends

Note that the classes described here are not meant to be used directly
by the end-users and the material in this man page is aimed at package
developers.

`MsBackend` is a virtual class that defines what each different backend
needs to provide. `MsBackend` objects provide access to mass
spectrometry data. Such backends can be classified into *in-memory* or
*on-disk* backends, depending on where the data, i.e spectra (m/z and
intensities) and spectra annotation (MS level, charge, polarity, ...)
are stored.

Typically, in-memory backends keep all data in memory ensuring fast data
access, while on-disk backends store (parts of) their data on disk and
retrieve it on demand.

The *Backend functions and implementation notes for new backend classes*
section documents the API that a backend must implement.

Currently available backends are:

- `MsBackendMemory` and `MsBackendDataFrame`: store all data in memory.
  The `MsBackendMemory` is optimized for accessing and processing the
  peak data (i.e. the numerical matrices with the m/z and intensity
  values) while the `MsBackendDataFrame` keeps all data in a
  `DataFrame`.

- `MsBackendMzR`: stores the m/z and intensities on-disk in raw data
  files (typically `mzML` or `mzXML`) and the spectra annotation
  information (header) in memory in a `DataFrame`. This backend requires
  the `mzR` package.

- `MsBackendHdf5Peaks`: stores the m/z and intensities on-disk in custom
  hdf5 data files and the remaining spectra variables in memory (in a
  `DataFrame`). This backend requires the `rhdf5` package.

See below for more details about individual backends.

## Usage

``` r
# S4 method for class 'MsBackend'
backendBpparam(object, BPPARAM = bpparam())

# S4 method for class 'MsBackend'
backendInitialize(object, ...)

# S4 method for class 'list'
backendMerge(object, ...)

# S4 method for class 'MsBackend'
backendMerge(object, ...)

# S4 method for class 'MsBackend'
backendParallelFactor(object, ...)

# S4 method for class 'MsBackend'
export(object, ...)

# S4 method for class 'MsBackend'
acquisitionNum(object)

# S4 method for class 'MsBackend'
peaksData(object, columns = c("mz", "intensity"))

# S4 method for class 'MsBackend'
peaksVariables(object)

# S4 method for class 'MsBackend,dataframeOrDataFrameOrmatrix'
cbind2(x, y = data.frame(), ...)

# S4 method for class 'MsBackend'
centroided(object)

# S4 method for class 'MsBackend'
centroided(object) <- value

# S4 method for class 'MsBackend'
collisionEnergy(object)

# S4 method for class 'MsBackend'
collisionEnergy(object) <- value

# S4 method for class 'MsBackend'
dataOrigin(object)

# S4 method for class 'MsBackend'
dataOrigin(object) <- value

# S4 method for class 'MsBackend'
dataStorage(object)

# S4 method for class 'MsBackend'
dataStorage(object) <- value

# S4 method for class 'MsBackend'
dropNaSpectraVariables(object)

# S4 method for class 'MsBackend,ANY'
extractByIndex(object, i)

# S4 method for class 'MsBackend,missing'
extractByIndex(object, i)

# S4 method for class 'MsBackend'
filterAcquisitionNum(object, n, file, ...)

# S4 method for class 'MsBackend'
filterDataOrigin(object, dataOrigin = character())

# S4 method for class 'MsBackend'
filterDataStorage(object, dataStorage = character())

# S4 method for class 'MsBackend'
filterEmptySpectra(object, ...)

# S4 method for class 'MsBackend'
filterIsolationWindow(object, mz = numeric(), ...)

# S4 method for class 'MsBackend'
filterMsLevel(object, msLevel = integer())

# S4 method for class 'MsBackend'
filterPolarity(object, polarity = integer())

# S4 method for class 'MsBackend'
filterPrecursorMzRange(object, mz = numeric())

# S4 method for class 'MsBackend'
filterPrecursorMz(object, mz = numeric())

# S4 method for class 'MsBackend'
filterPrecursorMzValues(object, mz = numeric(), ppm = 20, tolerance = 0)

# S4 method for class 'MsBackend'
filterPrecursorCharge(object, z = integer())

# S4 method for class 'MsBackend'
filterPrecursorScan(object, acquisitionNum = integer(), f = dataOrigin(object))

# S4 method for class 'MsBackend'
filterRanges(
  object,
  spectraVariables = character(),
  ranges = numeric(),
  match = c("all", "any")
)

# S4 method for class 'MsBackend'
filterRt(object, rt = numeric(), msLevel. = integer())

# S4 method for class 'MsBackend'
filterValues(
  object,
  spectraVariables = character(),
  values = numeric(),
  ppm = 0,
  tolerance = 0,
  match = c("all", "any")
)

# S4 method for class 'MsBackend'
intensity(object)

# S4 method for class 'MsBackend'
intensity(object) <- value

# S4 method for class 'MsBackend'
ionCount(object)

# S4 method for class 'MsBackend'
isCentroided(object, ...)

# S4 method for class 'MsBackend'
isEmpty(x)

# S4 method for class 'MsBackend'
isolationWindowLowerMz(object)

# S4 method for class 'MsBackend'
isolationWindowLowerMz(object) <- value

# S4 method for class 'MsBackend'
isolationWindowTargetMz(object)

# S4 method for class 'MsBackend'
isolationWindowTargetMz(object) <- value

# S4 method for class 'MsBackend'
isolationWindowUpperMz(object)

# S4 method for class 'MsBackend'
isolationWindowUpperMz(object) <- value

# S4 method for class 'MsBackend'
isReadOnly(object)

# S4 method for class 'MsBackend'
length(x)

# S4 method for class 'MsBackend'
msLevel(object)

# S4 method for class 'MsBackend'
msLevel(object) <- value

# S4 method for class 'MsBackend'
mz(object)

# S4 method for class 'MsBackend'
mz(object) <- value

# S4 method for class 'MsBackend'
lengths(x, use.names = FALSE)

# S4 method for class 'MsBackend'
polarity(object)

# S4 method for class 'MsBackend'
polarity(object) <- value

# S4 method for class 'MsBackend'
precScanNum(object)

# S4 method for class 'MsBackend'
precursorCharge(object)

# S4 method for class 'MsBackend'
precursorIntensity(object)

# S4 method for class 'MsBackend'
precursorMz(object)

# S4 method for class 'MsBackend'
precursorMz(object, ...) <- value

# S4 method for class 'MsBackend'
peaksData(object) <- value

# S4 method for class 'MsBackend'
reset(object)

# S4 method for class 'MsBackend'
rtime(object)

# S4 method for class 'MsBackend'
rtime(object) <- value

# S4 method for class 'MsBackend'
scanIndex(object)

# S4 method for class 'MsBackend'
selectSpectraVariables(object, spectraVariables = spectraVariables(object))

# S4 method for class 'MsBackend'
smoothed(object)

# S4 method for class 'MsBackend'
smoothed(object) <- value

# S4 method for class 'MsBackend'
spectraData(object, columns = spectraVariables(object))

# S4 method for class 'MsBackend'
spectraData(object) <- value

# S4 method for class 'MsBackend'
spectraNames(object)

# S4 method for class 'MsBackend'
spectraNames(object) <- value

# S4 method for class 'MsBackend'
spectraVariables(object)

# S4 method for class 'MsBackend,ANY'
split(x, f, drop = FALSE, ...)

# S4 method for class 'MsBackend'
supportsSetBackend(object, ...)

# S4 method for class 'MsBackend'
tic(object, initial = TRUE)

# S4 method for class 'MsBackend'
x[i, j, ..., drop = FALSE]

# S4 method for class 'MsBackend'
x$name

# S4 method for class 'MsBackend'
x$name <- value

# S4 method for class 'MsBackend'
x[[i, j, ...]]

# S4 method for class 'MsBackend'
x[[i, j, ...]] <- value

# S4 method for class 'MsBackend'
uniqueMsLevels(object, ...)

# S4 method for class 'MsBackend'
dataStorageBasePath(object)

# S4 method for class 'MsBackend'
dataStorageBasePath(object) <- value

# S4 method for class 'MsBackend'
longForm(object, columns = spectraVariables(object))

MsBackendDataFrame()

# S4 method for class 'MsBackendDataFrame'
backendInitialize(object, data, peaksVariables = c("mz", "intensity"), ...)

MsBackendHdf5Peaks()

MsBackendMemory()

# S4 method for class 'MsBackendMemory'
backendInitialize(object, data, peaksVariables = c("mz", "intensity"), ...)

MsBackendMzR()
```

## Arguments

- object:

  Object extending `MsBackend`.

- BPPARAM:

  for `backendBpparam()`: parameter object from the `BiocParallel`
  package defining the parallel processing setup. Defaults to
  `BPPARAM = bpparam()`. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

- ...:

  Additional arguments.

- columns:

  For
  [`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  accessor: optional `character` with column names (spectra variables)
  that should be included in the returned `DataFrame`. By default, all
  columns are returned. For
  [`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  accessor: optional `character` with requested columns in the
  individual `matrix` of the returned `list`. Defaults to
  `peaksVariables(object)` and depends on what *peaks variables* the
  backend provides. For `longForm()`: the spectra and peaks variables
  that should be included in the returned `data.frame`. Defaults to
  `spectraVariables(object)` and is thus the union of spectra and peaks
  variables.

- x:

  Object extending `MsBackend`.

- y:

  For
  [`cbind2()`](https://rformassspectrometry.github.io/Spectra/reference/combineSpectra.md):
  A `data.frame` or `DataFrame` with the spectra variables to be added
  to the backend. The number of rows of `y` and their order have to
  match the number of spectra and their order in `x`.

- value:

  replacement value for `<-` methods. See individual method description
  or expected data type.

- i:

  For `[`: `integer`, `logical` or `character` to subset the object.

- n:

  for
  [`filterAcquisitionNum()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `integer` with the acquisition numbers to filter for.

- file:

  For `filterFile()`: index or name of the file(s) to which the data
  should be subsetted. For
  [`export()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md):
  `character` of length 1 or equal to the number of spectra.

- dataOrigin:

  For
  [`filterDataOrigin()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `character` to define which spectra to keep. For
  [`filterAcquisitionNum()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  optionally specify if filtering should occur only for spectra of
  selected `dataOrigin`.

- dataStorage:

  For
  [`filterDataStorage()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `character` to define which spectra to keep. For
  [`filterAcquisitionNum()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  optionally specify if filtering should occur only for spectra of
  selected `dataStorage`.

- mz:

  For
  [`filterIsolationWindow()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `numeric(1)` with the m/z value to filter the object. For
  [`filterPrecursorMzRange()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `numeric(2)` with the lower and upper m/z boundary. For
  [`filterPrecursorMzValues()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `numeric` with the m/z value(s) to filter the object.

- msLevel:

  `integer` defining the MS level of the spectra to which the function
  should be applied. For
  [`filterMsLevel()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  the MS level to which `object` should be subsetted.

- polarity:

  For
  [`filterPolarity()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `integer` specifying the polarity to to subset `object`.

- ppm:

  For
  [`filterPrecursorMzValues()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `numeric(1)` with the m/z-relative maximal acceptable difference for a
  m/z to be considered matching. See
  [`MsCoreUtils::closest()`](https://rdrr.io/pkg/MsCoreUtils/man/matching.html)
  for details. For
  [`filterValues()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `numeric` of any length allowing to define a maximal accepted
  difference between user input `values` and the `spectraVariables`
  values. If it is not equal to the length of the value provided with
  parameter `spectraVariables`, `ppm[1]` will be recycled.

- tolerance:

  For
  [`filterPrecursorMzValues()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `numeric(1)` with the maximal absolute acceptable difference for a m/z
  value to be considered matching. See
  [`MsCoreUtils::closest()`](https://rdrr.io/pkg/MsCoreUtils/man/matching.html)
  for details. For
  [`filterValues()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `numeric` accepted tolerance between the `values` and the spectra
  variables. Defaults to `tolerance = 0`. If it is not equal to the
  length of the value provided with parameter `spectraVariables`,
  `tolerance[1]` will be recycled.

- z:

  For
  [`filterPrecursorCharge()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  [`integer()`](https://rdrr.io/r/base/integer.html) with the precursor
  charges to be used as filter.

- acquisitionNum:

  for
  [`filterPrecursorScan()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `integer` with the acquisition number of the spectra to which the
  object should be subsetted.

- f:

  `factor` defining the grouping to split `x`. See
  [`split()`](https://rformassspectrometry.github.io/Spectra/reference/combineSpectra.md).
  For
  [`filterPrecursorScan()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  factor defining from which original data files the spectra derive to
  avoid selecting spectra from different samples/files. Defaults to
  `f = dataOrigin(object)`.

- spectraVariables:

  For
  [`selectSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `character` with the names of the spectra variables to which the
  backend should be subsetted. For
  [`filterRanges()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
  and
  [`filterValues()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `character` vector specifying the column(s) from `spectraData(object)`
  on which to filter the data and that correspond to the the names of
  the spectra variables that should be used for the filtering.

- ranges:

  for
  [`filterRanges()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  A `numeric` vector of paired values (upper and lower boundary) that
  define the ranges to filter the `object`. These paired values need to
  be in the same order as the `spectraVariables` parameter (see below).

- match:

  For
  [`filterRanges()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
  and
  [`filterValues()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `character(1) ` defining whether the condition has to match for all
  provided `ranges`/`values` (`match = "all"`; the default), or for any
  of them (`match = "any"`) for spectra to be retained.

- rt:

  for
  [`filterRt()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `numeric(2)` defining the retention time range to be used to
  subset/filter `object`.

- msLevel.:

  same as `msLevel` above.

- values:

  For
  [`filterValues()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  A `numeric` vector that define the values to filter the `object`.
  `values` needs to be of same length than parameter `spectraVariables`
  and in the same order.

- use.names:

  For
  [`lengths()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  whether spectrum names should be used.

- drop:

  For `[`: not considered.

- initial:

  For
  [`tic()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  `logical(1)` whether the initially reported total ion current should
  be reported, or whether the total ion current should be (re)calculated
  on the actual data (`initial = FALSE`).

- j:

  For `[`: not supported.

- name:

  For `$` and `$<-`: the name of the spectra variable to return or set.

- data:

  For `backendInitialize()`: `DataFrame` with spectrum metadata/data.
  This parameter can be empty for `MsBackendMzR` backends but needs to
  be provided for `MsBackendDataFrame` backends.

- peaksVariables:

  For `backendInitialize()` for `MsBackendMemory`: `character`
  specifying which of the columns of the provided `data` contain *peaks
  variables* (i.e. information for individual mass peaks). Defaults to
  `peaksVariables = c("mz", "intensity")`. `"mz"` and `"intensity"`
  should **always** be specified.

## Value

See documentation of respective function.

## Implementation notes

Backends extending `MsBackend` **must** implement all of its methods
(listed above). Developers of new `MsBackend`s should follow the
`MsBackendMemory` implementation. To ensure a new implementation being
conform with the `MsBackend` definition, developers should included test
suites provided by this package in their unit test setup. For that a
variable `be` should be created in the package's `"testthat.R"` file
that represents a (initialized) instance of the developed backend. Then
the path to the test suites should be defined with
`test_suite <- system.file("test_backends", "test_MsBackend", package = "Spectra")`
followed by `test_dir(test_suite)` to run all test files in that
directory. Individual unit test files could be run with
`test_file(file.path(test_suite, "test_spectra_variables.R"), stop_on_failure = TRUE)`
(note that without `stop_on_failure = TRUE` tests would fail silently) .
Adding this code to the packages `"testthat.R"` file ensures that all
tests checking the validity of an `MsBackend` instance defined in the
`Spectra` package are also run on the newly develped backend class.

The `MsBackend` defines the following slots:

- `@readonly`: `logical(1)` whether the backend supports
  writing/replacing of m/z or intensity values.

Backends extending `MsBackend` **must** implement all of its methods
(listed above). Developers of new `MsBackend`s should follow the
`MsBackendDataFrame` implementation.

The
[`MsBackendCached()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackendCached.md)
backend provides a caching mechanism to allow *read only* backends to
add or change spectra variables. This backend shouldn't be used on its
own, but is meant to be extended. See
[`MsBackendCached()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackendCached.md)
for details.

The `MsBackend` defines the following slots:

- `@readonly`: `logical(1)` whether the backend supports
  writing/replacing of m/z or intensity values.

## Backend functions

New backend classes **must** extend the base `MsBackend` class will have
to implement some of the following methods (see the `MsBackend` vignette
for detailed description and examples):

- `[`: subset the backend. Only subsetting by element (*row*/`i`) is
  allowed. Parameter `i` should support `integer` indices and `logical`
  and should throw an error if `i` is out of bounds. The
  [`MsCoreUtils::i2index`](https://rdrr.io/pkg/MsCoreUtils/man/i2index.html)
  could be used to check the input `i`. For `i = integer()` an empty
  backend should be returned. Implementation of this method is optional,
  as the default calls the `extractByIndex()` method (which has to be
  implemented as the main subsetting method).

- `$`, `$<-`: access or set/add a single spectrum variable (column) in
  the backend. Using a `value` of `NULL` should allow deleting the
  specified spectra variable. An error should be thrown if the spectra
  variable is not available.

- `[[`, `[[<-`: access or set/add a single spectrum variable (column) in
  the backend. The default implementation uses `$`, thus these methods
  don't have to be implemented for new classes extending `MsBackend`.

- [`acquisitionNum()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  returns the acquisition number of each spectrum. Returns an `integer`
  of length equal to the number of spectra (with `NA_integer_` if not
  available).

- `backendBpparam()`: return the parallel processing setup supported by
  the backend class. This function can be used by any higher level
  function to evaluate whether the provided parallel processing setup
  (or the default one returned by `bpparam()`) is supported by the
  backend. Backends not supporting parallel processing (e.g. because
  they contain a connection to a database that can not be shared across
  processes) should extend this method to return only `SerialParam()`
  and hence disable parallel processing for (most) methods and
  functions. See also `backendParallelFactor()` for a function to
  provide a preferred splitting of the backend for parallel processing.

- `backendInitialize()`: initialises the backend. This method is
  supposed to be called rights after creating an instance of the backend
  class and should prepare the backend (e.g. set the data for the memory
  backend or read the spectra header data for the `MsBackendMzR`
  backend). Parameters can be defined freely for each backend, depending
  on what is needed to initialize the backend. It is however suggested
  to also support a parameter `data` that can be used to submit the full
  spectra data as a `DataFrame` to the backend. This would allow the
  backend to be also usable for the
  [`setBackend()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  function from `Spectra`. Note that eventually (for *read-only*
  backends) also the `supportsSetBackend` method would need to be
  implemented to return `TRUE`. The `backendInitialize()` method has
  also to ensure to correctly set spectra variable `dataStorage`.

- `backendMerge()`: merges (combines) `MsBackend` objects into a single
  instance. All objects to be merged have to be of the same type (e.g.
  `MsBackendDataFrame()`).

- `backendParallelFactor()`: returns a `factor` defining an optimal
  (preferred) way how the backend can be split for parallel processing
  used for all peak data accessor or data manipulation functions. The
  default implementation returns a factor of length 0
  ([`factor()`](https://rdrr.io/r/base/factor.html)) providing thus no
  default splitting. `backendParallelFactor()` for `MsBackendMzR` on the
  other hand returns `factor(dataStorage(object))` hence suggesting to
  split the object by data file.

- `backendRequiredSpectraVariables()`: returns a `character` with
  spectra variable names that are mandatory for a specific backend. The
  default returns an empty
  [`character()`](https://rdrr.io/r/base/character.html). The
  implementation for `MsBackendMzR` returns
  `c("dataStorage", "scanIndex")` as these two spectra variables are
  required to load the MS data on-the-fly. This method needs only to be
  implemented if a backend requires specific variables to be defined.

- [`cbind2()`](https://rformassspectrometry.github.io/Spectra/reference/combineSpectra.md):
  allows to appends multiple new spectra variables to the backend at
  once. The values for the new spectra variables have to be in the same
  order as the spectra in `x`. Replacing existing spectra variables is
  not supported through this function. For a more controlled way of
  adding spectra variables, the
  [`joinSpectraData()`](https://rformassspectrometry.github.io/Spectra/reference/combineSpectra.md)
  should be used.

- [`centroided()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  `centroided<-`: gets or sets the centroiding information of the
  spectra.
  [`centroided()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  returns a `logical` vector of length equal to the number of spectra
  with `TRUE` if a spectrum is centroided, `FALSE` if it is in profile
  mode and `NA` if it is undefined. See also
  [`isCentroided()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  for estimating from the spectrum data whether the spectrum is
  centroided. `value` for `centroided<-` is either a single `logical` or
  a `logical` of length equal to the number of spectra in `object`.

- [`collisionEnergy()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  `collisionEnergy<-`: gets or sets the collision energy for all spectra
  in `object`.
  [`collisionEnergy()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  returns a `numeric` with length equal to the number of spectra
  (`NA_real_` if not present/defined), `collisionEnergy<-` takes a
  `numeric` of length equal to the number of spectra in `object`.

- [`dataOrigin()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  gets a `character` of length equal to the number of spectra in
  `object` with the *data origin* of each spectrum. This could e.g. be
  the mzML file from which the data was read.

- [`dataStorage()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  gets a `character` of length equal to the number of spectra in
  `object` with the data storage of each spectrum. Note that missing
  values (`NA_character_`) are not supported for `dataStorage`.

- `dataStorageBasePath()`,
  `dataStorageBasePath<-: gets or sets the common *base* path of the directory containing all data files. If supported, the function is expected to return (or accept) a `character`of length 1. Most backends (such as for example the`MsBackendMemory`will not support this function and`dataStorageBasePath()`will return`NA_character\_`. For `MsBackendMzR`, this function allows to get or change the path to the directory containing the original data files, which is required if e.g. a serialized `MsBackendMzR\`
  instance gets copied to another computer or file system.

- [`dropNaSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  removes spectra variables (i.e. columns in the object's `spectraData`
  that contain only missing values (`NA`). Note that while columns with
  only `NA`s are removed, a
  [`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  call after
  [`dropNaSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
  might still show columns containing `NA` values for *core* spectra
  variables.

- [`export()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md):
  exports data from a `Spectra` class to a file. This method is called
  by the `export,Spectra` method that passes itself as a second argument
  to the function. The `export,MsBackend` implementation is thus
  expected to take a `Spectra` class as second argument from which all
  data is exported. Taking data from a `Spectra` class ensures that also
  all eventual data manipulations (cached in the `Spectra`'s lazy
  evaluation queue) are applied prior to export - this would not be
  possible with only a MsBackend class. An example implementation is the
  [`export()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  method for the `MsBackendMzR` backend that supports export of the data
  in *mzML* or *mzXML* format. See the documentation for the
  `MsBackendMzR` class below for more information.

- `extractByIndex()`: function to subset a backend to selected elements
  defined by the provided index. Similar to `[`, this method should
  allow extracting (or to subset) the data in any order. In contrast to
  `[`, however, `i` is expected to be an `integer` (while `[` should
  also support `logical` and eventually `character`). While being
  apparently redundant to `[`, this methods avoids package namespace
  errors/problems that can result in implementations of `[` being not
  found by R (which can happen sometimes in parallel processing using
  the
  [`BiocParallel::SnowParam()`](https://rdrr.io/pkg/BiocParallel/man/SnowParam-class.html)).
  This method is used internally by `Spectra` to extract/subset its
  backend. Implementation of this method is mandatory.

- [`filterAcquisitionNum()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  filters the object keeping only spectra matching the provided
  acquisition numbers (argument `n`). If `dataOrigin` or `dataStorage`
  is also provided, `object` is subsetted to the spectra with an
  acquisition number equal to `n` **in spectra with matching dataOrigin
  or dataStorage values** retaining all other spectra.

- [`filterDataOrigin()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  filters the object retaining spectra matching the provided
  `dataOrigin`. Parameter `dataOrigin` has to be of type `character` and
  needs to match exactly the data origin value of the spectra to subset.
  [`filterDataOrigin()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
  should return the data ordered by the provided `dataOrigin` parameter,
  i.e. if `dataOrigin = c("2", "1")` was provided, the spectra in the
  resulting object should be ordered accordingly (first spectra from
  data origin `"2"` and then from `"1"`). Implementation of this method
  is optional since a default implementation for `MsBackend` is
  available.

- [`filterDataStorage()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  filters the object retaining spectra matching the provided
  `dataStorage`. Parameter `dataStorage` has to be of type `character`
  and needs to match exactly the data storage value of the spectra to
  subset.
  [`filterDataStorage()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
  should return the data ordered by the provided `dataStorage`
  parameter, i.e. if `dataStorage = c("2", "1")` was provided, the
  spectra in the resulting object should be ordered accordingly (first
  spectra from data storage `"2"` and then from `"1"`). Implementation
  of this method is optional since a default implementation for
  `MsBackend` is available.

- [`filterEmptySpectra()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  removes empty spectra (i.e. spectra without peaks). Implementation of
  this method is optional since a default implementation for `MsBackend`
  is available.

- `filterFile()`: retains data of files matching the file index or file
  name provided with parameter `file`.

- [`filterIsolationWindow()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  retains spectra that contain `mz` in their isolation window m/z range
  (i.e. with an `isolationWindowLowerMz` `<=` `mz` and
  `isolationWindowUpperMz` `>=` `mz`. Implementation of this method is
  optional since a default implementation for `MsBackend` is available.

- [`filterMsLevel()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  retains spectra of MS level `msLevel`. Implementation of this method
  is optional since a default implementation for `MsBackend` is
  available.

- [`filterPolarity()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  retains spectra of polarity `polarity`. Implementation of this method
  is optional since a default implementation for `MsBackend` is
  available.

- [`filterPrecursorMzRange()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
  (previously `filterPrecursorMz`): retains spectra with a precursor m/z
  within the provided m/z range. Implementation of this method is
  optional since a default implementation for `MsBackend` is available.

- [`filterPrecursorMzValues()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  retains spectra with a precursor m/z matching any of the provided m/z
  values (given `ppm` and `tolerance`). Implementation of this method is
  optional since a default implementation for `MsBackend` is available.

- [`filterPrecursorCharge()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  retains spectra with the defined precursor charge(s). Implementation
  of this method is optional since a default implementation for
  `MsBackend` is available.

- [`filterPrecursorScan()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  retains parent (e.g. MS1) and children scans (e.g. MS2) of acquisition
  number `acquisitionNum`. Parameter `f` is supposed to define the
  origin of the spectra (i.e. the original data file) to ensure related
  spectra from the same file/sample are selected and retained.
  Implementation of this method is optional since a default
  implementation for `MsBackend` is available.

- [`filterRanges()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  allows filtering of the `Spectra` object based on user defined
  *numeric* ranges (parameter `ranges`) for one or more available
  spectra variables in object (spectra variable names can be specified
  with parameter `spectraVariables`). Spectra for which the value of a
  spectra variable is within it's defined range are retained. If
  multiple ranges/spectra variables are defined, the `match` parameter
  can be used to specify whether all conditions (`match = "all"`; the
  default) or if any of the conditions must match (`match = "any"`; all
  spectra for which values are within any of the provided ranges are
  retained). Implementation of this method is optional since a default
  implementation for `MsBackend` is available.

- [`filterRt()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  retains spectra of MS level `msLevel` with retention times within
  (`>=`) `rt[1]` and (`<=`) `rt[2]`. The filter is applied to all
  spectra if no MS level is specified (the default,
  `msLevel. = integer()`). Implementation of this method is optional
  since a default implementation for `MsBackend` is available.

- [`filterValues()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  allows filtering of the `Spectra` object based on similarities of
  *numeric* values of one or more `spectraVariables(object)` (parameter
  `spectraVariables`) to provided values (parameter `values`) given
  acceptable differences (parameters tolerance and ppm). If multiple
  values/spectra variables are defined, the `match` parameter can be
  used to specify whether all conditions (`match = "all"`; the default)
  or if any of the conditions must match (`match = "any"`; all spectra
  for which values are within any of the provided ranges are retained).
  Implementation of this method is optional since a default
  implementation for `MsBackend` is available.

- [`intensity()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  gets the intensity values from the spectra. Returns a
  [`IRanges::NumericList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  of `numeric` vectors (intensity values for each spectrum). The length
  of the `list` is equal to the number of `spectra` in `object`.

- `intensity<-`: replaces the intensity values. `value` has to be a
  `list` (or
  [`IRanges::NumericList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html))
  of length equal to the number of spectra and the number of values
  within each list element identical to the number of peaks in each
  spectrum (i.e. the `lengths(x)`). Note that just writeable backends
  support this method.

- [`ionCount()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  returns a `numeric` with the sum of intensities for each spectrum. If
  the spectrum is empty (see
  [`isEmpty()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)),
  `NA_real_` is returned.

- [`isCentroided()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  a heuristic approach assessing if the spectra in `object` are in
  profile or centroided mode. The function takes the `qtl` th quantile
  top peaks, then calculates the difference between adjacent m/z value
  and returns `TRUE` if the first quartile is greater than `k`. (See
  `Spectra:::.peaks_is_centroided` for the code.)

- [`isEmpty()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  checks whether a spectrum in `object` is empty (i.e. does not contain
  any peaks). Returns a `logical` vector of length equal number of
  spectra.

- [`isolationWindowLowerMz()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  `isolationWindowLowerMz<-`: gets or sets the lower m/z boundary of the
  isolation window.

- [`isolationWindowTargetMz()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  `isolationWindowTargetMz<-`: gets or sets the target m/z of the
  isolation window.

- [`isolationWindowUpperMz()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  `isolationWindowUpperMz<-`: gets or sets the upper m/z boundary of the
  isolation window.

- `isReadOnly()`: returns a `logical(1)` whether the backend is *read
  only* or does allow also to write/update data.

- [`length()`](https://rdrr.io/r/base/length.html): returns the number
  of spectra in the object.

- [`lengths()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  gets the number of peaks (m/z-intensity values) per spectrum. Returns
  an `integer` vector (length equal to the number of spectra). For empty
  spectra, `0` is returned.

- `longForm()`: extract the MS data in *long form*, i.e., as a
  `data.frame` with columns being requested spectra and peak variables
  and one row per mass peak. Parameter `columns` can be used to specify
  the columns (i.e., spectra or peaks variables) that should be
  returned. The default is `columns = spectraVariables(object)` and
  **all** spectra and peak variables are returned. It is strongly
  suggested to extract only selected columns and not the full data to
  avoid potential out-of-memory problems. Implementation of this method
  is optional as a default implementation for `MsBackend` is available
  which converts the `DataFrame` returned by
  [`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  into long form.

- [`msLevel()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  gets the spectra's MS level. Returns an `integer` vector (of length
  equal to the number of spectra) with the MS level for each spectrum
  (or `NA_integer_` if not available).

- `msLevel<-`: replaces the spectra's MS level.

- [`mz()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  gets the mass-to-charge ratios (m/z) from the spectra. Returns a
  [`IRanges::NumericList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  or length equal to the number of spectra, each element a `numeric`
  vector with the m/z values of one spectrum.

- `mz<-`: replaces the m/z values. `value` has to be a `list` of length
  equal to the number of spectra and the number of values within each
  list element identical to the number of peaks in each spectrum (i.e.
  the `lengths(x)`). Note that just writeable backends support this
  method.

- [`polarity()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  `polarity<-`: gets or sets the polarity for each spectrum.
  [`polarity()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  returns an `integer` vector (length equal to the number of spectra),
  with `0` and `1` representing negative and positive polarities,
  respectively. `polarity<-` expects an integer vector of length 1 or
  equal to the number of spectra.

- [`precursorCharge()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  [`precursorIntensity()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  [`precursorMz()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  `precScanNum()`, `precAcquisitionNum()`: get the charge (`integer`),
  intensity (`numeric`), m/z (`numeric`), scan index (`integer`) and
  acquisition number (`interger`) of the precursor for MS level 2 and
  above spectra from the object. Returns a vector of length equal to the
  number of spectra in `object`. `NA` are reported for MS1 spectra of if
  no precursor information is available.

- [`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  returns a `list` with the spectras' peak data, i.e. m/z and intensity
  values or other *peak variables*. The length of the list is equal to
  the number of spectra in `object`. Each element of the list has to be
  a two-dimensional array (`matrix` or `data.frame`) with columns
  depending on the provided `columns` parameter (by default `"mz"` and
  `"intensity"`, but depends on the backend's available
  `peaksVariables`). For an empty spectrum, a `matrix` (`data.frame`)
  with 0 rows and columns according to `columns` is returned. The
  optional parameter `columns`, if supported by the backend, allows to
  define which peak variables should be returned in the `numeric` peak
  `matrix`. As a default `c("mz", "intensity")` should be used.

- `peaksData<-` replaces the peak data (m/z and intensity values) of the
  backend. This method expects a `list` of two dimensional arrays
  (`matrix` or `data.frame`) with columns representing the peak
  variables. All existing peaks data is expected to be replaced with
  these new values. The length of the `list` has to match the number of
  spectra of `object`. Note that only writeable backends need to support
  this method.

- [`peaksVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  lists the available variables for mass peaks. Default peak variables
  are `"mz"` and `"intensity"` (which all backends need to support and
  provide), but some backends might provide additional variables. All
  these variables are expected to be returned (if requested) by the
  [`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  function.

- [`reset()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
  a backend (if supported). This method will be called on the backend by
  the `reset,Spectra` method that is supposed to restore the data to its
  original state (see `reset,Spectra` for more details). The function
  returns the *reset* backend. The default implementation for
  `MsBackend` returns the backend as-is.

- [`rtime()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  `rtime<-`: gets or sets the retention times for each spectrum (in
  seconds).
  [`rtime()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  returns a `numeric` vector (length equal to the number of spectra)
  with the retention time for each spectrum. `rtime<-` expects a numeric
  vector with length equal to the number of spectra.

- [`scanIndex()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  returns an `integer` vector with the *scan index* for each spectrum.
  This represents the relative index of the spectrum within each file.
  Note that this can be different to the
  [`acquisitionNum()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  of the spectrum which is the index of the spectrum as reported in the
  mzML file.

- [`selectSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  reduces the information within the backend to the selected spectra
  variables. It is suggested to **not** remove values for the
  `"dataStorage"` variable, since this might be required for some
  backends to work properly (such as the `MsBackendMzR`).

- [`smoothed()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),`smoothed<-`:
  gets or sets whether a spectrum is *smoothed*.
  [`smoothed()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  returns a `logical` vector of length equal to the number of spectra.
  `smoothed<-` takes a `logical` vector of length 1 or equal to the
  number of spectra in `object`.

- [`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  `spectraData<-`: gets or sets general spectrum metadata (annotation,
  also called header).
  [`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  returns a `DataFrame`, `spectraData<-` expects a `DataFrame` with the
  same number of rows as there are spectra in `object`. Note that
  [`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  has to return the full data, i.e. also the m/z and intensity values
  (as a `list` or `SimpleList` in columns `"mz"` and `"intensity"`. See
  also
  [`fillCoreSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/fillCoreSpectraVariables.md)
  for a function that can *complete* a spectra data data frame with
  eventually missing *core* spectra variables.

- [`spectraNames()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  returns a `character` vector with the names of the spectra in `object`
  or `NULL` if not set. `spectraNames<-` allows to set spectra names (if
  the object is not read-only).

- [`spectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  returns a `character` vector with the available spectra variables
  (columns, fields or attributes) available in `object`. This should
  return **all** spectra variables which are present in `object`, also
  `"mz"` and `"intensity"` (which are by default not returned by the
  `spectraVariables,Spectra` method).

- [`split()`](https://rformassspectrometry.github.io/Spectra/reference/combineSpectra.md):
  splits the backend into a `list` of backends (depending on parameter
  `f`). The default method for `MsBackend` uses
  [`split.default()`](https://rdrr.io/r/base/split.html), thus backends
  extending `MsBackend` don't necessarily need to implement this method.

- `supportsSetBackend()`: whether a `MsBackend` supports the `Spectra`
  [`setBackend()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  function. For a `MsBackend` to support
  [`setBackend()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  it needs to have a parameter called `data` in its
  `backendInitialize()` method that support receiving all spectra data
  as a `DataFrame` from another backend and to initialize the backend
  with this data. In general *read-only* backends do not support
  [`setBackend()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  hence, the default implementation of `supportsSetBackend()` returns
  `!isReadOnly(object)`. If a read-only backend would support the
  [`setBackend()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  and being initialized with a `DataFrame` an implementation of this
  method for that backend could be defined that returns `TRUE` (see also
  the `MsBackend` vignette for details and examples).

- [`tic()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  gets the total ion current/count (sum of signal of a spectrum) for all
  spectra in `object`. By default, the value reported in the original
  raw data file is returned. For an empty spectrum, `NA_real_` is
  returned.

- [`uniqueMsLevels()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  gets the unique MS levels of all spectra in `object`. The default
  implementation calls `unique(msLevel(object))` but more efficient
  implementations could be defined for specific backends.

## Subsetting and merging backend classes

Backend classes must support (implement) the `[` method to subset the
object. This method should only support subsetting by spectra (rows,
`i`) and has to return a `MsBackend` class.

Backends extending `MsBackend` should also implement the
`backendMerge()` method to support combining backend instances (only
backend classes of the same type should be merged). Merging should
follow the following rules:

- The whole spectrum data of the various objects should be merged. The
  resulting merged object should contain the union of the individual
  objects' spectra variables (columns/fields), with eventually missing
  variables in one object being filled with `NA`.

## In-memory data backends

`MsBackendMemory` and `MsBackendDataFrame`:

The `MsBackendMemory` and `MsBackendDataFrame` objects keep all MS data
in memory are thus ideal for fast data processing. Due to their large
memory footprint they are however not suited for large scale
experiments. The two backends store the data different. The
`MsBackendDataFrame` stores all data in a `DataFrame` and thus supports
also S4-classes as spectra variables. Also, sepratate access to m/z or
intensity values (i.e. using the
[`mz()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
and
[`intensity()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
methods) is faster for the `MsBackendDataFrame`. The `MsBackendMemory`
on the other hand, due to the way the data is organized internally,
provides much faster access to the full peak data (i.e. the numerical
matrices of m/z and intensity values). Also subsetting and access to any
spectra variable (except `"mz"` and `"intensity"`) is fastest for the
`MsBackendMemory`.

Thus, for most use cases, the `MsBackendMemory` provides a higher
performance and flexibility than the `MsBackendDataFrame` and should
thus be preferred. See also issue
[246](https://github.com/rformassspectrometry/Spectra/issues/246) for a
performance comparison.

New objects can be created with the `MsBackendMemory()` and
`MsBackendDataFrame()` function, respectively. Both backends can be
subsequently initialized with the `backendInitialize()` method, taking a
`DataFrame` (or `data.frame`) with the (full) MS data as first parameter
`data`. The second parameter `peaksVariables` allows to define which
columns in `data` contain *peak variables* such as the m/z and intensity
values of individual peaks per spectrum. The default for this parameter
is `peaksVariables = c("mz", "intensity")`. Note that it is not
supported to provide either `"mz"` or `"intensity"`, if provided, both
need to be present in the data frame. Alternatively, the function also
supports a data frame without m/z and intensity values, in which case a
`Spectra` without mass peaks is created.

Suggested columns of this `DataFrame` are:

- `"msLevel"`: `integer` with MS levels of the spectra.

- `"rt"`: `numeric` with retention times of the spectra.

- `"acquisitionNum"`: `integer` with the acquisition number of the
  spectrum.

- `"scanIndex"`: `integer` with the index of the scan/spectrum within
  the *mzML*/*mzXML*/*CDF* file.

- `"dataOrigin"`: `character` defining the *data origin*.

- `"dataStorage"`: `character` indicating grouping of spectra in
  different e.g. input files. Note that missing values are not
  supported.

- `"centroided"`: `logical` whether the spectrum is centroided.

- `"smoothed"`: `logical` whether the spectrum was smoothed.

- `"polarity"`: `integer` with the polarity information of the spectra.

- `"precScanNum"`: `integer` specifying the index of the (MS1) spectrum
  containing the precursor of a (MS2) spectrum.

- `"precursorMz"`: `numeric` with the m/z value of the precursor.

- `"precursorIntensity"`: `numeric` with the intensity value of the
  precursor.

- `"precursorCharge"`: `integer` with the charge of the precursor.

- `"collisionEnergy"`: `numeric` with the collision energy.

- `"mz"`:
  [`IRanges::NumericList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  of `numeric` vectors representing the m/z values for each spectrum.

- `"intensity"`:
  [`IRanges::NumericList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  of `numeric` vectors representing the intensity values for each
  spectrum.

Additional columns are allowed too.

The
[`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
function for `MsBackendMemory` and `MsBackendDataFrame` returns a `list`
of `numeric` `matrix` by default (with parameter
`columns = c("mz", "intensity")`). If other peak variables are
requested, a `list` of `data.frame` is returned (ensuring m/z and
intensity values are always `numeric`).

## `MsBackendMzR`, on-disk MS data backend

The `MsBackendMzR` keeps only a limited amount of data in memory, while
the spectra data (m/z and intensity values) are fetched from the raw
files on-demand. This backend uses the `mzR` package for data import and
retrieval and hence requires that package to be installed. Also, it can
only be used to import and represent data stored in *mzML*, *mzXML* and
*CDF* files.

The `MsBackendMzR` backend extends the `MsBackendDataFrame` backend
using its `DataFrame` to keep spectra variables (except m/z and
intensity) in memory.

New objects can be created with the `MsBackendMzR()` function which can
be subsequently filled with data by calling `backendInitialize()`
passing the file names of the input data files with argument `files`.

This backend provides an
[`export()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
method to export data from a `Spectra` in *mzML* or *mzXML* format. The
definition of the function is:

`export(object, x, file = tempfile(), format = c("mzML", "mzXML"), copy = FALSE)`

The parameters are:

- `object`: an instance of the `MsBackendMzR` class.

- `x`: the
  [Spectra](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  object to be exported.

- `file`: `character` with the (full) output file name(s). Should be of
  length 1 or equal `length(x)`. If a single file is specified, all
  spectra are exported to that file. Alternatively it is possible to
  specify for each spectrum in `x` the name of the file to which it
  should be exported (and hence `file` has to be of length equal
  `length(x)`).

- `format`: `character(1)`, either `"mzML"` or `"mzXML"` defining the
  output file format.

- `copy`: `logical(1)` whether general file information should be copied
  from the original MS data files. This only works if `x` uses a
  `MsBackendMzR` backend and if `dataOrigin(x)` contains the original MS
  data file names.

- `BPPARAM`: parallel processing settings.

See examples in
[Spectra](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
or the vignette for more details and examples.

The `MsBackendMzR` ignores parameter `columns` of the
[`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
function and returns **always** m/z and intensity values.

## `MsBackendHdf5Peaks`, on-disk MS data backend

The `MsBackendHdf5Peaks` keeps, similar to the `MsBackendMzR`, peak data
(i.e. m/z and intensity values) in custom data files (in HDF5 format) on
disk while the remaining spectra variables are kept in memory. This
backend supports updating and writing of manipulated peak data to the
data files.

New objects can be created with the `MsBackendHdf5Peaks()` function
which can be subsequently filled with data by calling the object's
`backendInitialize()` method passing the desired file names of the HDF5
data files along with the spectra variables in form of a `DataFrame`
(see `MsBackendDataFrame` for the expected format). An optional
parameter `hdf5path` allows to specify the folder where the HDF5 data
files should be stored to. If provided, this is added as the path to the
submitted file names (parameter `files`).

By default `backendInitialize()` will store all peak data into a single
HDF5 file which name has to be provided with the parameter `files`. To
store peak data across several HDF5 files `data` has to contain a column
`"dataStorage"` that defines the grouping of spectra/peaks into files:
peaks for spectra with the same value in `"dataStorage"` are saved into
the same HDF5 file. If parameter `files` is omitted, the value in
`dataStorage` is used as file name (replacing any file ending with
`".h5"`. To specify the file names, `files`' length has to match the
number of unique elements in `"dataStorage"`.

For details see examples on the
[`Spectra()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
help page.

The `MsBackendHdf5Peaks` ignores parameter `columns` of the
[`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
function and returns **always** m/z and intensity values.

## Author

Johannes Rainer, Sebastian Gibb, Laurent Gatto, Philippine Louail

## Examples

``` r

## The MsBackend class is a virtual class and can not be instantiated
## directly. Below we define a new backend class extending this virtual
## class
MsBackendDummy <- setClass("MsBackendDummy", contains = "MsBackend")
MsBackendDummy()
#> An object of class "MsBackendDummy"
#> Slot "readonly":
#> [1] FALSE
#> 
#> Slot "version":
#> [1] "0.1"
#> 

## This class inherits now all methods from `MsBackend`, all of which
## however throw an error. These methods would have to be implemented
## for the new backend class.
try(mz(MsBackendDummy()))
#> Error in .local(object, ...) : 
#>   'mz()' is not implemented for MsBackendDummy.FALSE

## See `MsBackendDataFrame` as a reference implementation for a backend
## class (in the *R/MsBackendDataFrame.R* file).

## MsBackendDataFrame
##
## The `MsBackendDataFrame` uses a `S4Vectors::DataFrame` to store all MS
## data. Below we create such a backend by passing a `DataFrame` with all
## data to it.
data <- DataFrame(msLevel = c(1L, 2L, 1L), scanIndex = 1:3)
data$mz <- list(c(1.1, 1.2, 1.3), c(1.4, 54.2, 56.4, 122.1), c(15.3, 23.2))
data$intensity <- list(c(3, 2, 3), c(45, 100, 12.2, 1), c(123, 12324.2))

## Backends are supposed to be created with their specific constructor
## function
be <- MsBackendDataFrame()

be
#> MsBackendDataFrame with 0 spectra

## The `backendInitialize()` method initializes the backend filling it with
## data. This method can take any parameters needed for the backend to
## get loaded with the data (e.g. a file name from which to load the data,
## a database connection or, in this case, a data frame containing the data).
be <- backendInitialize(be, data)

be
#> MsBackendDataFrame with 3 spectra
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1        NA         1
#> 2         2        NA         2
#> 3         1        NA         3
#>  ... 16 more variables/columns.

## Data can be accessed with the accessor methods
msLevel(be)
#> [1] 1 2 1

mz(be)
#> NumericList of length 3
#> [[1]] 1.1 1.2 1.3
#> [[2]] 1.4 54.2 56.4 122.1
#> [[3]] 15.3 23.2

## Even if no data was provided for all spectra variables, its accessor
## methods are supposed to return a value.
precursorMz(be)
#> [1] NA NA NA

## The `peaksData()` method is supposed to return the peaks of the spectra as
## a `list`.
peaksData(be)
#> [[1]]
#>       mz intensity
#> [1,] 1.1         3
#> [2,] 1.2         2
#> [3,] 1.3         3
#> 
#> [[2]]
#>         mz intensity
#> [1,]   1.4      45.0
#> [2,]  54.2     100.0
#> [3,]  56.4      12.2
#> [4,] 122.1       1.0
#> 
#> [[3]]
#>        mz intensity
#> [1,] 15.3     123.0
#> [2,] 23.2   12324.2
#> 

## List available peaks variables
peaksVariables(be)
#> [1] "mz"        "intensity"

## Use columns to extract specific peaks variables. Below we extract m/z and
## intensity values, but in reversed order to the default.
peaksData(be, columns = c("intensity", "mz"))
#> [[1]]
#>      intensity  mz
#> [1,]         3 1.1
#> [2,]         2 1.2
#> [3,]         3 1.3
#> 
#> [[2]]
#>      intensity    mz
#> [1,]      45.0   1.4
#> [2,]     100.0  54.2
#> [3,]      12.2  56.4
#> [4,]       1.0 122.1
#> 
#> [[3]]
#>      intensity   mz
#> [1,]     123.0 15.3
#> [2,]   12324.2 23.2
#> 

## List available spectra variables (i.e. spectrum metadata)
spectraVariables(be)
#>  [1] "msLevel"                 "rtime"                  
#>  [3] "acquisitionNum"          "scanIndex"              
#>  [5] "mz"                      "intensity"              
#>  [7] "dataStorage"             "dataOrigin"             
#>  [9] "centroided"              "smoothed"               
#> [11] "polarity"                "precScanNum"            
#> [13] "precursorMz"             "precursorIntensity"     
#> [15] "precursorCharge"         "collisionEnergy"        
#> [17] "isolationWindowLowerMz"  "isolationWindowTargetMz"
#> [19] "isolationWindowUpperMz" 

## Extract precursor m/z, rtime, MS level spectra variables
spectraData(be, c("precursorMz", "rtime", "msLevel"))
#> DataFrame with 3 rows and 3 columns
#>   precursorMz     rtime   msLevel
#>     <numeric> <numeric> <integer>
#> 1          NA        NA         1
#> 2          NA        NA         2
#> 3          NA        NA         1

## MsBackendMemory
##
## The `MsBackendMemory` uses a more efficient internal data organization
## and allows also adding arbitrary additional peaks variables (annotations)
## Below we thus add a column "peak_ann" with arbitrary names/ids for each
## peak and add the name of this column to the `peaksVariables` parameter
## of the `backendInitialize()` method (in addition to `"mz"` and
## `"intensity"` that should **always** be specified.
data$peak_ann <- list(c("a", "", "d"), c("", "d", "e", "f"), c("h", "i"))
be <- backendInitialize(MsBackendMemory(), data,
    peaksVariables = c("mz", "intensity", "peak_ann"))
be
#> MsBackendMemory with 3 spectra
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1        NA         1
#> 2         2        NA         2
#> 3         1        NA         3
#>  ... 17 more variables/columns.

spectraVariables(be)
#>  [1] "msLevel"                 "rtime"                  
#>  [3] "acquisitionNum"          "scanIndex"              
#>  [5] "mz"                      "intensity"              
#>  [7] "dataStorage"             "dataOrigin"             
#>  [9] "centroided"              "smoothed"               
#> [11] "polarity"                "precScanNum"            
#> [13] "precursorMz"             "precursorIntensity"     
#> [15] "precursorCharge"         "collisionEnergy"        
#> [17] "isolationWindowLowerMz"  "isolationWindowTargetMz"
#> [19] "isolationWindowUpperMz"  "peak_ann"               

## peak_ann is also listed as a peaks variable
peaksVariables(be)
#> [1] "mz"        "intensity" "peak_ann" 

## The additional peaks variable can be accessed using the `peaksData()`
## function
peaksData(be, "peak_ann")
#> [[1]]
#>   peak_ann
#> 1        a
#> 2         
#> 3        d
#> 
#> [[2]]
#>   peak_ann
#> 1         
#> 2        d
#> 3        e
#> 4        f
#> 
#> [[3]]
#>   peak_ann
#> 1        h
#> 2        i
#> 

## The $<- method can be used to replace values of an existing peaks
## variable. It is important that the number of elements matches the
## number of peaks per spectrum.
be$peak_ann <- list(1:3, 1:4, 1:2)

## A peaks variable can again be removed by setting it to NULL
be$peak_ann <- NULL

peaksVariables(be)
#> [1] "mz"        "intensity"
```
