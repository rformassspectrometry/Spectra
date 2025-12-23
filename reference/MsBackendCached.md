# Base MsBackend class providing data caching mechanism

The `MsBackendCached` class is a rudimentary implementation of the
[MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
providing a simple mechanism to cache spectra data locally. This class
is thought to be used as a base class for other `MsBackend`
implementations to reuse its caching mechanism and avoid having to
re-implement commonly used methods. This class is thus not thought to be
used directly by a user.

The `MsBackendCached` caching mechanism allows `MsBackend` instances to
add or replace spectra variables even if the backend used by them does
not allow to alter values (e.g. if a SQL database is used as a backend).
Any replacement operation with `$<-` will add the specified values to a
local `data.frame` within the `MsBackendCached` class that allows to
*cache* these values (increasing obviously the memory demand of the
object).

Any data accessor functions of the extending `MsBackend` class (such as
`$` or
[`msLevel()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
or
[`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md))
should first use `callNextMethod()` to call the respective accessor of
`MsBackendCached` that will evaluate if the requested spectra
variable(s) are in the local cache and return these. If the requested
spectra variables are neither in the local cache, nor listed in the
`@spectraVariables` slot (which defines all spectra variables that can
be requested from the extending `MsBackend` class) but are *core spectra
variables* then missing values of the correct data type are returned.

## Usage

``` r
MsBackendCached()

# S4 method for class 'MsBackendCached'
backendInitialize(
  object,
  data = data.frame(),
  nspectra = 0L,
  spectraVariables = character(),
  ...
)

# S4 method for class 'MsBackendCached'
dataStorage(object)

# S4 method for class 'MsBackendCached,ANY'
extractByIndex(object, i)

# S4 method for class 'MsBackendCached'
length(x)

# S4 method for class 'MsBackendCached'
spectraVariables(object)

# S4 method for class 'MsBackendCached'
spectraData(object, columns = spectraVariables(object))

# S4 method for class 'MsBackendCached'
spectraData(object) <- value

# S4 method for class 'MsBackendCached'
x[i, j, ..., drop = FALSE]

# S4 method for class 'MsBackendCached'
x$name

# S4 method for class 'MsBackendCached'
x$name <- value

# S4 method for class 'MsBackendCached'
selectSpectraVariables(object, spectraVariables = spectraVariables(object))

# S4 method for class 'MsBackendCached'
show(object)

# S4 method for class 'MsBackendCached'
intensity(object)

# S4 method for class 'MsBackendCached'
ionCount(object)

# S4 method for class 'MsBackendCached'
mz(object)
```

## Arguments

- object:

  A `MsBackendCached` object.

- data:

  For
  [`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md):
  (optional) `data.frame` with cached values. The number of rows (and
  their order) has to match the number of spectra.

- nspectra:

  For
  [`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md):
  `integer` with the number of spectra.

- spectraVariables:

  For
  [`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md):
  `character` with the names of the spectra variables that are provided
  by the extending backend. For
  [`selectSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md):
  `character` specifying the spectra variables to keep.

- ...:

  ignored

- i:

  For `[`: `integer` with the indices to subset the object.

- x:

  A `MsBackendCached` object.

- columns:

  For
  [`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  `character` with the names of the spectra variables to retrieve.

- value:

  replacement value for `<-` methods. See individual method description
  or expected data type.

- j:

  For `[`: ignored.

- drop:

  For `[`: not considered.

- name:

  For `$<-`: the name of the spectra variable to set.

## Value

See documentation of respective function.

## Implementation notes

Classes extending the `MsBackendCached` need to

- call the
  [`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  method of this class in their own
  [`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  method and set at least the number of spectra with the `nspectra`
  parameter and the `spectraVariables` that are available to the
  (extending) backend class.

- implement the
  [`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  method that also calls the
  [`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  method from `MsBackendCached` to also retrieve cached values (e.g.
  using `res <- callNextMethod()` at the beginning of the `spectraData`
  function). The `spectraData,MsBackendCached` method will return `NULL`
  if the selected spectra variables were not cached and are not *core
  spectra variables* not being provided by the extending backend. Thus,
  the extending backend can then proceed to retrieve the respective
  values from its own backend/data storage.

- implement eventually the `[` method that calls in addition the `[`
  from the `MsBackendCached`.

All other methods accessing or setting spectra variables don't need to
be implemented by the extending backend class (the default
implementations of the `MsBackendCached` will then be used instead;
these ensure that cached values are returned first). Spectra variables
can be modified or added using the `$<-` method of the
`MsBackendCached`. Replacing or adding multiple variables using the
`spectraData<-` is not supported by `MsBackendCached`. The extending
backend might however implement such a method that internally uses `$<-`
to add/replace single variables.

The `MsBackendCached` has the following slots:

- `nspectra`: `integer(1)` defining the number of spectra of the
  backend. This variable needs to be set and must match the number of
  rows of `localData` and the actual number of spectra in the
  (extending) backend.

- `localData`: `data.frame` with the cached local data. Any replacement
  operation with `$<-` will set/add a column with the respective values.

- `spectraVariables`: `character` defining the spectra variables that
  are provided by the extending `MsBackend` class (e.g. all spectra
  variables that can be retrieved from the data base or original data
  files).

## Available methods

- [`acquisitionNum()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  returns the acquisition number of each spectrum. Returns an `integer`
  of length equal to the number of spectra (with `NA_integer_` if not
  available).

- [`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md):
  *initializes* the backend. The method takes parameters `data`
  (`data.frame` with cached data), `nspectra` (`integer` defining the
  number of spectra) and `spectraVariables` (`character` with the
  spectra variables that are provided by the extending backend.

- [`centroided()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  `centroided<-`: gets or sets the centroiding information of the
  spectra. `centroided` returns a `logical` vector of length equal to
  the number of spectra with `TRUE` if a spectrum is centroided, `FALSE`
  if it is in profile mode and `NA` if it is undefined. See also
  `isCentroided` for estimating from the spectrum data whether the
  spectrum is centroided. `value` for `centroided<-` is either a single
  `logical` or a `logical` of length equal to the number of spectra in
  `object`.

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

- [`intensity()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  gets the intensity values from the spectra. Returns a
  [`IRanges::NumericList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  of `numeric` vectors (intensity values for each spectrum). The length
  of the `list` is equal to the number of `spectra` in `object`.

- [`ionCount()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  returns a `numeric` with the sum of intensities for each spectrum. If
  the spectrum is empty (see
  [`isEmpty()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)),
  `NA_real_` is returned.

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

- [`length()`](https://rdrr.io/r/base/length.html): returns the number
  of spectra (i.e. the `@nspectra`).

- [`lengths()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  gets the number of peaks (m/z-intensity values) per spectrum. Returns
  an `integer` vector (length equal to the number of spectra). For empty
  spectra, `0` is returned.

- [`msLevel()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  gets the spectra's MS level. Returns an `integer` vector (of length
  equal to the number of spectra) with the MS level for each spectrum
  (or `NA_integer_` if not available).

- [`mz()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  gets the mass-to-charge ratios (m/z) from the spectra. Returns a
  [`IRanges::NumericList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  or length equal to the number of spectra, each element a `numeric`
  vector with the m/z values of one spectrum.

- [`polarity()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
  `polarity<-`: gets or sets the polarity for each spectrum. `polarity`
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
  subset the object to specified spectra variables. This will eventually
  remove spectra variables listed in `@spectraVariables` and will also
  drop columns from the local cache if not among `spectraVariables`.

- [`smoothed()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),`smoothed<-`:
  gets or sets whether a spectrum is *smoothed*.
  [`smoothed()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  returns a `logical` vector of length equal to the number of spectra.
  `smoothed<-` takes a `logical` vector of length 1 or equal to the
  number of spectra in `object`.

- [`spectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  returns the available spectra variables, i.e. the unique set of *core
  spectra variables*, cached spectra variables and spectra variables
  defined in the `@spectraVariables` slot (i.e. spectra variables
  thought to be provided by the extending `MsBackend` instance).

- [`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md):
  returns a `DataFrame` with cached spectra variablers or initialized
  *core spectra variables*. Parameter `spectraVariables` allows to
  specify the variables to retrieve. The function returns `NULL` if the
  requested variables are not cached and are not provided by the
  extending backend. Note that this method **only** returns cached
  spectra variables or core spectra variables **not** provided by the
  extending backend. It is the responsibility of the extending backend
  to add/provide these.

- `[`: subsets the cached data. Parameter `i` needs to be an `integer`
  vector.

- `$`, `$<-`: access or set/add a single spectrum variable (column) in
  the backend.

## See also

[MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
for the documentation of MS backends.

## Author

Johannes Rainer
