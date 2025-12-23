# Accessing mass spectrometry data

As detailed in the documentation of the
[Spectra](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
class, a `Spectra` object is a container for mass spectrometry (MS) data
that includes both the mass peaks data (or *peaks data*, generally *m/z*
and intensity values) as well as spectra metadata (so called *spectra
variables*). Spectra variables generally define one value per spectrum,
while for peaks variables one value per mass peak is defined and hence
multiple values per spectrum (depending on the number of mass peaks of a
spectrum).

Data can be extracted from a `Spectra` object using dedicated accessor
functions or also using the `$` operator. Depending on the backend class
used by the `Spectra` to represent the data, data can also be added or
replaced (again, using dedicated functions or using `$<-`).

## Usage

``` r
asDataFrame(
  object,
  i = seq_along(object),
  spectraVars = spectraVariables(object)
)

# S4 method for class 'Spectra'
acquisitionNum(object)

# S4 method for class 'Spectra'
centroided(object)

# S4 method for class 'Spectra'
centroided(object) <- value

# S4 method for class 'Spectra'
collisionEnergy(object)

# S4 method for class 'Spectra'
collisionEnergy(object) <- value

coreSpectraVariables()

# S4 method for class 'Spectra'
dataOrigin(object)

# S4 method for class 'Spectra'
dataOrigin(object) <- value

# S4 method for class 'Spectra'
dataStorage(object)

# S4 method for class 'Spectra'
intensity(object, f = processingChunkFactor(object), ...)

# S4 method for class 'Spectra'
ionCount(object)

# S4 method for class 'Spectra'
isCentroided(object, ...)

# S4 method for class 'Spectra'
isEmpty(x)

# S4 method for class 'Spectra'
isolationWindowLowerMz(object)

# S4 method for class 'Spectra'
isolationWindowLowerMz(object) <- value

# S4 method for class 'Spectra'
isolationWindowTargetMz(object)

# S4 method for class 'Spectra'
isolationWindowTargetMz(object) <- value

# S4 method for class 'Spectra'
isolationWindowUpperMz(object)

# S4 method for class 'Spectra'
isolationWindowUpperMz(object) <- value

# S4 method for class 'Spectra'
length(x)

# S4 method for class 'Spectra'
lengths(x, use.names = FALSE)

# S4 method for class 'Spectra'
longForm(
  object,
  columns = union(spectraVariables(object), peaksVariables(object))
)

# S4 method for class 'Spectra'
msLevel(object)

# S4 method for class 'Spectra'
mz(object, f = processingChunkFactor(object), ...)

# S4 method for class 'Spectra'
peaksData(
  object,
  columns = c("mz", "intensity"),
  f = processingChunkFactor(object),
  return.type = c("SimpleList", "list"),
  ...,
  BPPARAM = bpparam()
)

# S4 method for class 'Spectra'
peaksVariables(object)

# S4 method for class 'Spectra'
polarity(object)

# S4 method for class 'Spectra'
polarity(object) <- value

# S4 method for class 'Spectra'
precScanNum(object)

# S4 method for class 'Spectra'
precursorCharge(object)

# S4 method for class 'Spectra'
precursorIntensity(object)

# S4 method for class 'Spectra'
precursorMz(object)

# S4 method for class 'Spectra'
precursorMz(object, ...) <- value

# S4 method for class 'Spectra'
rtime(object)

# S4 method for class 'Spectra'
rtime(object) <- value

# S4 method for class 'Spectra'
scanIndex(object)

# S4 method for class 'Spectra'
smoothed(object)

# S4 method for class 'Spectra'
smoothed(object) <- value

# S4 method for class 'Spectra'
spectraData(object, columns = spectraVariables(object))

# S4 method for class 'Spectra'
spectraData(object) <- value

# S4 method for class 'Spectra'
spectraNames(object)

# S4 method for class 'Spectra'
spectraNames(object) <- value

# S4 method for class 'Spectra'
spectraVariables(object)

# S4 method for class 'Spectra'
tic(object, initial = TRUE)

# S4 method for class 'Spectra'
uniqueMsLevels(object, ...)

# S4 method for class 'Spectra'
x$name

# S4 method for class 'Spectra'
x$name <- value

# S4 method for class 'Spectra'
x[[i, j, ...]]

# S4 method for class 'Spectra'
x[[i, j, ...]] <- value
```

## Arguments

- object:

  A `Spectra` object.

- i:

  For `asDataFrame()`: A `numeric` indicating which scans to coerce to a
  `DataFrame` (default is `seq_along(object)`).

- spectraVars:

  [`character()`](https://rdrr.io/r/base/character.html) indicating what
  spectra variables to add to the `DataFrame`. Default is
  `spectraVariables(object)`, i.e. all available variables.

- value:

  A vector with values to replace the respective spectra variable. Needs
  to be of the correct data type for the spectra variable.

- f:

  For `intensity()`, `mz()` and `peaksData()`: factor defining how data
  should be chunk-wise loaded an processed. Defaults to
  [`processingChunkFactor()`](https://rformassspectrometry.github.io/Spectra/reference/processingChunkSize.md).

- ...:

  Additional arguments.

- x:

  A `Spectra` object.

- use.names:

  For `lengths()`: ignored.

- columns:

  For `spectraData()` accessor: optional `character` with column names
  (spectra variables) that should be included in the returned
  `DataFrame`. By default, all columns are returned. For `peaksData()`
  accessor: optional `character` with requested columns in the
  individual `matrix` of the returned `list`. Defaults to
  `c("mz", "value")` but any values returned by `peaksVariables(object)`
  with `object` being the `Spectra` object are supported. For
  `longForm()`: `character` with the spectra and peaks variables to
  include in the returned `data.frame`. Defaults to
  `union(spectraVariables(object), peaksVariables(object))`.

- return.type:

  For `peaksData()`: `character(1)` allowing to specify if the results
  should be returned as a `SimpleList` or as a `list`. Defaults to
  `return.type = "SimpleList"`.

- BPPARAM:

  Parallel setup configuration. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information. See also
  [`processingChunkSize()`](https://rformassspectrometry.github.io/Spectra/reference/processingChunkSize.md)
  for more information on parallel processing.

- initial:

  For `tic()`: `logical(1)` whether the initially reported total ion
  current should be reported, or whether the total ion current should be
  (re)calculated on the actual data (`initial = FALSE`, same as
  `ionCount()`).

- name:

  For `$` and `$<-`: the name of the spectra variable to return or set.

- j:

  For `[`: not supported.

## Spectra variables

A common set of *core spectra variables* are defined for `Spectra`.
These have a pre-defined data type and each `Spectra` will return a
value for these if requested. If no value for a spectra variable is
defined, a missing value (of the correct data type) is returned. The
list of core spectra variables and their respective data type is:

- *acquisitionNum* `integer(1)`: the index of acquisition of a spectrum
  during an MS run.

- *centroided* `logical(1)`: whether the spectrum is in profile or
  centroid mode.

- *collisionEnergy* `numeric(1)`: collision energy used to create an MSn
  spectrum.

- *dataOrigin* `character(1)`: the *origin* of the spectrum's data, e.g.
  the mzML file from which it was read.

- *dataStorage* `character(1)`: the (current) storage location of the
  spectrum data. This value depends on the backend used to handle and
  provide the data. For an *in-memory* backend like the
  `MsBackendDataFrame` this will be `"<memory>"`, for an on-disk backend
  such as the `MsBackendHdf5Peaks` it will be the name of the HDF5 file
  where the spectrum's peak data is stored.

- *isolationWindowLowerMz* `numeric(1)`: lower m/z for the isolation
  window in which the (MSn) spectrum was measured.

- *isolationWindowTargetMz* `numeric(1)`: the target m/z for the
  isolation window in which the (MSn) spectrum was measured.

- *isolationWindowUpperMz* `numeric(1)`: upper m/z for the isolation
  window in which the (MSn) spectrum was measured.

- *msLevel* `integer(1)`: the MS level of the spectrum.

- *polarity* `integer(1)`: the polarity of the spectrum (`0` and `1`
  representing negative and positive polarity, respectively).

- *precScanNum* `integer(1)`: the scan (acquisition) number of the
  precursor for an MSn spectrum.

- *precursorCharge* `integer(1)`: the charge of the precursor of an MSn
  spectrum.

- *precursorIntensity* `numeric(1)`: the intensity of the precursor of
  an MSn spectrum.

- *precursorMz* `numeric(1)`: the m/z of the precursor of an MSn
  spectrum.

- *rtime* `numeric(1)`: the retention time of a spectrum.

- *scanIndex* `integer(1)`: the index of a spectrum within a (raw) file.

- *smoothed* `logical(1)`: whether the spectrum was smoothed.

For each of these spectra variable a dedicated accessor function is
defined (such as `msLevel()` or `rtime()`) that allows to extract the
values of that spectra variable for all spectra in a `Spectra` object.
Also, replacement functions are defined, but not all backends might
support replacing values for spectra variables. As described above,
additional spectra variables can be defined or added. The
`spectraVariables()` function can be used to

Values for multiple spectra variables, or all spectra vartiables\* can
be extracted with the `spectraData()` function.

## Peaks variables

`Spectra` also provide mass peak data with the *m/z* and intensity
values being the *core* peaks variables:

- *intensity* `numeric`: intensity values for the spectrum's peaks.

- *mz* `numeric`: the m/z values for the spectrum's peaks.

Values for these can be extracted with the `mz()` and `intensity()`
functions, or the `peaksData()` function. The former functions return a
`NumericList` with the respective values, while the latter returns a
`List` with `numeric` two-column matrices. The list of peaks matrices
can also be extracted using `as(x, "list")` or `as(x, "SimpleList")`
with `x` being a `Spectra` object.

Some `Spectra`/backends provide also values for additional peaks
variables. The set of available peaks variables can be extracted with
the `peaksVariables()` function.

## Functions to access MS data

The set of available functions to extract data from, or set data in, a
`Spectra` object are (in alphabetical order) listed below. Note that
there are also other functions to extract information from a `Spectra`
object documented in
[`addProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md).

- `$`, `$<-`: gets (or sets) a spectra variable for all spectra in
  `object`. See examples for details. Note that replacing values of a
  peaks variable is not supported with a non-empty processing queue,
  i.e. if any filtering or data manipulations on the peaks data was
  performed. In these cases
  [`applyProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
  needs to be called first to apply all cached data operations.

- `[[`, `[[<-`: access or set/add a single spectrum variable (column) in
  the backend.

- `acquisitionNum()`: returns the acquisition number of each spectrum.
  Returns an `integer` of length equal to the number of spectra (with
  `NA_integer_` if not available).

- `asDataFrame()`: converts the `Spectra` to a `DataFrame` (in long
  format) contining all data. Returns a `DataFrame`. See also
  `longForm()` for a potentially more efficient implementation that
  returns a `data.frame` in long form.

- `centroided()`, `centroided<-`: gets or sets the centroiding
  information of the spectra. `centroided()` returns a `logical` vector
  of length equal to the number of spectra with `TRUE` if a spectrum is
  centroided, `FALSE` if it is in profile mode and `NA` if it is
  undefined. See also `isCentroided()` for estimating from the spectrum
  data whether the spectrum is centroided. `value` for `centroided<-` is
  either a single `logical` or a `logical` of length equal to the number
  of spectra in `object`.

- `collisionEnergy()`, `collisionEnergy<-`: gets or sets the collision
  energy for all spectra in `object`. `collisionEnergy()` returns a
  `numeric` with length equal to the number of spectra (`NA_real_` if
  not present/defined), `collisionEnergy<-` takes a `numeric` of length
  equal to the number of spectra in `object`.

- `coreSpectraVariables()`: returns the *core* spectra variables along
  with their expected data type.

- `dataOrigin()`, `dataOrigin<-`: gets or sets the *data origin* for
  each spectrum. `dataOrigin()` returns a `character` vector (same
  length than `object`) with the origin of the spectra. `dataOrigin<-`
  expects a `character` vector (same length than `object`) with the
  replacement values for the data origin of each spectrum.

- `dataStorage()`: returns a `character` vector (same length than
  `object`) with the data storage location of each spectrum.

- `intensity()`: gets the intensity values from the spectra. Returns a
  [`IRanges::NumericList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  of `numeric` vectors (intensity values for each spectrum). The length
  of the list is equal to the number of `spectra` in `object`.

- `ionCount()`: returns a `numeric` with the sum of intensities for each
  spectrum. If the spectrum is empty (see `isEmpty()`), `NA_real_` is
  returned.

- `isCentroided()`: a heuristic approach assessing if the spectra in
  `object` are in profile or centroided mode. The function takes the
  `qtl`th quantile top peaks, then calculates the difference between
  adjacent m/z value and returns `TRUE` if the first quartile is greater
  than `k`. (See `Spectra:::.isCentroided()` for the code.)

- `isEmpty()`: checks whether a spectrum in `object` is empty (i.e. does
  not contain any peaks). Returns a `logical` vector of length equal
  number of spectra.

- `isolationWindowLowerMz()`, `isolationWindowLowerMz<-`: gets or sets
  the lower m/z boundary of the isolation window.

- `isolationWindowTargetMz()`, `isolationWindowTargetMz<-`: gets or sets
  the target m/z of the isolation window.

- `isolationWindowUpperMz()`, `isolationWindowUpperMz<-`: gets or sets
  the upper m/z boundary of the isolation window.

- [`length()`](https://rdrr.io/r/base/length.html): gets the number of
  spectra in the object.

- `lengths()`: gets the number of peaks (m/z-intensity values) per
  spectrum. Returns an `integer` vector (length equal to the number of
  spectra). For empty spectra, `0` is returned.

- `longForm()`: extract the MS data as a `data.frame` in *long form*
  with columns being spectra and peaks variables and one row per mass
  peak. Parameter `columns` allows to define the spectra and peaks
  variables that should be included in the returned `data.frame` (with
  the default being
  `columns = union(spectraVariables(object), peaksVariables(object)))`.

- `msLevel()`: gets the spectra's MS level. Returns an integer vector
  (names being spectrum names, length equal to the number of spectra)
  with the MS level for each spectrum.

- `mz()`: gets the mass-to-charge ratios (m/z) from the spectra. Returns
  a
  [`IRanges::NumericList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  or length equal to the number of spectra, each element a `numeric`
  vector with the m/z values of one spectrum.

- `peaksData()`: gets the *peaks* data for all spectra in `object`.
  Peaks data consist of the m/z and intensity values as well as possible
  additional annotations (variables) of all peaks of each spectrum. The
  function returns a
  [`S4Vectors::SimpleList()`](https://rdrr.io/pkg/S4Vectors/man/SimpleList-class.html)
  of two dimensional arrays (either `matrix` or `data.frame`), with each
  array providing the values for the requested *peak variables* (by
  default `"mz"` and `"intensity"`). Optional parameter `columns` is
  passed to the backend's `peaksData()` function to allow the selection
  of specific (or additional) peaks variables (columns) that should be
  extracted (if available). Importantly, it is **not** guaranteed that
  each backend supports this parameter (while each backend must support
  extraction of `"mz"` and `"intensity"` columns). Parameter `columns`
  defaults to `c("mz", "intensity")` but any value returned by
  `peaksVariables(object)` is supported. Note also that it is possible
  to extract the peak data with `as(x, "list")` and
  `as(x, "SimpleList")` as a `list` and `SimpleList`, respectively. Note
  however that, in contrast to `peaksData()`, `as()` does not support
  the parameter `columns`.

- `peaksVariables()`: lists the available variables for mass peaks
  provided by the backend. Default peak variables are `"mz"` and
  `"intensity"` (which all backends need to support and provide), but
  some backends might provide additional variables. These variables
  correspond to the column names of the peak data array returned by
  `peaksData()`.

- `polarity()`, `polarity<-`: gets or sets the polarity for each
  spectrum. `polarity()` returns an `integer` vector (length equal to
  the number of spectra), with `0` and `1` representing negative and
  positive polarities, respectively. `polarity<-` expects an `integer`
  vector of length 1 or equal to the number of spectra.

- `precursorCharge()`, `precursorIntensity()`, `precursorMz()`,
  `precScanNum()`, `precAcquisitionNum()`: gets the charge (`integer`),
  intensity (`numeric`), m/z (`numeric`), scan index (`integer`) and
  acquisition number (`interger`) of the precursor for MS level \> 2
  spectra from the object. Returns a vector of length equal to the
  number of spectra in `object`. `NA` are reported for MS1 spectra of if
  no precursor information is available.

- `rtime()`, `rtime<-`: gets or sets the retention times (in seconds)
  for each spectrum. `rtime()` returns a `numeric` vector (length equal
  to the number of spectra) with the retention time for each spectrum.
  `rtime<-` expects a numeric vector with length equal to the number of
  spectra.

- `scanIndex()`: returns an `integer` vector with the *scan index* for
  each spectrum. This represents the relative index of the spectrum
  within each file. Note that this can be different to the
  `acquisitionNum` of the spectrum which represents the index of the
  spectrum during acquisition/measurement (as reported in the mzML
  file).

- `smoothed()`,`smoothed<-`: gets or sets whether a spectrum is
  *smoothed*. `smoothed()` returns a `logical` vector of length equal to
  the number of spectra. `smoothed<-` takes a `logical` vector of length
  1 or equal to the number of spectra in `object`.

- `spectraData()`: gets general spectrum metadata (annotation, also
  called header). `spectraData()` returns a `DataFrame`. Note that this
  method does by default **not** return m/z or intensity values.

- `spectraData<-`: **replaces** the full spectra data of the `Spectra`
  object with the one provided with `value`. The `spectraData<-`
  function expects a `DataFrame` to be passed as value with the same
  number of rows as there a spectra in `object`. Note that replacing
  values of peaks variables is not supported with a non-empty processing
  queue, i.e. if any filtering or data manipulations on the peaks data
  was performed. In these cases
  [`applyProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
  needs to be called first to apply all cached data operations and empty
  the processing queue.

- `spectraNames()`, `spectraNames<-`: gets or sets the spectra names.

- `spectraVariables()`: returns a `character` vector with the available
  spectra variables (columns, fields or attributes of each spectrum)
  available in `object`. Note that `spectraVariables()` does not list
  the *peak variables* (`"mz"`, `"intensity"` and eventual additional
  annotations for each MS peak). Peak variables are returned by
  `peaksVariables()`.

- `tic()`: gets the total ion current/count (sum of signal of a
  spectrum) for all spectra in `object`. By default, the value reported
  in the original raw data file is returned. For an empty spectrum, `0`
  is returned.

- `uniqueMsLevels()`: get the unique MS levels available in `object`.
  This function is supposed to be more efficient than
  `unique(msLevel(object))`.

## See also

- [`addProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
  for functions to analyze `Spectra`.

- [Spectra](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  for a general description of the `Spectra` object.

## Author

Sebastian Gibb, Johannes Rainer, Laurent Gatto, Philippine Louail

## Examples

``` r

## Create a Spectra from mzML files and use the `MsBackendMzR` on-disk
## backend.
sciex_file <- dir(system.file("sciex", package = "msdata"),
    full.names = TRUE)
sciex <- Spectra(sciex_file, backend = MsBackendMzR())
sciex
#> MSn data (Spectra) with 1862 spectra in a MsBackendMzR backend:
#>        msLevel     rtime scanIndex
#>      <integer> <numeric> <integer>
#> 1            1     0.280         1
#> 2            1     0.559         2
#> 3            1     0.838         3
#> 4            1     1.117         4
#> 5            1     1.396         5
#> ...        ...       ...       ...
#> 1858         1   258.636       927
#> 1859         1   258.915       928
#> 1860         1   259.194       929
#> 1861         1   259.473       930
#> 1862         1   259.752       931
#>  ... 34 more variables/columns.
#> 
#> file(s):
#> 20171016_POOL_POS_1_105-134.mzML
#> 20171016_POOL_POS_3_105-134.mzML

## Get the number of spectra in the data set
length(sciex)
#> [1] 1862

## Get the number of mass peaks per spectrum - limit to the first 6
lengths(sciex) |> head()
#> [1]  578 1529 1600 1664 1417 1602

## Get the MS level for each spectrum - limit to the first 6 spectra
msLevel(sciex) |> head()
#> [1] 1 1 1 1 1 1

## Alternatively, we could also use $ to access a specific spectra variable.
## This could also be used to add additional spectra variables to the
## object (see further below).
sciex$msLevel |> head()
#> [1] 1 1 1 1 1 1

## Get the intensity and m/z values.
intensity(sciex)
#> NumericList of length 1862
#> [[1]] 0 412 0 0 412 0 0 412 0 0 412 0 0 ... 0 412 0 0 412 0 0 412 0 0 412 412 0
#> [[2]] 0 140 0 0 140 0 0 419 0 0 140 0 0 ... 0 140 0 0 140 0 0 140 0 0 279 140 0
#> [[3]] 0 132 263 263 132 132 0 0 132 132 0 0 ... 0 0 132 0 0 132 0 0 132 0 132 0
#> [[4]] 0 139 139 0 0 139 0 0 139 139 0 139 0 ... 0 0 139 0 0 277 0 0 139 0 139 0
#> [[5]] 0 164 0 0 328 0 164 0 0 164 0 0 164 ... 164 0 0 164 0 0 164 0 164 0 328 0
#> [[6]] 0 146 146 146 0 0 146 0 0 146 0 0 ... 146 0 0 146 146 0 0 146 0 0 146 0
#> [[7]] 0 296 0 296 0 0 148 0 0 148 0 0 148 ... 0 0 148 0 0 148 0 0 148 0 0 148 0
#> [[8]] 0 170 0 170 170 170 0 170 0 0 170 0 ... 170 0 0 170 0 0 170 0 0 170 170 0
#> [[9]] 0 157 0 314 0 0 157 0 0 157 0 0 314 ... 0 0 157 0 0 157 0 0 157 0 157 0
#> [[10]] 0 151 302 302 604 0 302 0 0 151 0 0 ... 151 0 0 151 0 151 0 151 0 151 0
#> ...
#> <1852 more elements>
mz(sciex)
#> NumericList of length 1862
#> [[1]] 105.043454833354 105.044900379521 ... 133.982027457992 133.983660012089
#> [[2]] 105.027517517902 105.028962955892 ... 133.982017159657 133.98364971537
#> [[3]] 105.037635723077 105.03908123069 ... 133.988547442185 133.990180037683
#> [[4]] 105.037635723077 105.03908123069 ... 133.98364971537 133.985282281029
#> [[5]] 105.034744757582 105.036190245303 ... 133.986914856634 133.988547442185
#> [[6]] 105.041972245917 105.043417783368 ... 133.982017159657 133.98364971537
#> [[7]] 105.037635723077 105.03908123069 ... 133.996710519135 133.998343164363
#> [[8]] 105.034744757582 105.036190245303 ... 133.980384613891 133.982017159657
#> [[9]] 105.040526728357 105.041972255863 ... 133.97875207807 133.980384613891
#> [[10]] 105.036190235357 105.037635733023 ... 133.985282281029 133.986914856634
#> ...
#> <1852 more elements>

## Convert a subset of the Spectra object to a long DataFrame.
asDataFrame(sciex, i = 1:3, spectraVars = c("rtime", "msLevel"))
#> DataFrame with 3707 rows and 4 columns
#>             mz intensity     rtime   msLevel
#>      <numeric> <numeric> <numeric> <integer>
#> 1      105.043         0      0.28         1
#> 2      105.045       412      0.28         1
#> 3      105.046         0      0.28         1
#> 4      107.055         0      0.28         1
#> 5      107.057       412      0.28         1
#> ...        ...       ...       ...       ...
#> 3703   133.984         0     0.838         1
#> 3704   133.985       132     0.838         1
#> 3705   133.987         0     0.838         1
#> 3706   133.989       132     0.838         1
#> 3707   133.990         0     0.838         1

## Create a Spectra providing a `DataFrame` containing the spectrum data.

spd <- DataFrame(msLevel = c(1L, 2L), rtime = c(1.1, 1.2))
spd$mz <- list(c(100, 103.2, 104.3, 106.5), c(45.6, 120.4, 190.2))
spd$intensity <- list(c(200, 400, 34.2, 17), c(12.3, 15.2, 6.8))

s <- Spectra(spd)
s
#> MSn data (Spectra) with 2 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1       1.1        NA
#> 2         2       1.2        NA
#>  ... 16 more variables/columns.

## List all available spectra variables (i.e. spectrum data and metadata).
spectraVariables(s)
#>  [1] "msLevel"                 "rtime"                  
#>  [3] "acquisitionNum"          "scanIndex"              
#>  [5] "dataStorage"             "dataOrigin"             
#>  [7] "centroided"              "smoothed"               
#>  [9] "polarity"                "precScanNum"            
#> [11] "precursorMz"             "precursorIntensity"     
#> [13] "precursorCharge"         "collisionEnergy"        
#> [15] "isolationWindowLowerMz"  "isolationWindowTargetMz"
#> [17] "isolationWindowUpperMz" 

## For all *core* spectrum variables accessor functions are available. These
## return NA if the variable was not set.
centroided(s)
#> [1] NA NA
dataStorage(s)
#> [1] "<memory>" "<memory>"
rtime(s)
#> [1] 1.1 1.2
precursorMz(s)
#> [1] NA NA

## The core spectra variables are:
coreSpectraVariables()
#>                 msLevel                   rtime          acquisitionNum 
#>               "integer"               "numeric"               "integer" 
#>               scanIndex                      mz               intensity 
#>               "integer"           "NumericList"           "NumericList" 
#>             dataStorage              dataOrigin              centroided 
#>             "character"             "character"               "logical" 
#>                smoothed                polarity             precScanNum 
#>               "logical"               "integer"               "integer" 
#>             precursorMz      precursorIntensity         precursorCharge 
#>               "numeric"               "numeric"               "integer" 
#>         collisionEnergy  isolationWindowLowerMz isolationWindowTargetMz 
#>               "numeric"               "numeric"               "numeric" 
#>  isolationWindowUpperMz 
#>               "numeric" 

## Add an additional metadata column.
s$spectrum_id <- c("sp_1", "sp_2")

## List spectra variables, "spectrum_id" is now also listed
spectraVariables(s)
#>  [1] "msLevel"                 "rtime"                  
#>  [3] "acquisitionNum"          "scanIndex"              
#>  [5] "dataStorage"             "dataOrigin"             
#>  [7] "centroided"              "smoothed"               
#>  [9] "polarity"                "precScanNum"            
#> [11] "precursorMz"             "precursorIntensity"     
#> [13] "precursorCharge"         "collisionEnergy"        
#> [15] "isolationWindowLowerMz"  "isolationWindowTargetMz"
#> [17] "isolationWindowUpperMz"  "spectrum_id"            

## Get the values for the new spectra variable
s$spectrum_id
#> [1] "sp_1" "sp_2"

## Extract specific spectra variables.
spectraData(s, columns = c("spectrum_id", "msLevel"))
#> DataFrame with 2 rows and 2 columns
#>   spectrum_id   msLevel
#>   <character> <integer>
#> 1        sp_1         1
#> 2        sp_2         2


##  --------  PEAKS VARIABLES AND DATA  --------

## Get the peak data (m/z and intensity values).
pks <- peaksData(s)
pks
#> List of length 2
pks[[1]]
#>         mz intensity
#> [1,] 100.0     200.0
#> [2,] 103.2     400.0
#> [3,] 104.3      34.2
#> [4,] 106.5      17.0
pks[[2]]
#>         mz intensity
#> [1,]  45.6      12.3
#> [2,] 120.4      15.2
#> [3,] 190.2       6.8

## Note that we could get the same resulb by coercing the `Spectra` to
## a `list` or `SimpleList`:
as(s, "list")
#> [[1]]
#>         mz intensity
#> [1,] 100.0     200.0
#> [2,] 103.2     400.0
#> [3,] 104.3      34.2
#> [4,] 106.5      17.0
#> 
#> [[2]]
#>         mz intensity
#> [1,]  45.6      12.3
#> [2,] 120.4      15.2
#> [3,] 190.2       6.8
#> 
as(s, "SimpleList")
#> List of length 2

## Or use `mz()` and `intensity()` to extract the m/z and intensity values
## separately
mz(s)
#> NumericList of length 2
#> [[1]] 100 103.2 104.3 106.5
#> [[2]] 45.6 120.4 190.2
intensity(s)
#> NumericList of length 2
#> [[1]] 200 400 34.2 17
#> [[2]] 12.3 15.2 6.8

## Some `MsBackend` classes provide support for arbitrary peaks variables
## (in addition to the mandatory `"mz"` and `"intensity"` values. Below
## we create a simple data frame with an additional peak variable `"pk_ann"`
## and create a `Spectra` with a `MsBackendMemory` for that data.
## Importantly the number of values (per spectrum) need to be the same
## for all peak variables.

tmp <- data.frame(msLevel = c(2L, 2L), rtime = c(123.2, 123.5))
tmp$mz <- list(c(103.1, 110.4, 303.1), c(343.2, 453.1))
tmp$intensity <- list(c(130.1, 543.1, 40), c(0.9, 0.45))
tmp$pk_ann <- list(c(NA_character_, "A", "P"), c("B", "P"))

## Create the Spectra. With parameter `peaksVariables` we can define
## the columns in `tmp` that contain peaks variables.
sps <- Spectra(tmp, source = MsBackendMemory(),
    peaksVariables = c("mz", "intensity", "pk_ann"))
peaksVariables(sps)
#> [1] "mz"        "intensity" "pk_ann"   

## Extract just the m/z and intensity values
peaksData(sps)[[1L]]
#>         mz intensity
#> [1,] 103.1     130.1
#> [2,] 110.4     543.1
#> [3,] 303.1      40.0

## Extract the full peaks data
peaksData(sps, columns = peaksVariables(sps))[[1L]]
#>      mz intensity pk_ann
#> 1 103.1     130.1   <NA>
#> 2 110.4     543.1      A
#> 3 303.1      40.0      P

## Access just the pk_ann variable
sps$pk_ann
#> [[1]]
#> [1] NA  "A" "P"
#> 
#> [[2]]
#> [1] "B" "P"
#> 

```
