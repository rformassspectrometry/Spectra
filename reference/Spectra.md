# The Spectra class to manage and access MS data

The `Spectra` class encapsules spectral mass spectrometry (MS) data and
related metadata. The MS data is represented by a *backend* extending
the virual
[MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
class which provides the data to the `Spectra` object. The `Spectra`
class implements only data accessor, filtering and analysis methods for
the MS data and relies on its *backend* to provide the MS data. This
allows to change data representations of a `Spectra` object depending on
the user's needs and properties of the data. Different backends and
their properties are explained in the
[MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
documentation.

Documentation on other topics and functionality of `Spectra`can be found
in:

- [`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  for accessing and using MS data through `Spectra` objects.

- [`filterMsLevel()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
  to subset and filter `Spectra` objects.

- [`plotSpectra()`](https://rformassspectrometry.github.io/Spectra/reference/spectra-plotting.md)
  for visualization of `Spectra` objects.

- [`processingChunkSize()`](https://rformassspectrometry.github.io/Spectra/reference/processingChunkSize.md)
  for information on parallel and chunk-wise data processing.

- [`combineSpectra()`](https://rformassspectrometry.github.io/Spectra/reference/combineSpectra.md)
  for merging, aggregating and splitting of `Spectra` objects.

- [`combinePeaks()`](https://rformassspectrometry.github.io/Spectra/reference/combinePeaks.md)
  for merging and aggregating `Spectra`'s mass peaks data.

- [`addProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
  for data analysis functions.

- [`compareSpectra()`](https://rformassspectrometry.github.io/Spectra/reference/compareSpectra.md)
  for spectra similarity calculations.

## Usage

``` r
# S4 method for class 'missing'
Spectra(
  object,
  processingQueue = list(),
  metadata = list(),
  ...,
  backend = MsBackendMemory(),
  BPPARAM = bpparam()
)

# S4 method for class 'MsBackend'
Spectra(object, processingQueue = list(), metadata = list(), ...)

# S4 method for class 'character'
Spectra(
  object,
  processingQueue = list(),
  metadata = list(),
  source = MsBackendMzR(),
  backend = source,
  ...,
  BPPARAM = bpparam()
)

# S4 method for class 'ANY'
Spectra(
  object,
  processingQueue = list(),
  metadata = list(),
  source = MsBackendMemory(),
  backend = source,
  ...,
  BPPARAM = bpparam()
)

# S4 method for class 'Spectra,MsBackend'
setBackend(
  object,
  backend,
  f = processingChunkFactor(object),
  ...,
  BPPARAM = bpparam()
)

# S4 method for class 'Spectra'
export(object, backend, ...)

# S4 method for class 'Spectra'
dataStorageBasePath(object)

# S4 method for class 'Spectra'
dataStorageBasePath(object) <- value
```

## Arguments

- object:

  For `Spectra()`: an object to instantiate the `Spectra` object and
  initialize the with data.. See section on creation of `Spectra`
  objects for details. For all other methods a `Spectra` object.

- processingQueue:

  For `Spectra()`: optional `list` of
  [ProtGenerics::ProcessingStep](https://rdrr.io/pkg/ProtGenerics/man/ProcessingStep.html)
  objects.

- metadata:

  For `Spectra()`: optional `list` with metadata information.

- ...:

  Additional arguments.

- backend:

  For `Spectra()`:
  [MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  to be used as backend. See section on creation of `Spectra` objects
  for details. For `setBackend()`: instance of
  [MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  that supports `setBackend()` (i.e. for which
  [`supportsSetBackend()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  returns `TRUE`). Such backends have a parameter `data` in their
  [`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  function that support passing the full spectra data to the initialize
  method. See section on creation of `Spectra` objects for details. For
  `export()`:
  [MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  to be used to export the data.

- BPPARAM:

  Parallel setup configuration. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information. This is passed directly to the
  [`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  method of the
  [MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md).

- source:

  For `Spectra()`: instance of
  [MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  that can be used to import spectrum data from the provided files. See
  section *Creation of objects* for more details.

- f:

  For `setBackend()`: factor defining how to split the data for
  parallelized copying of the spectra data to the new backend. For some
  backends changing this parameter can lead to errors. Defaults to
  [`processingChunkFactor()`](https://rformassspectrometry.github.io/Spectra/reference/processingChunkSize.md).

- value:

  For
  [`dataStorageBasePath()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md):
  A `character` vector that defines the base directory where the data
  storage files can be found.

## Details

The `Spectra` class uses by default a lazy data manipulation strategy,
i.e. data manipulations such as performed with
[`replaceIntensitiesBelow()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
are not applied immediately to the data, but applied on-the-fly to the
spectrum data once it is retrieved. This enables data manipulation
operations also for *read only* data representations. For some backends
that allow to write data back to the data storage (such as the
[`MsBackendMemory()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md),
[`MsBackendDataFrame()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
and
[`MsBackendHdf5Peaks()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md))
it is possible to apply to queue with the
[`applyProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
function (see the
[`applyProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
function for details).

Clarifications regarding scan/acquisition numbers and indices:

- A `spectrumId` (or `spectrumID`) is a vendor specific field in the
  mzML file that contains some information about the run/spectrum, e.g.:
  `controllerType=0 controllerNumber=1 scan=5281 file=2`

- `acquisitionNum` is a more a less sanitize spectrum id generated from
  the `spectrumId` field by `mzR` (see
  [here](https://github.com/sneumann/mzR/blob/master/src/pwiz/data/msdata/MSData.cpp#L552-L580)).

- `scanIndex` is the `mzR` generated sequence number of the spectrum in
  the raw file (which doesn't have to be the same as the
  `acquisitionNum`)

See also [this issue](https://github.com/lgatto/MSnbase/issues/525).

## Data stored in a `Spectra` object

The `Spectra` object is a container for MS data that includes mass peak
data (*m/z* and related intensity values, also referred to as *peaks
data* in the context of `Spectra`) and metadata of individual spectra
(so called *spectra variables*). While a core set of spectra variables
(the
[`coreSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md))
are guaranteed to be provided by a `Spectra`, it is possible to add
arbitrary additional spectra variables to a `Spectra` object.

The `Spectra` object is designed to contain MS data of a (large) set of
mass spectra. The data is organized *linearly* and can be thought of a
list of mass spectra, i.e. each element in the `Spectra` is one
spectrum.

## Creation of objects

`Spectra` classes can be created with the `Spectra()` constructor
function which supports the following formats:

- parameter `object` is a `data.frame` or `DataFrame` containing the
  full spectrum data (spectra variables in columns as well as columns
  with the individual MS peak data, *m/z* and intensity). The provided
  `backend` (by default a
  [MsBackendMemory](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md))
  will be initialized with that data.

- parameter `object` is a
  [MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  (assumed to be already initialized).

- parameter `object` is missing, in which case it is supposed that the
  data is provided by the
  [MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  class passed along with the `backend` argument.

- parameter `object` is of type `character` and is expected to be the
  file names(s) from which spectra should be imported. Parameter
  `source` allows to define a
  [MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  that is able to import the data from the provided source files. The
  default value for `source` is
  [`MsBackendMzR()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  which allows to import spectra data from mzML, mzXML or CDF files.

With `...` additional arguments can be passed to the backend's
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method. Parameter `backend` allows to specify which
[MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
should be used for data representation and storage.

## Data representation of a `Spectra`

The MS data which can be accessed through the `Spectra` object is
*represented* by its backend, which means that this backend defines how
and where the data is stored (e.g. in memory or on disk). The `Spectra`
object relies on the backend to provide the MS data whenever it needs it
for data processing. Different backends with different properties, such
as minimal memory requirement or fast data access, are defined in the
*Spectra* package or one of the MsBackend\* packages. More information
on backends and their properties is provided in the documentation of
[MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md).

On-disk backends keep only a limited amount of data in memory retrieving
most of the data (usually the MS peak data) upon request on-the-fly from
their on-disk data representations. Moving the on-disk data storage of
such a backend or a serialized object to a different location in the
file system will cause data corruption. The
[`dataStorageBasePath()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
and `dataStorageBasePath<-` functions allow in such cases (and if
thebackend classes support this operation), to get or change the *base*
path to the directory of the backend's data storage. In-memory backends
such as
[MsBackendMemory](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
or
[MsBackendDataFrame](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
keeping all MS data in memory don't support, and need, this function,
but for
[MsBackendMzR](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
this function can be used to update/adapt the path to the directory
containing the original data files. Thus, for `Spectra` objects (using
this backend) that were moved to another file system or computer, these
functions allow to adjust/adapt the base file path.

## Changing data representation of a `Spectra`

The data representation, i.e. the backend of a `Spectra` object can be
changed with the `setBackend()` method that takes an instance of the new
backend as second parameter `backend`. A call to
`setBackend(sps, backend = MsBackendDataFrame())` would for example
change the backend of `sps` to the *in-memory* `MsBackendDataFrame`.
Changing to a backend is only supported if that backend has a `data`
parameter in its
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method and if
[`supportsSetBackend()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
returns `TRUE` for that backend. `setBackend()` will transfer the full
spectra data from the originating backend as a `DataFrame` to the new
backend.

Generally, it is not possible to change **to** a read-only backend such
as the
[`MsBackendMzR()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
backend.

The definition of the function is:
`setBackend(object, backend, ..., f = dataStorage(object), BPPARAM = bpparam())`
and its parameters are:

- `object`: the `Spectra` object.

- `backend`: an instance of the new backend, e.g. `[MsBackendMemory()]`.

- `f`: factor allowing to parallelize the change of the backends. By
  default the process of copying the spectra data from the original to
  the new backend is performed separately (and in parallel) for each
  file. Users are advised to use the default setting.

- `...`: optional additional arguments passed to the
  [`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  method of the new `backend`.

- `BPPARAM`: setup for the parallel processing. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for details.

## Exporting data from a `Spectra` object

Data from a `Spectra` object can be **exported** to a file with the
`export()` function. The actual export of the data is performed by the
`export` method of the
[MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
class defined with the mandatory parameter `backend` which defines also
the format in which the data is exported. Note however that not all
backend classes support export of data. From the `MsBackend` classes in
the `Spectra` package currently only the `MsBackendMzR` backend supports
data export (to mzML/mzXML file(s)); see the help page of the
[MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
for information on its arguments or the examples below or the vignette
for examples.

The definition of the function is `export(object, backend, ...)` and its
parameters are:

- `object`: the `Spectra` object to be exported.

- `backend`: instance of a class extending
  [MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  which supports export of the data (i.e. which has a defined `export`
  method).

- `...`: additional parameters specific for the `MsBackend` passed with
  parameter `backend`.

## Author

Sebastian Gibb, Johannes Rainer, Laurent Gatto, Philippine Louail

## Examples

``` r

##  --------  CREATION OF SPECTRA OBJECTS  --------

## Create a Spectra providing a `DataFrame` containing the spectrum data.

spd <- DataFrame(msLevel = c(1L, 2L), rtime = c(1.1, 1.2))
spd$mz <- list(c(100, 103.2, 104.3, 106.5), c(45.6, 120.4, 190.2))
spd$intensity <- list(c(200, 400, 34.2, 17), c(12.3, 15.2, 6.8))

data <- Spectra(spd)
data
#> MSn data (Spectra) with 2 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1       1.1        NA
#> 2         2       1.2        NA
#>  ... 16 more variables/columns.

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


##  --------  CHANGING DATA REPRESENTATIONS  --------

## The MS data is on disk and will be read into memory on-demand. We can
## however change the backend to a MsBackendMemory backend which will
## keep all of the data in memory.
sciex_im <- setBackend(sciex, MsBackendMemory())
sciex_im
#> MSn data (Spectra) with 1862 spectra in a MsBackendMemory backend:
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
#> Processing:
#>  Switch backend from MsBackendMzR to MsBackendMemory [Tue Dec 23 13:20:10 2025] 

## The `MsBackendMemory()` supports the `setBackend()` method:
supportsSetBackend(MsBackendMemory())
#> [1] TRUE

## Thus, it is possible to change to that backend with `setBackend()`. Most
## read-only backends however don't support that, such as the
## `MsBackendMzR` and `setBackend()` would fail to change to that backend.
supportsSetBackend(MsBackendMzR())
#> [1] FALSE

## The on-disk object `sciex` is light-weight, because it does not keep the
## MS peak data in memory. The `sciex_im` object in contrast keeps all the
## data in memory and its size is thus much larger.
object.size(sciex)
#> 411016 bytes
object.size(sciex_im)
#> 55817368 bytes

## The spectra variable `dataStorage` returns for each spectrum the location
## where the data is stored. For in-memory objects:
head(dataStorage(sciex_im))
#> [1] "<memory>" "<memory>" "<memory>" "<memory>" "<memory>" "<memory>"

## While objects that use an on-disk backend will list the files where the
## data is stored.
head(dataStorage(sciex))
#> [1] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [2] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [3] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [4] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [5] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [6] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"

## The spectra variable `dataOrigin` returns for each spectrum the *origin*
## of the data. If the data is read from e.g. mzML files, this will be the
## original mzML file name:
head(dataOrigin(sciex))
#> [1] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [2] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [3] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [4] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [5] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [6] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
head(dataOrigin(sciex_im))
#> [1] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [2] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [3] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [4] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [5] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"
#> [6] "/__w/_temp/Library/msdata/sciex/20171016_POOL_POS_1_105-134.mzML"


##  -------- DATA EXPORT  --------

## Some `MsBackend` classes provide an `export()` method to export the data
## to the file format supported by the backend.
## The `MsBackendMzR` for example allows to export MS data to mzML or
## mzXML file(s), the `MsBackendMgf` (defined in the MsBackendMgf R package)
## would allow to export the data in mgf file format.
## Below we export the MS data in `data`. We call the `export()` method on
## this object, specify the backend that should be used to export the data
## (and which also defines the output format) and provide a file name.
fl <- tempfile()
export(data, MsBackendMzR(), file = fl)

## This exported our data in mzML format. Below we read the first 6 lines
## from that file.
readLines(fl, n = 6)
#> [1] "<?xml version=\"1.0\" encoding=\"utf-8\"?>"                                                                                                                                                                                            
#> [2] "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd\">"                 
#> [3] "  <mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\" id=\"\" version=\"1.1.0\">"
#> [4] "    <cvList count=\"2\">"                                                                                                                                                                                                              
#> [5] "      <cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"4.1.197\" URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"/>"                                      
#> [6] "      <cv id=\"UO\" fullName=\"Unit Ontology\" version=\"09:04:2014\" URI=\"https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo\"/>"                                                           

## If only a single file name is provided, all spectra are exported to that
## file. To export data with the `MsBackendMzR` backend to different files, a
## file name for each individual spectrum has to be provided.
## Below we export each spectrum to its own file.
fls <- c(tempfile(), tempfile())
export(data, MsBackendMzR(), file = fls)

## Reading the data from the first file
res <- Spectra(backendInitialize(MsBackendMzR(), fls[1]))

mz(res)
#> NumericList of length 1
#> [[1]] 100 103.2 104.3 106.5
mz(data)
#> NumericList of length 2
#> [[1]] 100 103.2 104.3 106.5
#> [[2]] 45.6 120.4 190.2
```
