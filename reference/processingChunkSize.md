# Parallel and chunk-wise processing of `Spectra`

Many operations on `Spectra` objects, specifically those working with
the actual MS data (peaks data), allow a chunk-wise processing in which
the `Spectra` is splitted into smaller parts (chunks) that are
iteratively processed. This enables parallel processing of the data (by
data chunk) and also reduces the memory demand since only the MS data of
the currently processed subset is loaded into memory and processed. This
chunk-wise processing, which is by default disabled, can be enabled by
setting the processing chunk size of a `Spectra` with the
`processingChunkSize()` function to a value which is smaller than the
length of the `Spectra` object. Setting
`processingChunkSize(sps) <- 1000` will cause any data manipulation
operation on the `sps`, such as
[`filterIntensity()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
or
[`bin()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md),
to be performed eventually in parallel for sets of 1000 spectra in each
iteration.

Such chunk-wise processing is specifically useful for `Spectra` objects
using an *on-disk* backend or for very large experiments. For small data
sets or `Spectra` using an in-memory backend, a direct processing might
however be more efficient. Setting the chunk size to `Inf` will disable
the chunk-wise processing.

For some backends a certain type of splitting and chunk-wise processing
might be preferable. The `MsBackendMzR` backend for example needs to
load the MS data from the original (mzML) files, hence chunk-wise
processing on a per-file basis would be ideal. The
[`backendParallelFactor()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
function for `MsBackend` allows backends to suggest a preferred
splitting of the data by returning a `factor` defining the respective
data chunks. The `MsBackendMzR` returns for example a `factor` based on
the *dataStorage* spectra variable. A `factor` of length 0 is returned
if no particular preferred splitting should be performed. The suggested
chunk definition will be used if no finite `processingChunkSize()` is
defined. Setting the `processingChunkSize` overrides
`backendParallelFactor`.

See the *Large-scale data handling and processing with Spectra* for more
information and examples.

Functions to configure parallel or chunk-wise processing:

- `processingChunkSize()`: allows to get or set the size of the chunks
  for parallel processing or chunk-wise processing of a `Spectra` in
  general. With a value of `Inf` (the default) no chunk-wise processing
  will be performed.

- `processingChunkFactor()`: returns a `factor` defining the chunks into
  which a `Spectra` will be split for chunk-wise (parallel) processing.
  A `factor` of length 0 indicates that no chunk-wise processing will be
  performed.

## Usage

``` r
# S4 method for class 'Spectra'
processingChunkSize(object)

# S4 method for class 'Spectra'
processingChunkSize(object) <- value

# S4 method for class 'Spectra'
processingChunkFactor(object)

# S4 method for class 'Spectra'
backendBpparam(object, BPPARAM = bpparam())
```

## Arguments

- object:

  `Spectra` object.

- value:

  `integer(1)` defining the chunk size.

- BPPARAM:

  Parallel setup configuration. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

## Value

`processingChunkSize()` returns the currently defined processing chunk
size (or `Inf` if it is not defined). `processingChunkFactor()` returns
a `factor` defining the chunks into which `x` will be split for
(parallel) chunk-wise processing or a `factor` of length 0 if no
splitting is defined.

## Note

Some backends might not support parallel processing at all. For these,
the
[`backendBpparam()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
function will always return a
[`SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
independently on how parallel processing was defined.

## Author

Johannes Rainer
