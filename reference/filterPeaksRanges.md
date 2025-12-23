# Filter peaks based on spectra and peaks variable ranges

The `filterPeaksRanges()` function allows to filter the peaks matrices
of a
[Spectra](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
object using any set of range-based filters on numeric spectra variables
or peaks variables. These ranges can be passed to the function using the
`...` as `<variable name> = <range>` pairs. `<variable name>` has to be
an available spectra or peaks variable. `<range>` can be a `numeric` of
length 2 defining the lower and upper boundary, or a `numeric`
two-column matrix (multi-row matrices are also supported, see further
below). `filterPeaksRanges(s, mz = c(200, 300))` would for example
reduce the peaks matrices of the `Spectra` object `s` to mass peaks with
an m/z value between 200 and 300. `filterPeaksRanges()` returns the
original `Spectra` object with the filter operation added to the
processing queue. Thus, the filter gets **only** applied when the peaks
data gets extracted with
[`mz()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
[`intensity()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
or
[`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md).
If ranges for both spectra **and** peaks variables are defined, the
function evaluates first whether the spectra variable value for a
spectrum is within the provided range and, if so, applies also the peaks
variable-based filter (otherwise an empty peaks matrix is returned).

If more than one spectra variable and/or peaks variable are defined,
their filter results are combined with a logical AND: a peak matrix is
only returned for a spectrum if all values of spectra variables are
within the provided (respective) ranges for spectra variables, and this
matrix is further filtered to contain only those peaks which values are
within the provided peaks variable ranges.

**Filtering with multiple ranges** per spectra and peaks variables is
also supported: ranges can also be provided as multi-row numeric
(two-column) matrices. In this case, the above described procedure is
applied for each row separately and their results are combined with a
logical OR, i.e. peaks matrices are returned that match any of the
conditions/filters of a row. The number of rows of the provided ranges
(being it for spectra or peaks variables) have to match.

**Missing value handling**: any comparison which involves a missing
value (being it a spectra variable value, a peaks variable value or a
value in one of the provided ranges) is treated as a logical `FALSE`.
For example, if the retention time of a spectrum is `NA` and the data is
filtered using a retention time range, an empty peaks matrix is returned
(for `keep = TRUE`, for `keep = FALSE` the full peaks matrix is
returned).

## Usage

``` r
filterPeaksRanges(object, ..., keep = TRUE)
```

## Arguments

- object:

  A
  [Spectra](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  object.

- ...:

  the ranges for the spectra and/or peaks variables. Has to be provided
  as `<name> = <range>` pairs with `<name>` being the name of a spectra
  or peaks variable (of numeric data type) and `<range>` being either a
  `numeric` of length 2 or a `numeric` two column matrix (see function
  desription above for details),

- keep:

  `logical(1)` whether to keep (default) or remove peaks that match the
  provided range(s).

## Note

In contrast to some other *filter* functions, this function does not
provide a `msLevel` parameter that allows to define the MS level of
spectra on which the filter should be applied. The filter(s) will always
be applied to **all** spectra (irrespectively of their MS level).
Through combination of multiple filter ranges it is however possible to
apply MS level-dependent filters (see examples below for details).

The filter will not be applied immediately to the data but only executed
when the mass peak data is accessed (through
[`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md),
[`mz()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
or
[`intensity()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md))
or by calling
[`applyProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md).

## Author

Johannes Rainer

## Examples

``` r

## Define a test Spectra
d <- data.frame(rtime = c(123.2, 134.2), msLevel = c(1L, 2L))
d$mz <- list(c(100.1, 100.2, 100.3, 200.1, 200.2, 300.3),
    c(100.3, 100.4, 200.2, 400.3, 400.4))
## Use the index of the mass peak within the spectrum as index for
## better illustration of filtering results
d$intensity <- list(c(1:6), 1:5)
s <- Spectra(d)
s
#> MSn data (Spectra) with 2 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1     123.2        NA
#> 2         2     134.2        NA
#>  ... 16 more variables/columns.

## Filter peaks removing all mass peaks with an m/z between 200 and 300
res <- filterPeaksRanges(s, mz = c(200, 300), keep = FALSE)
res
#> MSn data (Spectra) with 2 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1     123.2        NA
#> 2         2     134.2        NA
#>  ... 16 more variables/columns.
#> Lazy evaluation queue: 1 processing step(s)
#> Processing:
#>  Filter: remove peaks based on user-provided ranges for 1 variables [Tue Dec 23 13:21:25 2025] 

## The Spectra object has still the same length and spectra variables
length(res)
#> [1] 2
res$rtime
#> [1] 123.2 134.2

## The filter gets applied when mass peak data gets extracted, using either
## `mz()`, `intensity()` or `peaksData()`. The filtered peaks data does
## not contain any mass peaks with m/z values between 200 and 300:
peaksData(res)[[1L]]
#>         mz intensity
#> [1,] 100.1         1
#> [2,] 100.2         2
#> [3,] 100.3         3
#> [4,] 300.3         6
peaksData(res)[[2L]]
#>         mz intensity
#> [1,] 100.3         1
#> [2,] 100.4         2
#> [3,] 400.3         4
#> [4,] 400.4         5

## We next combine spectra and filter variables. We want to keep only mass
## peaks of MS2 spectra that have an m/z between 100 and 110.
res <- filterPeaksRanges(s, mz = c(100, 110), msLevel = c(2, 2))
res
#> MSn data (Spectra) with 2 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1     123.2        NA
#> 2         2     134.2        NA
#>  ... 16 more variables/columns.
#> Lazy evaluation queue: 1 processing step(s)
#> Processing:
#>  Filter: select peaks based on user-provided ranges for 2 variables [Tue Dec 23 13:21:26 2025] 
length(res)
#> [1] 2

## Only data for peaks are returned for which the spectra's MS level is
## between 2 and 2 and with an m/z between 100 and 110. The peaks data for
## the first spectrum, that has MS level 1, is thus empty:
peaksData(res)[[1L]]
#>      mz intensity

## While the peaks matrix for the second spectrum (with MS level 2) contains
## the mass peaks with m/z between 100 and 110.
peaksData(res)[[2L]]
#>         mz intensity
#> [1,] 100.3         1
#> [2,] 100.4         2

## To keep also the peaks data for the first spectrum, we need to define
## an additional set of ranges, which we define using a second row in each
## ranges matrix. We use the same filter as above, i.e. keeping only mass
## peaks with an m/z between 100 and 110 for spectra with MS level 2, but
## add an additional row for MS level 1 spectra keeping mass peaks with an
## m/z between 0 and 2000. Filter results of different rows are combined
## using a logical OR, i.e. peaks matrices with mass peaks are returned
## matching either the first, or the second row.
res <- filterPeaksRanges(s, mz = rbind(c(100, 110), c(0, 1000)),
    msLevel = rbind(c(2, 2), c(1, 1)))

## The results for the MS level 2 spectrum are the same as before, but with
## the additional row we keep the full peaks matrix of the MS1 spectrum:
peaksData(res)[[1L]]
#>         mz intensity
#> [1,] 100.1         1
#> [2,] 100.2         2
#> [3,] 100.3         3
#> [4,] 200.1         4
#> [5,] 200.2         5
#> [6,] 300.3         6
peaksData(res)[[2L]]
#>         mz intensity
#> [1,] 100.3         1
#> [2,] 100.4         2

## As a last example we define a filter that keeps all mass peaks with an
## m/z either between 100 and 200, or between 300 and 400.
res <- filterPeaksRanges(s, mz = rbind(c(100, 200), c(300, 400)))
peaksData(res)[[1L]]
#>         mz intensity
#> [1,] 100.1         1
#> [2,] 100.2         2
#> [3,] 100.3         3
#> [4,] 300.3         6
peaksData(res)[[2L]]
#>         mz intensity
#> [1,] 100.3         1
#> [2,] 100.4         2

## Such filters could thus be defined to restrict/filter the MS data to
## specific e.g. retention time and m/z ranges.
```
