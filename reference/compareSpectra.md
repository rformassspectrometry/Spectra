# Spectra similarity calculations

`compareSpectra()` compares each spectrum in `x` with each spectrum in
`y` using the function provided with `FUN` (defaults to
[`MsCoreUtils::ndotproduct()`](https://rdrr.io/pkg/MsCoreUtils/man/distance.html)).
If `y` is missing, each spectrum in `x` is compared with each other
spectrum in `x`. The matching/mapping of peaks between the compared
spectra is done with the `MAPFUN` function. The default
[`joinPeaks()`](https://rformassspectrometry.github.io/Spectra/reference/joinPeaks.md)
matches peaks of both spectra and allows to keep all peaks from the
first spectrum (`type = "left"`), from the second (`type = "right"`),
from both (`type = "outer"`) and to keep only matching peaks
(`type = "inner"`); see
[`joinPeaks()`](https://rformassspectrometry.github.io/Spectra/reference/joinPeaks.md)
for more information and examples). The `MAPFUN` function should have
parameters `x`, `y`, `xPrecursorMz` and `yPrecursorMz` as these values
are passed to the function.

In addition to
[`joinPeaks()`](https://rformassspectrometry.github.io/Spectra/reference/joinPeaks.md)
also
[`joinPeaksGnps()`](https://rformassspectrometry.github.io/Spectra/reference/joinPeaks.md)
is supported for GNPS-like similarity score calculations. Note that
[`joinPeaksGnps()`](https://rformassspectrometry.github.io/Spectra/reference/joinPeaks.md)
should only be used in combination with `FUN = MsCoreUtils::gnps` (see
[`joinPeaksGnps()`](https://rformassspectrometry.github.io/Spectra/reference/joinPeaks.md)
for more information and details). Use `MAPFUN = joinPeaksNone` to
disable internal peak matching/mapping if a similarity scoring function
is used that performs the matching internally.

`FUN` is supposed to be a function to compare intensities of (matched)
peaks of the two spectra that are compared. The function needs to take
two matrices with columns `"mz"` and `"intensity"` as input and is
supposed to return a single numeric as result. In addition to the two
peak matrices the spectra's precursor m/z values are passed to the
function as parameters `xPrecursorMz` (precursor m/z of the `x` peak
matrix) and `yPrecursorMz` (precursor m/z of the `y` peak matrix).
Additional parameters to functions `FUN` and `MAPFUN` can be passed with
`...`. Parameters `ppm` and `tolerance` are passed to both `MAPFUN` and
`FUN`. The function returns a `matrix` with the results of `FUN` for
each comparison, number of rows equal to `length(x)` and number of
columns equal `length(y)` (i.e. element in row 2 and column 3 is the
result from the comparison of `x[2]` with `y[3]`). If `SIMPLIFY = TRUE`
the `matrix` is *simplified* to a `numeric` if length of `x` or `y` is
one. See also the vignette for additional examples, such as using
spectral entropy similarity in the scoring.

## Usage

``` r
# S4 method for class 'Spectra,Spectra'
compareSpectra(
  x,
  y,
  MAPFUN = joinPeaks,
  tolerance = 0,
  ppm = 20,
  FUN = ndotproduct,
  ...,
  matchedPeaksCount = FALSE,
  SIMPLIFY = TRUE
)

# S4 method for class 'Spectra,missing'
compareSpectra(
  x,
  y = NULL,
  MAPFUN = joinPeaks,
  tolerance = 0,
  ppm = 20,
  FUN = ndotproduct,
  ...,
  matchedPeaksCount = FALSE,
  SIMPLIFY = TRUE
)
```

## Arguments

- x:

  A `Spectra` object.

- y:

  A `Spectra` object.

- MAPFUN:

  For `compareSpectra()`: function to map/match peaks between the two
  compared spectra. See
  [`joinPeaks()`](https://rformassspectrometry.github.io/Spectra/reference/joinPeaks.md)
  for more information and possible functions. Defaults to
  [`joinPeaks()`](https://rformassspectrometry.github.io/Spectra/reference/joinPeaks.md).

- tolerance:

  `numeric(1)` allowing to define a constant maximal accepted difference
  between m/z values for peaks to be matched. This parameter is directly
  passed to `MAPFUN`.

- ppm:

  `numeric(1)` defining a relative, m/z-dependent, maximal accepted
  difference between m/z values for peaks to be matched. This parameter
  is directly passed to `MAPFUN`.

- FUN:

  function to compare intensities of peaks between two spectra. Defaults
  to
  [`MsCoreUtils::ndotproduct()`](https://rdrr.io/pkg/MsCoreUtils/man/distance.html).

- ...:

  Additional arguments passed to the internal functions.

- matchedPeaksCount:

  `logical(1)` whether the number of matching peaks between the compared
  spectra should be returned in addition to the similarity scores. Note
  that with `matchedPeaksCount = TRUE` a 3-dimensional `array` is
  returned. See examples below for details.

- SIMPLIFY:

  `logical(1)` defining whether the result matrix should be *simplified*
  to a `numeric` if possible (i.e. if either `x` or `y` is of length 1).

## Value

For `matchedPeaksCount = FALSE` (the default) a `matrix` with the
similarity scores. With `matchedPeaksCount = FALSE` and
`SIMPLIFY = TRUE` a `numeric` vector. For `matchedPeaksCount = TRUE` a
3-dimensional array with the scores reported in the first matrix in z
dimension (`[, , 1]`) and the number of matching peaks in the second
matrix in z dimension (`[, , 2]`).

## Author

Sebastian Gibb, Johannes Rainer, Laurent Gatto

## Examples

``` r

## Load a `Spectra` object with LC-MS/MS data.
fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML",
    package = "msdata")
sps_dda <- Spectra(fl)
sps_dda
#> MSn data (Spectra) with 7602 spectra in a MsBackendMzR backend:
#>        msLevel     rtime scanIndex
#>      <integer> <numeric> <integer>
#> 1            1     0.231         1
#> 2            1     0.351         2
#> 3            1     0.471         3
#> 4            1     0.591         4
#> 5            1     0.711         5
#> ...        ...       ...       ...
#> 7598         1   899.491      7598
#> 7599         1   899.613      7599
#> 7600         1   899.747      7600
#> 7601         1   899.872      7601
#> 7602         1   899.993      7602
#>  ... 34 more variables/columns.
#> 
#> file(s):
#> PestMix1_DDA.mzML

## Restrict to MS2 (fragment) spectra:
sps_ms2 <- filterMsLevel(sps_dda, msLevel = 2L)

## Compare spectra: comparing spectra 2 and 3 against spectra 10:20 using
## the normalized dotproduct method.
res <- compareSpectra(sps_ms2[2:3], sps_ms2[10:20])
## first row contains comparisons of spectrum 2 with spectra 10 to 20 and
## the second row comparisons of spectrum 3 with spectra 10 to 20
res
#>      [,1] [,2]       [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#> [1,]    0    0 0.05181365    0    0    0    0    0    0     0     0
#> [2,]    0    0 0.20331585    0    0    0    0    0    0     0     0

## Setting parameter `matchedPeaksCount = TRUE` returns in addition to the
## simialrity score also the number of matching peaks between the compared
## spectra. The results are then returned as a 3-dimensional array, with the
## first matrix in z dimension (`[, , 1]`) containing the scores and the
## second matrix in z dimention (`[, , 2]`) the number of matching peaks:
res <- compareSpectra(sps_ms2[2:3], sps_ms2[10:20], matchedPeaksCount = TRUE)

## the scores
res[, , 1L]
#>      [,1] [,2]       [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#> [1,]    0    0 0.05181365    0    0    0    0    0    0     0     0
#> [2,]    0    0 0.20331585    0    0    0    0    0    0     0     0

## the number of matching peaks
res[, , 2L]
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#> [1,]    0    0    1    0    0    0    0    0    0     0     0
#> [2,]    0    0    1    0    0    0    0    0    0     0     0

## We next calculate the pairwise similarity for the first 10 spectra
compareSpectra(sps_ms2[1:10])
#>             [,1]       [,2]       [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#>  [1,] 1.00000000 0.03261724 0.12798950    0    0    0    0  NaN    0     0
#>  [2,] 0.03261724 1.00000000 0.07398183    0    0    0    0  NaN    0     0
#>  [3,] 0.12798950 0.07398183 1.00000000    0    0    0    0  NaN    0     0
#>  [4,] 0.00000000 0.00000000 0.00000000    1    0    0    0  NaN    0     0
#>  [5,] 0.00000000 0.00000000 0.00000000    0    1    0    0  NaN    0     0
#>  [6,] 0.00000000 0.00000000 0.00000000    0    0    1    0  NaN    0     0
#>  [7,] 0.00000000 0.00000000 0.00000000    0    0    0    1  NaN    0     0
#>  [8,]        NaN        NaN        NaN  NaN  NaN  NaN  NaN  NaN  NaN   NaN
#>  [9,] 0.00000000 0.00000000 0.00000000    0    0    0    0  NaN    1     0
#> [10,] 0.00000000 0.00000000 0.00000000    0    0    0    0  NaN    0     1

## Use compareSpectra to determine the number of common (matching) peaks
## with a ppm of 10:
## type = "inner" uses a *inner join* to match peaks, i.e. keeps only
## peaks that can be mapped betwen both spectra. The provided FUN returns
## simply the number of matching peaks.
compareSpectra(sps_ms2[2:3], sps_ms2[10:20], ppm = 10, type = "inner",
    FUN = function(x, y, ...) nrow(x))
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#> [1,]    0    0    0    0    0    0    0    0    0     0     0
#> [2,]    0    0    0    0    0    0    0    0    0     0     0

## We repeat this calculation between all pairwise combinations
## of the first 20 spectra
compareSpectra(sps_ms2[1:20], ppm = 10, type = "inner",
    FUN = function(x, y, ...) nrow(x))
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>  [1,]    3    0    0    0    0    0    0    0    0     0     0     1     0
#>  [2,]    0    3    1    0    0    0    0    0    0     0     0     0     0
#>  [3,]    0    1    3    0    0    0    0    0    0     0     0     0     0
#>  [4,]    0    0    0    3    0    0    0    0    0     0     0     0     0
#>  [5,]    0    0    0    0    1    0    0    0    0     0     0     0     0
#>  [6,]    0    0    0    0    0   11    0    0    0     0     0     0     0
#>  [7,]    0    0    0    0    0    0    8    0    0     0     0     0     0
#>  [8,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>  [9,]    0    0    0    0    0    0    0    0    7     0     0     0     0
#> [10,]    0    0    0    0    0    0    0    0    0     4     0     0     0
#> [11,]    0    0    0    0    0    0    0    0    0     0     3     0     0
#> [12,]    1    0    0    0    0    0    0    0    0     0     0     4     0
#> [13,]    0    0    0    0    0    0    0    0    0     0     0     0     2
#> [14,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [15,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [16,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [17,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [18,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [19,]    0    0    0    0    0    0    0    0    1     0     0     0     0
#> [20,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>       [,14] [,15] [,16] [,17] [,18] [,19] [,20]
#>  [1,]     0     0     0     0     0     0     0
#>  [2,]     0     0     0     0     0     0     0
#>  [3,]     0     0     0     0     0     0     0
#>  [4,]     0     0     0     0     0     0     0
#>  [5,]     0     0     0     0     0     0     0
#>  [6,]     0     0     0     0     0     0     0
#>  [7,]     0     0     0     0     0     0     0
#>  [8,]     0     0     0     0     0     0     0
#>  [9,]     0     0     0     0     0     1     0
#> [10,]     0     0     0     0     0     0     0
#> [11,]     0     0     0     0     0     0     0
#> [12,]     0     0     0     0     0     0     0
#> [13,]     0     0     0     0     0     0     0
#> [14,]    16     0     2     0     0     2     0
#> [15,]     0     1     0     0     0     0     0
#> [16,]     2     0     5     0     0     2     0
#> [17,]     0     0     0     3     1     0     0
#> [18,]     0     0     0     1     1     0     0
#> [19,]     2     0     2     0     0    11     0
#> [20,]     0     0     0     0     0     0     4
```
