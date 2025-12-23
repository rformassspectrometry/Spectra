# Join (map) peaks of two spectra

These functions map peaks from two spectra with each other if the
difference between their m/z values is smaller than defined with
parameters `tolerance` and `ppm`. All functions take two matrices

- `joinPeaks()`: maps peaks from two spectra allowing to specify the
  type of *join* that should be performed: `type = "outer"` each peak in
  `x` will be matched with each peak in `y`, for peaks that do not match
  any peak in the other spectra an `NA` intensity is returned. With
  `type = "left"` all peaks from the left spectrum (`x`) will be matched
  with peaks in `y`. Peaks in `y` that do not match any peak in `x` are
  omitted. `type = "right"` is the same as `type = "left"` only for `y`.
  Only peaks that can be matched between `x` and `y` are returned by
  `type = "inner"`, i.e. only peaks present in both spectra are
  reported.

- `joinPeaksGnps()`: matches/maps peaks between spectra with the same
  approach used in GNPS: peaks are considered matching if a) the
  difference in their m/z values is smaller than defined by `tolerance`
  and `ppm` (this is the same as `joinPeaks`) **and** b) the difference
  of their m/z *adjusted* for the difference of the spectras' precursor
  is smaller than defined by `tolerance` and `ppm`. Based on this
  definition, peaks in `x` can match up to two peaks in `y` hence peaks
  in the returned matrices might be reported multiple times. Note that
  if one of `xPrecursorMz` or `yPrecursorMz` are `NA` or if both are the
  same, the results are the same as with `joinPeaks()`. To calculate
  GNPS similarity scores,
  [`MsCoreUtils::gnps()`](https://rdrr.io/pkg/MsCoreUtils/man/gnps.html)
  should be called on the aligned peak matrices (i.e. `compareSpectra`
  should be called with `MAPFUN = joinPeaksGnps` and
  `FUN = MsCoreUtils::gnps`).

- `joinPeaksNone()`: does not perform any peak matching but simply
  returns the peak matrices in a `list`. This function should be used
  with the `MAPFUN` parameter of
  [`compareSpectra()`](https://rformassspectrometry.github.io/Spectra/reference/compareSpectra.md)
  if the spectra similarity function used (parameter `FUN` of
  [`compareSpectra()`](https://rformassspectrometry.github.io/Spectra/reference/compareSpectra.md))
  performs its own peak matching and does hence not expect matched peak
  matrices as an input.

## Usage

``` r
joinPeaks(x, y, type = "outer", tolerance = 0, ppm = 10, ...)

joinPeaksGnps(
  x,
  y,
  xPrecursorMz = NA_real_,
  yPrecursorMz = NA_real_,
  tolerance = 0,
  ppm = 0,
  type = "outer",
  ...
)

joinPeaksNone(x, y, ...)
```

## Arguments

- x:

  `matrix` with two columns `"mz"` and `"intensity"` containing the m/z
  and intensity values of the mass peaks of a spectrum.

- y:

  `matrix` with two columns `"mz"` and `"intensity"` containing the m/z
  and intensity values of the mass peaks of a spectrum.

- type:

  For `joinPeaks()` and `joinPeaksGnps()`: `character(1)` specifying the
  type of join that should be performed. See function description for
  details.

- tolerance:

  `numeric(1)` defining a constant maximal accepted difference between
  m/z values of peaks from the two spectra to be matched/mapped.

- ppm:

  `numeric(1)` defining a relative, m/z-dependent, maximal accepted
  difference between m/z values of peaks from the two spectra to be
  matched/mapped.

- ...:

  optional parameters passed to the
  [`MsCoreUtils::join()`](https://rdrr.io/pkg/MsCoreUtils/man/matching.html)
  function.

- xPrecursorMz:

  for `joinPeaksGnps()`: `numeric(1)` with the precursor m/z of the
  spectrum `x`.

- yPrecursorMz:

  for `joinPeaksGnps()`: `numeric(1)` with the precursor m/z of the
  spectrum `y`.

## Value

All functions return a `list` of elements `"x"` and `"y"` each being a
two column matrix with m/z (first column) and intensity values (second
column). The two matrices contain the matched peaks between input
matrices `x` and `y` and hence have the same number of rows. Peaks
present in `x` but not in the `y` input matrix have m/z and intensity
values of `NA` in the result matrix for `y` (and *vice versa*).

## Implementation notes

A mapping function must take two numeric matrices `x` and `y` as input
and must return `list` with two elements named `"x"` and `"y"` that
represent the aligned input matrices. The function should also have
`...` in its definition. Parameters `ppm` and `tolerance` are suggested
but not required.

## See also

- [`compareSpectra()`](https://rformassspectrometry.github.io/Spectra/reference/compareSpectra.md)
  for the function to calculate similarities between spectra.

- [`MsCoreUtils::gnps()`](https://rdrr.io/pkg/MsCoreUtils/man/gnps.html)
  in the *MsCoreUtils* package for more information on the GNPS
  similarity score.

## Author

Johannes Rainer, Michael Witting

## Examples

``` r

x <- cbind(c(31.34, 50.14, 60.3, 120.9, 230, 514.13, 874.1),
    1:7)
y <- cbind(c(12, 31.35, 70.3, 120.9 + ppm(120.9, 5),
    230 + ppm(230, 10), 315, 514.14, 901, 1202),
    1:9)

## No peaks with identical m/z
joinPeaks(x, y, ppm = 0, type = "inner")
#> $x
#>      [,1] [,2]
#> 
#> $y
#>      [,1] [,2]
#> 

## With ppm 10 two peaks are overlapping
joinPeaks(x, y, ppm = 10, type = "inner")
#> $x
#>       [,1] [,2]
#> [1,] 120.9    4
#> [2,] 230.0    5
#> 
#> $y
#>          [,1] [,2]
#> [1,] 120.9006    4
#> [2,] 230.0023    5
#> 

## Outer join: contain all peaks from x and y
joinPeaks(x, y, ppm = 10, type = "outer")
#> $x
#>         [,1] [,2]
#>  [1,]     NA   NA
#>  [2,]  31.34    1
#>  [3,]     NA   NA
#>  [4,]  50.14    2
#>  [5,]  60.30    3
#>  [6,]     NA   NA
#>  [7,] 120.90    4
#>  [8,] 230.00    5
#>  [9,]     NA   NA
#> [10,] 514.13    6
#> [11,]     NA   NA
#> [12,] 874.10    7
#> [13,]     NA   NA
#> [14,]     NA   NA
#> 
#> $y
#>            [,1] [,2]
#>  [1,]   12.0000    1
#>  [2,]        NA   NA
#>  [3,]   31.3500    2
#>  [4,]        NA   NA
#>  [5,]        NA   NA
#>  [6,]   70.3000    3
#>  [7,]  120.9006    4
#>  [8,]  230.0023    5
#>  [9,]  315.0000    6
#> [10,]        NA   NA
#> [11,]  514.1400    7
#> [12,]        NA   NA
#> [13,]  901.0000    8
#> [14,] 1202.0000    9
#> 

## Left join: keep all peaks from x and those from y that match
joinPeaks(x, y, ppm = 10, type = "left")
#> $x
#>        [,1] [,2]
#> [1,]  31.34    1
#> [2,]  50.14    2
#> [3,]  60.30    3
#> [4,] 120.90    4
#> [5,] 230.00    5
#> [6,] 514.13    6
#> [7,] 874.10    7
#> 
#> $y
#>          [,1] [,2]
#> [1,]       NA   NA
#> [2,]       NA   NA
#> [3,]       NA   NA
#> [4,] 120.9006    4
#> [5,] 230.0023    5
#> [6,]       NA   NA
#> [7,]       NA   NA
#> 

## Right join: keep all peaks from y and those from x that match. Using
## a constant tolerance of 0.01
joinPeaks(x, y, tolerance = 0.01, type = "right")
#> $x
#>         [,1] [,2]
#>  [1,]     NA   NA
#>  [2,]  31.34    1
#>  [3,]     NA   NA
#>  [4,] 120.90    4
#>  [5,] 230.00    5
#>  [6,]     NA   NA
#>  [7,] 514.13    6
#>  [8,]     NA   NA
#>  [9,]     NA   NA
#> 
#> $y
#>            [,1] [,2]
#>  [1,]   12.0000    1
#>  [2,]   31.3500    2
#>  [3,]   70.3000    3
#>  [4,]  120.9006    4
#>  [5,]  230.0023    5
#>  [6,]  315.0000    6
#>  [7,]  514.1400    7
#>  [8,]  901.0000    8
#>  [9,] 1202.0000    9
#> 

## GNPS-like peak matching

## Define spectra
x <- cbind(mz = c(10, 36, 63, 91, 93), intensity = c(14, 15, 999, 650, 1))
y <- cbind(mz = c(10, 12, 50, 63, 105), intensity = c(35, 5, 16, 999, 450))
## The precursor m/z
pmz_x <- 91
pmz_y <- 105

## Plain joinPeaks identifies only 2 matching peaks: 1 and 5
joinPeaks(x, y)
#> $x
#>      mz intensity
#> [1,] 10        14
#> [2,] NA        NA
#> [3,] 36        15
#> [4,] NA        NA
#> [5,] 63       999
#> [6,] 91       650
#> [7,] 93         1
#> [8,] NA        NA
#> 
#> $y
#>       mz intensity
#> [1,]  10        35
#> [2,]  12         5
#> [3,]  NA        NA
#> [4,]  50        16
#> [5,]  63       999
#> [6,]  NA        NA
#> [7,]  NA        NA
#> [8,] 105       450
#> 

## joinPeaksGnps finds 4 matches
joinPeaksGnps(x, y, pmz_x, pmz_y)
#> $x
#>       mz intensity
#>  [1,] 10        14
#>  [2,] 36        15
#>  [3,] 36        15
#>  [4,] 63       999
#>  [5,] 91       650
#>  [6,] 91       650
#>  [7,] 93         1
#>  [8,] NA        NA
#>  [9,] NA        NA
#> [10,] NA        NA
#> 
#> $y
#>        mz intensity
#>  [1,]  10        35
#>  [2,]  NA        NA
#>  [3,]  50        16
#>  [4,]  63       999
#>  [5,]  NA        NA
#>  [6,] 105       450
#>  [7,]  NA        NA
#>  [8,]  12         5
#>  [9,]  50        16
#> [10,] 105       450
#> 

## with one of the two precursor m/z being NA, the result are the same as
## with joinPeaks (with type = "left").
joinPeaksGnps(x, y, pmz_x, yPrecursorMz = NA)
#> $x
#>      mz intensity
#> [1,] 10        14
#> [2,] NA        NA
#> [3,] 36        15
#> [4,] NA        NA
#> [5,] 63       999
#> [6,] 91       650
#> [7,] 93         1
#> [8,] NA        NA
#> 
#> $y
#>       mz intensity
#> [1,]  10        35
#> [2,]  12         5
#> [3,]  NA        NA
#> [4,]  50        16
#> [5,]  63       999
#> [6,]  NA        NA
#> [7,]  NA        NA
#> [8,] 105       450
#> 
```
