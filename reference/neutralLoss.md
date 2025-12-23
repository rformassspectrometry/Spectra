# Calculate Neutral Loss Spectra

This help page lists functions that convert MS/MS spectra to neutral
loss spectra. The main function for this is `neutralLoss` and the
specific algorithm to be used is defined (and configured) with dedicated
*parameter* objects (paramer `param` of the `neutralLoss()` function).

The parameter objects for the different algorithms are:

- `PrecursorMzParam()`: calculates neutral loss spectra as in Aisporna
  *et al.* 2022 by subtracting the (fragment's) peak m/z value from the
  precursor m/z value of each spectrum (precursor m/z - fragment m/z).
  Parameter `msLevel` allows to restrict calculation of neutral loss
  spectra to specified MS level(s). Spectra from other MS level(s) are
  returned as-is. Parameter `filterPeaks` allows to remove certain peaks
  from the neutral loss spectra. By default (`filterPeaks = "none"`) no
  filtering takes place. With `filterPeaks = "removePrecursor"` all
  fragment peaks with an m/z value matching the precursor m/z
  (considering also `ppm` and `tolerance` are removed. With
  `filterPeaks = "abovePrecursor"`, all fragment peaks with an m/z
  larger than the precursor m/z (m/z \> precursor m/z - `tolerance` -
  `ppm` of the precursor m/z) are removed (thus removing also in most
  cases the fragment peaks representing the precursor). Finally, with
  `filterPeaks = "belowPrecursor"` all fragment peaks with an m/z
  smaller than the precursor m/z (m/z \< precursor m/z + `tolerance` +
  `ppm` of the precursor m/z) are removed. Also in this case the
  precursor fragment peak is (depending on the values of `ppm` and
  `tolerance`) removed.

## Usage

``` r
PrecursorMzParam(
  filterPeaks = c("none", "abovePrecursor", "belowPrecursor", "removePrecursor"),
  msLevel = c(2L, NA_integer_),
  ppm = 10,
  tolerance = 0
)

# S4 method for class 'Spectra,PrecursorMzParam'
neutralLoss(object, param, ...)
```

## Arguments

- filterPeaks:

  For `PrecursorMzParam()`: `character(1)` or `function` defining if and
  how fragment peaks should be filtered before calculation. Pre-defined
  options are: `"none"` (keep all peaks), `"abovePrecursor"` (removes
  all fragment peaks with an m/z \>= precursor m/z), `"belowPrecursor"`
  (removes all fragment peaks with an m/z \<= precursor m/z). In
  addition, it is possible to pass a custom function with this parameter
  with arguments `x` (two column peak matrix) and `precursorMz` (the
  precursor m/z) that returns the sub-setted two column peak matrix.

- msLevel:

  `integer` defining for which MS level(s) the neutral loss spectra
  should be calculated. Defaults to `msLevel = c(2L, NA)` thus, neutral
  loss spectra will be calculated for all spectra with MS level equal to
  2 or with missing/undefined MS level. All spectra with a MS level
  different than `msLevel` will be returned unchanged.

- ppm:

  `numeric(1)` with m/z-relative acceptable difference in m/z values to
  filter peaks. Defaults to `ppm = 10`. See function description for
  details.

- tolerance:

  `numeric(1)` with absolute acceptable difference in m/z values to
  filter peaks. Defaults to `tolerance = 0`. See function description
  for details.

- object:

  [`Spectra()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  object with the fragment spectra for which neutral loss spectra should
  be calculated.

- param:

  One of the *parameter* objects discussed below.

- ...:

  Currently ignored.

## Value

A
[`Spectra()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
object with calculated neutral loss spectra.

## Note

By definition, mass peaks in a `Spectra` object need to be ordered by
their m/z value (in increasing order). Thus, the order of the peaks in
the calculated neutral loss spectra might not be the same than in the
original `Spectra` object.

Note also that for spectra with a missing precursor m/z empty spectra
are returned (i.e. spectra without peaks) since it is not possible to
calcualte the neutral loss spectra.

## References

Aisporna A, Benton PH, Chen A, Derks RJE, Galano JM, Giera M and Siuzdak
G (2022). Neutral Loss Mass Spectral Data Enhances Molecular Similarity
Analysis in METLIN. Journal of the American Society for Mass
Spectrometry.
[doi:10.1021/jasms.1c00343](https://doi.org/10.1021/jasms.1c00343)

## See also

[`addProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
for other data analysis and manipulation functions.

## Author

Johannes Rainer

## Examples

``` r

## Create a simple example Spectra object with some MS1, MS2 and MS3 spectra.
DF <- DataFrame(msLevel = c(1L, 2L, 3L, 1L, 2L, 3L),
                precursorMz = c(NA, 40, 20, NA, 300, 200))
DF$mz <- IRanges::NumericList(
                      c(3, 12, 14, 15, 16, 200),
                      c(13, 23, 39, 86),
                      c(5, 7, 20, 34, 50),
                      c(5, 7, 9, 20, 100),
                      c(15, 53, 299, 300),
                      c(34, 56, 100, 200, 204, 309)
                  , compress = FALSE)
DF$intensity <- IRanges::NumericList(1:6, 1:4, 1:5, 1:5, 1:4, 1:6,
                                     compress = FALSE)
sps <- Spectra(DF, backend = MsBackendDataFrame())

## Calculate neutral loss spectra for all MS2 spectra, keeping MS1 and MS3
## spectra unchanged.
sps_nl <- neutralLoss(sps, PrecursorMzParam(msLevel = 2L))
mz(sps)
#> NumericList of length 6
#> [[1]] 3 12 14 15 16 200
#> [[2]] 13 23 39 86
#> [[3]] 5 7 20 34 50
#> [[4]] 5 7 9 20 100
#> [[5]] 15 53 299 300
#> [[6]] 34 56 100 200 204 309
mz(sps_nl)
#> NumericList of length 6
#> [[1]] 3 12 14 15 16 200
#> [[2]] -46 1 17 27
#> [[3]] 5 7 20 34 50
#> [[4]] 5 7 9 20 100
#> [[5]] 0 1 247 285
#> [[6]] 34 56 100 200 204 309

## Calculate neutral loss spectra for MS2 and MS3 spectra, removing peaks
## with an m/z >= precursorMz
sps_nl <- neutralLoss(sps, PrecursorMzParam(
    filterPeaks = "abovePrecursor", msLevel = 2:3))
mz(sps_nl)
#> NumericList of length 6
#> [[1]] 3 12 14 15 16 200
#> [[2]] 1 17 27
#> [[3]] 13 15
#> [[4]] 5 7 9 20 100
#> [[5]] 1 247 285
#> [[6]] 100 144 166
## This removed also the peak with m/z 39 from the second spectrum

## Removing all fragment peaks matching the precursor m/z with a tolerance
## of 1 and ppm 10
sps_nl <- neutralLoss(sps, PrecursorMzParam(
    filterPeaks = "removePrecursor", tolerance = 1, ppm = 10, msLevel = 2:3))
mz(sps_nl)
#> NumericList of length 6
#> [[1]] 3 12 14 15 16 200
#> [[2]] -46 17 27
#> [[3]] -30 -14 13 15
#> [[4]] 5 7 9 20 100
#> [[5]] 247 285
#> [[6]] -109 -4 100 144 166

## Empty spectra are returned for MS 2 spectra with undefined precursor m/z.
sps$precursorMz <- NA_real_
sps_nl <- neutralLoss(sps, PrecursorMzParam())
mz(sps_nl)
#> NumericList of length 6
#> [[1]] 3 12 14 15 16 200
#> [[2]] numeric(0)
#> [[3]] 5 7 20 34 50
#> [[4]] 5 7 9 20 100
#> [[5]] numeric(0)
#> [[6]] 34 56 100 200 204 309
```
