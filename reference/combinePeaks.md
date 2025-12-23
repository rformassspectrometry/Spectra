# Aggregating and combining mass peaks data

In addition to aggregating content of spectra variables (describe in
[`combineSpectra()`](https://rformassspectrometry.github.io/Spectra/reference/combineSpectra.md))
it is also possible to aggregate and combine mass peaks data from
individual spectra within a `Spectra`. These `combinePeaks()` function
combines mass peaks **within each spectrum** with a difference in their
m/z values that is smaller than the maximal acceptable difference
defined by `ppm` and `tolerance`. Parameters `intensityFun` and `mzFun`
allow to define functions to aggregate the intensity and m/z values for
each such group of peaks. With `weighted = TRUE` (the default), the m/z
value of the combined peak is calculated using an intensity-weighted
mean and parameter `mzFun` is ignored. The
[`MsCoreUtils::group()`](https://rdrr.io/pkg/MsCoreUtils/man/group.html)
function is used for the grouping of mass peaks. Parameter `msLevel.`
allows to define selected MS levels for which peaks should be combined.
This function returns a `Spectra` with the same number of spectra than
the input object, but with possibly combined peaks within each spectrum.
Additional peak variables (other than `"mz"` and `"intensity"`) are
dropped (i.e. their values are replaced with `NA`) for combined peaks
unless they are constant across the combined peaks. See also
[`reduceSpectra()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
for a function to select a single *representative* mass peak for each
peak group.

## Usage

``` r
# S4 method for class 'Spectra'
combinePeaks(
  object,
  tolerance = 0,
  ppm = 20,
  intensityFun = base::mean,
  mzFun = base::mean,
  weighted = TRUE,
  msLevel. = uniqueMsLevels(object),
  ...
)
```

## Arguments

- object:

  A `Spectra` object.

- tolerance:

  `numeric(1)` allowing to define a constant maximal accepted difference
  between m/z values for peaks to be grouped. Default is
  `tolerance = 0`.

- ppm:

  `numeric(1)` defining a relative, m/z-dependent, maximal accepted
  difference between m/z values for peaks to be grouped. Default is
  `ppm = 20`.

- intensityFun:

  Function to aggregate intensities for all peaks in each peak group
  into a single intensity value.

- mzFun:

  Function to aggregate m/z values for all mass peaks within each peak
  group into a single m/z value. This parameter is ignored if
  `weighted = TRUE` (the default).

- weighted:

  `logical(1)` whether m/z values of peaks within each peak group should
  be aggregated into a single m/z value using an intensity-weighted
  mean. Defaults to `weighted = TRUE`.

- msLevel.:

  `integer` defining the MS level(s) of the spectra to which the
  function should be applied (defaults to all MS levels of `object`.

- ...:

  ignored.

## See also

- [`combineSpectra()`](https://rformassspectrometry.github.io/Spectra/reference/combineSpectra.md)
  for functions to combine or aggregate `Spectra`'s spectra data.

- [`combinePeaksData()`](https://rformassspectrometry.github.io/Spectra/reference/combinePeaksData.md)
  for the function to combine the mass peaks data.

- [`reduceSpectra()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
  and similar functions to filter mass peaks data.

- [Spectra](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  for a general description of the `Spectra` object.

## Author

Sebastian Gibb, Johannes Rainer, Laurent Gatto

## Examples

``` r

## Create a Spectra from mzML files and use the `MsBackendMzR` on-disk
## backend.
sciex_file <- dir(system.file("sciex", package = "msdata"),
    full.names = TRUE)
sciex <- Spectra(sciex_file, backend = MsBackendMzR())

## Combine mass peaks per spectrum with a difference in their m/z value
## that is smaller than 20 ppm. The intensity values of such peaks are
## combined by summing their values, while for the m/z values the median
## is reported
sciex_comb <- combinePeaks(sciex, ppm = 20,
    intensityFun = sum, mzFun = median)

## Comparing the number of mass peaks before and after aggregation
lengths(sciex) |> head()
#> [1]  578 1529 1600 1664 1417 1602
lengths(sciex_comb) |> head()
#> [1] 149 366 379 374 344 378

## Plotting the first spectrum before and after aggregation
par(mfrow = c(1, 2))
plotSpectra(sciex[2L])
plotSpectra(sciex_comb[2L])


## Using `reduceSpectra()` to keep for each group of mass peaks with a
## difference in their m/z values < 20ppm the one with the highest intensity.
sciex_red <- reduceSpectra(sciex, ppm = 20)

## Comparing the number of mass peaks before and after the operation
lengths(sciex) |> head()
#> [1]  578 1529 1600 1664 1417 1602
lengths(sciex_red) |> head()
#> [1] 149 366 379 374 344 378
```
