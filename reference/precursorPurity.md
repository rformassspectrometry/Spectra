# Calculating Precursor Purity for MS2 spectra

MS instruments generally collect precursor ions in a discrete *m/z*
*isolation window* before fragmenting them and recording the respective
fragment (MS2) spectrum. Ideally, only a single ion species is
fragmented, depending also on the size of the isolation window,
different ions (with slightly different *m/z*) might be fragmented. The
resulting MS2 spectrum might thus contain fragments from different ions
and hence be less *pure*.

The `precursorPurity()` function calculates the **precursor purity** of
MS2 (fragment) spectra expressed as the ratio between the itensity of
the highest signal in the isolation window to the sum of intensities of
all MS1 peaks in the isolation window. This is similar to the
calculation performed in the
[*msPurity*](https://www.bioconductor.org/packages/release/bioc/html/msPurity.html)
Bioconductor package.

The peak intensities within the isolation window is extracted from the
last MS1 spectrum before the respective MS2 spectrum. The spectra are
thus expected to be ordered by retention time. For the isolation window
either the isolation window reported in the `Spectra` object is used, or
it is calculated based on the MS2 spectra's precursor m/z. By default,
the isolation window is calculated based on the precursor m/z and
parameters `tolerance` and `ppm`: precursorMz +/- (`tolerance` +
`ppm(precursorMz, ppm)`). If the actually used precursor isolation
window is defined and available in the `Spectra` object, it can be used
instead by setting `useReportedIsolationWindow = TRUE` (default is
`useReportedIsolationWindow = FALSE`). Note that parameters `tolerance`
and `ppm` are ignored for `useReportedIsolationWindow = TRUE`.

## Usage

``` r
precursorPurity(
  object,
  tolerance = 0.05,
  ppm = 0,
  useReportedIsolationWindow = FALSE,
  BPPARAM = SerialParam()
)
```

## Arguments

- object:

  [`Spectra()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  object with LC-MS/MS data.

- tolerance:

  `numeric(1)` defining an absolute value (in Da) to be used to define
  the isolation window. For the precursor purity calculation of an MS2
  spectrum, all MS1 peaks from the previous MS1 scan with an m/z between
  the fragment spectrum's precursorMz +/- (tolerance + ppm(precursorMz,
  ppm)) are considered.

- ppm:

  `numeric(1)` defining the m/z dependent acceptable difference in m/z.
  See documentation of parameter `tolerance` for more information.

- useReportedIsolationWindow:

  `logical(1)` whether the reported isolation window, defined by spectra
  variables `isolationWindowLowerMz` and `isolationWindowUpperMz` in the
  input
  [Spectra](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  object, should be used instead of calculating the isolation window
  from the reported precursor m/z and parameters `tolerance` and `ppm`.
  Only few manufacturers report the isolation window with the spectra
  variables `isolationWindowLowerMz` and `isolationWindowTargetMz`, thus
  the default for this parameter is `FALSE`.

- BPPARAM:

  parallel processing setup. Defaults to `BPPARAM = SerialParam()`. See
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  for more information.

## Value

`numeric` vector of length equal to the number of spectra in `object`,
with values representing the calculated precursor purity for each
spectrum. For MS1 spectra, `NA_real_` is returned. For MS2 spectra, the
purity is defined as the proportion of maximum signal to the total ion
current within the isolation window that is attributable to the selected
precursor ion. If no matching MS1 scan is found or the precursor peak is
missing, `NA_real_` is returned.

## Note

This approach is applicable only when fragment spectra are obtained
through data-dependent acquisition (DDA), as it assumes that the peak
with the highest intensity within the given isolation m/z window (from
the previous MS1 spectrum) corresponds to the precursor ion.

The spectra in `object` have to be ordered by their retention time.

## See also

[`addProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
for other data analysis and manipulation functions.

## Author

Ahlam Mentag, Johannes Rainer

## Examples

``` r

## Load a test DDA file
library(msdata)
fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML",
                 package = "msdata")
sps_dda <- Spectra(fl)

## Define the isolation window based on the MS2 spectra's precursor m/z
## and parameter `tolerance`: isolation window with size 1Da:
pp <- precursorPurity(sps_dda, tolerance = 0.5)

## values for MS1 spectra are NA
head(pp[msLevel(sps_dda) == 1])
#> [1] NA NA NA NA NA NA

head(pp[msLevel(sps_dda) == 2])
#> [1] 1.0000000 0.5894016 1.0000000 1.0000000 0.9223336 1.0000000

## Use the reported isolation window (if defined in the `Spectra`):
filterMsLevel(sps_dda, 2L) |>
    isolationWindowLowerMz() |>
    head()
#> [1] 137.46394  56.44191  89.44488 206.52937  54.50829 120.59900
filterMsLevel(sps_dda, 2L) |>
    isolationWindowUpperMz() |>
    head()
#> [1] 138.46394  57.44191  90.44488 207.52937  55.50829 121.59900

pp_2 <- precursorPurity(sps_dda, useReportedIsolationWindow = TRUE)

head(pp_2[msLevel(sps_dda) == 2])
#> [1] 1.0000000 0.5894016 1.0000000 1.0000000 0.9223336 1.0000000
```
