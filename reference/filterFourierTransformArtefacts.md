# Fast fourier transform artefact filter

The `filterFourierTransformArtefacts()` function removes (Orbitrap) fast
fourier artefact peaks from spectra. Such artefacts (also referred to as
*ripples*) seem to be related to the
[*ringing*](https://en.wikipedia.org/wiki/Ringing_artifacts) phenomenon
and are frequently seen in Orbitrap data as small random mass peaks ~
0.01 Da from a main peak with a very large intensity. See also
[here](https://www.shimadzu.com/an/service-support/technical-support/analysis-basics/tips-ftir/apodization.html)
for more details and information. The data set `fft_spectrum` represents
a
[`Spectra()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
object with a single Orbitrap spectrum with such artefacts (see examples
below).

See also
[`Spectra()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
(section \*Data subsetting, filtering and merging) for the definition of
the function.

## Details

The current implementation iterates through all intensity ordered peaks
in a spectrum and removes all peaks with an m/z within +/-
`halfWindowSize` of the current peak if their intensity is lower than
`threshold` times the current peak's intensity. Additional parameters
`keepIsotopes`, `maxCharge` and `isotopeTolerance` allow to avoid
removing of potential `[13]C` isotope peaks (`maxCharge` being the
maximum charge that should be considered and `isotopeTolerance` the
absolute acceptable tolerance for matching their m/z).

## Author

Jan Stanstrup, Johannes Rainer

## Examples

``` r
library(Spectra)
data(fft_spectrum)

plotSpectra(fft_spectrum)


## Focus on an artefact
plotSpectra(fft_spectrum, xlim = c(264.5, 265.5))

plotSpectra(fft_spectrum, xlim = c(264.5, 265.5), ylim = c(0, 5e6))


fft_spectrum <- filterFourierTransformArtefacts(fft_spectrum)
fft_spectrum
#> MSn data (Spectra) with 1 spectra in a MsBackendDataFrame backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1   367.665       195
#>  ... 33 more variables/columns.
#> Lazy evaluation queue: 1 processing step(s)
#> Processing:
#>  Switch backend from MsBackendMzR to MsBackendDataFrame [Mon Nov 22 14:14:45 2021]
#>  Remove fast fourier artefacts. [Tue Dec 23 13:21:20 2025] 
plotSpectra(fft_spectrum, xlim = c(264.5, 265.5), ylim = c(0, 5e6))


## R code to download/extract the data.

if (FALSE) { # \dontrun{
library(Spectra)
# get orbitrap data
download.file("https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS469/AV_01_v2_male_arm1_juice.mzXML", "AV_01_v2_male_arm1_juice.mzXML")
data <- Spectra("AV_01_v2_male_arm1_juice.mzXML")
extracted_spectrum <- data[195]
} # }
```
