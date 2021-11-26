#' @title Oribtrap spectrum
#'
#' @description
#'
#' The `fft_spectrum` variable is a [Spectra()] with a single spectrum
#' measured on an Orbitrap.
#'
#' The code to create/extract the spectrum is shown in the exampleship
#'
#' @name fft_spectrum
#'
#' @examples
#'
#' library(Spectra)
#' data(fft_spectrum)
#'
#' plotSpectra(fft_spectrum)
#'
#' ## R code to download/extract the data.
#'
#' \dontrun{
#' library(Spectra)
#' # get orbitrap data
#' download.file("https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS469/download/4cc5d820-dc5d-4766-8112-7a05f74acef4?file=AV_01_v2_male_arm1_juice.mzXML", "AV_01_v2_male_arm1_juice.mzXML")
#' data <- Spectra("AV_01_v2_male_arm1_juice.mzXML")
#' extracted_spectrum <- data[195]
#' }
NULL

#' @title Fast fourier transform artefact filter
#'
#' @aliases filterFourierTransformArtefacts
#'
#' @description
#'
#' The `filterFourierTransformArtefacts` function removes (Orbitrap) fast
#' fourier artefact peaks from spectra. Such artefacts (also referred to as
#' *rippples*) seem to be related to the
#' [*ringing*](https://en.wikipedia.org/wiki/Ringing_artifacts) phenomenon and
#' are frequently seen in Orbitrap data as small random mass peaks ~ 0.01 Da
#' from a main peak with a very large intensity. See also
#' [here](https://www.shimadzu.com/an/service-support/technical-support/analysis-basics/tips-ftir/apodization.html)
#' for more details and information.
#'
#' See also [Spectra()] (section *Data subsetting, filtering and merging) for
#' the definition of the function.
#'
#' @details
#'
#' The current implementation iterates through all intensity ordered peaks in a
#' spectrum and removes all peaks with an m/z within +/- `halfWindowSize` of
#' the current peak if their intensity is lower than `threshold` times the
#' current peak's intensity. Additional parameters `keepIsotopes`, `maxCharge`
#' and `isotopeTolerance` allow to avoid removing of potential `[13]C` isotope
#' peaks (`maxCharge` being the maximum charge that should be considered
#' and `isotopeTolerance` the absolute acceptable tolerance for matching
#' their m/z).
#'
#' @author Jan Stanstrup, Johannes Rainer
#'
#' @name filterFourierTransformArtefacts
NULL
