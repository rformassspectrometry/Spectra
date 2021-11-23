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
