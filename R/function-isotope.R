#' Title Isotope identification in spectra
#'
#' @description
#' Given a spectrum (i.e. a peak matrix with m/z and intensity values)
#' the function identifies groups of peaks that correspond to isotopes
#'
#' @param x `matrix` with spectrum data (columns `mz` and `intensity`).
#'
#' @param isotopeDefinition `matrix` with isotopes definition (columns `mzd`,
#' `min_intensity` and `max_intensity`). This matrix has to have rows ordered 
#' according to the column `mzd` values in order for each returned group to be 
#' ordered.
#' 
#' @param tolerance numeric representing the tolerance for the relaxed matching 
#' of m/z values of peaks
#'
#' @param ppm numeric(1) representing a relative, value-specific
#' parts-per-million (PPM) tolerance for the relaxed matching of m/z values of
#' peaks
#'
#' @param seedMz numeric ordered vector containing m/z values. If provided,
#' the function checks if there are peaks in x whose m/z match them. If so,
#' it looks for isotope groups on this subsets of peaks.
#'
#' @return list of vectors. Each vector in the returned list contains the
#' indexes of the rows in `x` that match a isotope.
#'
#' @details The function iterates over the rows of `x`. For the `i`-th m/z
#' value in `x` it computes m/z of would-be isotope peaks by adding the `mzd`s
#' in `isotopeDefinition`. Possible matches of those with the other m/z values
#' in `x` are then searched. If any is found and the ratio between the
#' intensities of the corresponding peaks and the `i`-th peak is in between two
#' given values (specified in the `min_intensity` and `max_intensity` column of 
#' isotopeDefinition`), their indexes are grouped together. Those indexes are 
#' excluded from the set of indexes that are searched for further groups.
#'
#' @examples
#'
.isotope_peaks <- function(x, isotopeDefinition = isotopeDefinition(), 
                           tolerance = 0, ppm = 20, seedMz = numeric()) {
  lst <- list()
  to_test <- rep(TRUE, nrow(x))
  to_test[which(x[, 2] <= 0)] <- FALSE
  if (length(seedMz)) 
    idxs <- na.omit(MsCoreUtils::closest(seedMz, x[which(to_test), 1], 
                                         tolerance = tolerance, ppm = ppm, 
                                         duplicates = "closest"))
  else idxs <- which(to_test)
  for (i in idxs)
  {
    if(to_test[i])
    {
      to_test[i] <- FALSE
      tmp <- x[i, 1] + isotopeDefinition[, 1]
      wtt <- which(to_test)
      cls <- MsCoreUtils::closest(tmp, x[wtt, 1], tolerance = tolerance, 
                                  ppm = ppm)
      int_ok <- .is_isotope_intensity_range(x[, 2][wtt[cls]], x[i, 2], 
                                            isotopeDefinition)
      if (length(int_ok) > 0)
      {
        lst[[length(lst) + 1]] <- c(i, wtt[cls][int_ok])
        to_test[wtt[cls][int_ok]] <- FALSE
      }
    }
  }
  lst
}



#' Title Checking the intensity
#'
#' @param x intensity of the matching peaks. x has length equal to the number 
#' of rows of isotopeDefinition. The i-th element of x represent the intensity 
#' associated to the i-th mzd in isotopeDefinition if any or NA. 
#' @param intensity  intensity of the main peak
#' @param isotopeDefinition isotope definition data.frame
#'
#' @return indexes of the intensities in x that are part of a isotopic group
#' @export
#'
#' @examples
.is_isotope_intensity_range <- function(x, intensity, isotopeDefinition) {
  which(x >= isotopeDefinition[, "min_intensity"] * intensity & 
          x <= isotopeDefinition[, "max_intensity"] * intensity)
}


#' Title Isotope definition
#'
#' @description
#' Creates the definition of isotopes
#'
#' @param charge specifies whether single or double charged ions are expected
#'
#' @return `data.frame` with columns `mzd`, `min_intensity` and `max_intensity`

# isotopeDefinition <- function(charge = c(1L, 2L)) {
#   data.frame(mz = c(1.0033, 0.99...), intensity = ....)
#
# }
