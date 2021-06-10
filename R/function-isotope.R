#' Title Isotope identification in spectra
#'
#' @description
#' Given a spectrum (i.e. a peak matrix with m/z and intensity values)
#' the function identifies groups of peaks that correspond to isotopes
#'
#' @param x `matrix` with spectrum data (columns `mz` and `intensity`).
#'
#' @param isotopeDefinition `matrix` with isotopes definition (columns 
#' `subst_name`, `subst_degree`, `md`, `min_slope`, `max_slope`). This matrix 
#' has to have rows ordered in such a way that the column `md` is sorted in 
#' increasing order.
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
#' it looks for isotope groups related to this subset of peaks.
#' 
#' @param charge numeric(1) representing the charge of the ionized compounds 
#'
#' @return list of vectors. Each vector in the returned list contains the
#' indexes of the rows in `x` that match a certain isotope group found by 
#' the function.
#'
#' @details The function iterates over the peaks (rows) in x. Firstly, it checks 
#' the presence of peaks whose m/z difference with respect to the current peak 
#' matches (the matching can be relaxed via the parameters `ppm` and 
#' `tolerance`) certain m/z differences related to mass differences of 
#' certain substitutions in `substDefinition` and the chosen `charge`. Then, if 
#' any such peak is found, the function checks if also their intensity is 
#' compatible with them being part of a isotopic group and if so they are 
#' grouped together. When some peaks are grouped together their indexes are 
#' excluded from the set of indexes that are searched for further groups. 
#'
#' @examples
#'
.isotope_peaks <- function(x, substDefinition = substDefinition(), 
                           tolerance = 0, ppm = 20, seedMz = numeric(),
                           charge = 1) {
  lst <- list()
  to_test <- rep(TRUE, nrow(x))
  to_test[which(x[, 2] <= 0)] <- FALSE
  if (length(seedMz)) 
    idxs <- na.omit(MsCoreUtils::closest(seedMz, x[which(to_test), 1], 
                                         tolerance = tolerance, ppm = ppm, 
                                         duplicates = "closest"))
  else idxs <- which(to_test)
  mzd <- substDefinition[, "md"]/charge
  for (i in idxs)
  {
    if(to_test[i])
    {
      to_test[i] <- FALSE
      tmp <- x[i, 1] + mzd
      wtt <- which(to_test)
      cls <- MsCoreUtils::closest(tmp, x[wtt, 1], tolerance = tolerance, 
                                  ppm = ppm, duplicates = "closest")
      int_ok <- .is_isotope_intensity_range(x[, 2][wtt[cls]], x[i, 1] * charge,
                                            x[i, 2], substDefinition)
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
#' of rows of substDefinition. The i-th element of x represent the intensity 
#' associated to the peak whose m/z difference is associated to the i-th mzd in 
#' substDefinition if any or NA. 
#' @param m mass of the current (assumed monoisotopic) peak
#' @param intensity  intensity of the current (assumed monoisotopic) peak
#' @param substDefinition substitutions definition data.frame
#'
#' @return indexes of the intensities in x that are part of a isotopic group
#' @export
#'
#' @examples
.is_isotope_intensity_range <- function(x, m, intensity, substDefinition) {
  R_min <- (m * substDefinition[, "min_slope"]) ^ substDefinition[, "subst_degree"]
  R_max <- (m * substDefinition[, "max_slope"]) ^ substDefinition[, "subst_degree"]
  which(x >= R_min * intensity & x <= R_max * intensity)
}


#' Title Isotope definition
#'
#' @description
#' Creates the definition of isotopes
#'
#' @param charge specifies whether single or double charged ions are expected
#'
#' @return `data.frame` with columns `mzd`, `min_ratio` and `max_ratio`

# isotopeDefinition <- function(charge = c(1L, 2L)) {
#   data.frame(mz = c(1.0033, 0.99...), intensity = ....)
#
# }