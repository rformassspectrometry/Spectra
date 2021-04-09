#' Title Isotope identification in spectra
#'
#' @description
#' Given a spectrum (i.e. a peak matrix with m/z and intensity values)
#' the function identifies groups of peaks that correspond to isotopes
#'
#' @param x `matrix` with spectrum data (columns `mz` and `intensity`).
#'
#' @param isotopeDefinition `matrix` with isotopes definition (columns `mzd`
#' and `maxint`). This matrix has to have rows ordered according to the column
#' `mzd` values in order for each returned group to be ordered.
#'
#' @param ppm numeric(1) representing a relative, value-specific
#' parts-per-million (PPM) tolerance for the relaxed matching of m/z values of
#' peaks
#'
#' @param seed_mz numeric ordered vector containing m/z values. If provided,
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
#' intensities of the corresponding peaks and the `i`-th peak is smaller than a
#' threshold (specified in the `maxint` column of isotopeDefinition`), 
#' their indexes are grouped together. Those indexes are excluded from the set 
#' of indexes that are searched for further groups.
#'
#' @examples
#'
.isotope_peaks <- function(x, isotopeDefinition = isotopeDefinition(), ppm = 20,
                           seed_mz = NULL) {
  lst <- list()
  to_test <- rep(TRUE, nrow(x))
  to_test[which(x[, 2] <= 0)] <- FALSE
  if (is.null(seed_mz)) idxs <- which(to_test)
  else idxs <- na.omit(MsCoreUtils::closest(seed_mz, x[which(to_test), 1], 
                                            tolerance = 0, ppm = ppm, 
                                            duplicates = "closest"))
  for (i in idxs)
  {
    if(to_test[i])
    {
      to_test[i] <- FALSE
      tmp <- x[i, 1] + isotopeDefinition[, 1]
      wtt <- which(to_test)
      cls <- MsCoreUtils::closest(tmp, x[wtt, 1], tolerance = 0, ppm = ppm)
      int_ok <- which(x[, 2][wtt[cls]] < isotopeDefinition[, 2] * x[i, 2])
      if (length(int_ok) > 0)
      {
        lst[[length(lst) + 1]] <- c(i, wtt[cls][int_ok])
        to_test[wtt[cls][int_ok]] <- FALSE
      }
    }
  }
  lst
}

#' Title Isotope definition
#'
#' @description
#' Creates the definition of isotopes
#'
#' @param charge specifies whether single or double charged ions are expected
#'
#' @return `data.frame` with columns `mzd` and `maxint`

# isotopeDefinition <- function(charge = c(1L, 2L)) {
#   data.frame(mz = c(1.0033, 0.99...), intensity = ....)
#
# }
