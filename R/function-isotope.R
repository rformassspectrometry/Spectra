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
#' @param duplicates character(1), how to handle duplicated matches. Has to be
#' one of c("keep", "closest", "remove"). No abbreviations allowed.
#'
#' @return list of vectors. Each vector in the returned list contains the
#' indexes of the rows in `x` that match a isotope.
#'
#' @details The function iterates over the rows of `x`. For the `i`-th m/z
#' value in `x` it computes m/z of would-be isotope peaks by adding the `mzd`s
#' in `isotopeDefinition`. Possible matches of those with the other m/z values
#' in `x` are then searched. If any is found and the ratio between the
#' intensities of the corresponding peaks and the `i`-th peak is smaller than a
#' threshold (specified in the `maxint` column of
#' `isotopeDefinition`), their indexes are grouped together.
#'
#' It can happen that some index is present in more than one group. The user can
#' decide what to do in this case by choosing the parameter `duplicates`.
#' `duplicates = "keep"` leaves the groups as they are.
#' `duplicates = "closest"` removes from the list groups for which another group
#'  which shares some indexes with it exists and constitutes a better match.
#' `duplicates = "remove"` removes from the list groups for which another group
#'  which shares some indexes with it exists
#'
#' @examples
#'
.isotope_peaks <- function(x, isotopeDefinition = isotopeDefinition(), ppm = 20,
                           seed_mz = NULL, duplicates = "keep") {

  if(duplicates != "keep" && duplicates != "closest" && duplicates != "remove")
    stop("'duplicates' has to be one of \"keep\", \"closest\" or \"remove\".")

  lst <- list()
  idxs_o <- which(x[, 2] > 0)
  x_o <- x
  x <- x[idxs_o, ]
  if (is.null(seed_mz)) idxs <- 1:nrow(x)
  else idxs <- na.omit(MsCoreUtils::closest(seed_mz, x[, 1], tolerance = 0,
                                            ppm = ppm, duplicates = "closest"))
  sgd <- NULL
  for (i in idxs)
  {
    append_to_current <- FALSE
    for (j in 1:nrow(isotopeDefinition))
    {
      tmp <- x[i, 1] + isotopeDefinition[j, 1]
      cls <- MsCoreUtils::closest(tmp, x[-(1:i), 1], tolerance = 0, ppm = ppm)
      if(!is.na(cls) &  x[, 2][cls + i] < isotopeDefinition[j, 2]*x[i, 2])
      {
        len <- length(lst)
        if(!append_to_current)
        {
          lst[[len + 1]] <- c(idxs_o[i], idxs_o[cls + i])
          sgd[len + 1] <- tmp - x[cls + i, 1]
          append_to_current <- TRUE
        }
        else
        {
          lst[[len]] <- c(lst[[len]], idxs_o[cls + i])
          sgd[len] <- sgd[len] + tmp - x[cls + i, 1]
        }
      }
    }
  }

  if(length(lst) > 1 && duplicates == "closest")
  {
    maxmzd <- max(isotopeDefinition[,1])
    k <- 1
    while (k < length(lst))
    {
      kk <- k+1
      k_removed <- FALSE
      while(kk<=length(lst) && !k_removed &&
            x_o[max(lst[[k]]), 1] + maxmzd > x_o[lst[[kk]][1], 1]*(1 - ppm*10^(-6)))
      {
        if(length(intersect(lst[[k]], lst[[kk]])) > 0)
        {
          if(sgd[[k]]/length(lst[[k]]) > sgd[[kk]]/length(lst[[kk]]))
          {
            lst[[k]] <- NULL
            k_removed <- TRUE
          }
          else if(sgd[[k]]/length(lst[[k]]) < sgd[[kk]]/length(lst[[kk]]))
          {
            lst[[kk]] <- NULL
          }
          else if(length(lst[[kk]]) > length(lst[[k]]))
          {
            lst[[k]] <- NULL
            k_removed <- TRUE
          }
          else
          {
            lst[[kk]] <- NULL
          }
        }
        else kk <- kk + 1
      }
      if(k_removed == FALSE)
        k <- k + 1
    }
  }

  if (length(lst) > 1 && duplicates == "remove")
  {
    tabunlst <- table(unlist(lst))
    dupl <- as.numeric(names(tabunlst[which(tabunlst > 1)]))
    rem <- NULL
    for(k in seq_along(lst))
    {
      if (length(intersect(lst[[k]], dupl)) > 0)
        rem <- c(rem, k)
    }
    lst[rem] <- NULL
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
