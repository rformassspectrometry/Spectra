#' @include Spectra.R
NULL

setClassUnion("functionOrNull", c("function", "NULL"))

#' @title Calculate Neutral Loss Spectra
#'
#' @name neutralLoss
#'
#' @description
#'
#' This help page lists functions that convert MS/MS spectra to neutral loss
#' spectra. The main function for this is `neutralLoss` and the specific
#' algorithm to be used is defined (and configured) with dedicated *parameter*
#' objects (paramer `param` of the `neutralLoss` function).
#'
#' The parameter objects for the different algorithms are:
#'
#' - `PrecursorMzParam`: calculates neutral loss spectra as in Aisporna
#'   *et al.* 2022 by subtracting the (fragment's) peak m/z value from the
#'   precursor m/z value of each spectrum (precursor m/z - fragment m/z).
#'   Parameter `msLevel` allows to restrict calculation of neutral loss
#'   spectra to specified MS level(s). Spectra from other MS level(s) are
#'   returned as-is. Parameter `filterPeaks` allows to filter fragment peaks
#'   e.g. if their m/z value is larger than the precursor's
#'   (`filterPeaks = "abovePrecursor"`) or below
#'   (`filterPeaks = "belowPrecursor")`. By default, no peak filtering
#'   is performed. See also parameter `filterPeaks` below for details.
#'
#' @note
#'
#' By definition, mass peaks in a `Spectra` object need to ordered by their m/z
#' value (in increasing order). Thus, the order of the peaks in the calculated
#' neutral loss spectra might not be the same than in the original `Spectra`
#' object.
#'
#' Note also that for spectra with a missing precursor m/z empty spectra are
#' returned (i.e. spectra without peaks) since it is not possible to calcualte
#' the neutral loss spectra.
#'
#' @param filterPeaks For `PrecursorMzParam`: `character(1)` or `function`
#'     defining if and how fragment peaks should be filtered before calculation.
#'     Pre-defined options are: `"none"` (keep all peaks), `"abovePrecursor"`
#'     (removes all fragment peaks with an m/z >= precursor m/z),
#'     `"belowPrecursor"` (removes all fragment peaks with an m/z <= precursor
#'     m/z). In addition, it is possible to pass a custom function with this
#'     parameter with arguments `x` (two column peak matrix) and `precursorMz`
#'     (the precursor m/z) that returns the sub-setted two column peak matrix.
#'
#' @param msLevel `integer` defining for which MS level(s) the neutral loss
#'     spectra should be calculated. Defaults to `msLevel = c(2L, NA)` thus,
#'     neutral loss spectra will be calculated for all spectra with MS level
#'     equal to 2 or with missing/undefined MS level. All spectra with a MS
#'     level different than `msLevel` will be returned unchanged.
#'
#' @param object [Spectra()] object with the fragment spectra for which neutral
#'     loss spectra should be calculated.
#'
#' @param param One of the *parameter* objects discussed below.
#'
#' @param ... Currently ignored.
#'
#' @return A [Spectra()] object with calculated neutral loss spectra.
#'
#' @author Johannes Rainer
#'
#' @references
#'
#' Aisporna A, Benton PH, Chen A, Derks RJE, Galano JM, Giera M and Siuzdak G
#' (2022). Neutral Loss Mass Spectral Data Enhances Molecular Similarity
#' Analysis in METLIN. Journal of the American Society for Mass Spectrometry.
#' \doi{10.1021/jasms.1c00343}
#'
#' @examples
#'
#' ## Create a simple example Spectra object with some MS1, MS2 and MS3 spectra.
#' DF <- DataFrame(msLevel = c(1L, 2L, 3L, 1L, 2L, 3L),
#'                 precursorMz = c(NA, 40, 20, NA, 300, 200))
#' DF$mz <- IRanges::NumericList(
#'                       c(3, 12, 14, 15, 16, 200),
#'                       c(13, 23, 39, 86),
#'                       c(5, 7, 20, 34, 50),
#'                       c(5, 7, 9, 20, 100),
#'                       c(15, 53, 299, 300),
#'                       c(34, 56, 100, 200, 204, 309)
#'                   , compress = FALSE)
#' DF$intensity <- IRanges::NumericList(1:6, 1:4, 1:5, 1:5, 1:4, 1:6,
#'                                      compress = FALSE)
#' sps <- Spectra(DF, backend = MsBackendDataFrame())
#'
#' ## Calculate neutral loss spectra for all MS2 spectra, keeping MS1 and MS3
#' ## spectra unchanged.
#' sps_nl <- neutralLoss(sps, PrecursorMzParam(msLevel = 2L))
#' mz(sps)
#' mz(sps_nl)
#'
#' ## Calculate neutral loss spectra for MS2 and MS3 spectra, removing peaks
#' ## with an m/z >= precursorMz
#' sps_nl <- neutralLoss(sps, PrecursorMzParam(
#'     filterPeaks = "abovePrecursor", msLevel = 2:3))
#' mz(sps_nl)
#'
#' ## Empty spectra are returned for MS 2 spectra with undefined precursor m/z.
#' sps$precursorMz <- NA_real_
#' sps_nl <- neutralLoss(sps, PrecursorMzParam())
#' mz(sps_nl)
NULL

#' @importClassesFrom ProtGenerics Param
#'
#' @noRd
setClass("PrecursorMzParam",
         contains = "Param",
         slots = c(
             filterPeaks = "functionOrNull",
             msLevel = "integer"
         ),
         prototype = prototype(
             filterPeaks = NULL,
             msLevel = 2L,
             version = "0.1"
         ),
         validity = function(object) {
           msg <- NULL
           msg
         })

#' @rdname neutralLoss
#'
#' @export
PrecursorMzParam <- function(filterPeaks = c("none", "abovePrecursor",
                                             "belowPrecursor"),
                             msLevel = c(2L, NA_integer_)) {
    filterFun <- NULL
    if (is.character(filterPeaks)) {
        filterPeaks <- match.arg(filterPeaks)
        filterFun <- switch(
            filterPeaks,
            none = function(x, ...) x,
            abovePrecursor = function(x, precursorMz, ...) {
                x[which(x[, "mz"] < precursorMz), , drop = FALSE]
            },
            belowPrecursor = function(x, precursorMz, ...) {
                x[which(x[, "mz"] > precursorMz), , drop = FALSE]
            })
    }
    if (!is.function(filterFun))
        stop("'filterPeaks' should be either one of \"none\", ",
             "\"abovePrecursor\" or \"belowPrecursor\" or a function.")
    if (!is.numeric(msLevel))
        stop("'msLevel' is expected to be of type integer")
    new("PrecursorMzParam", filterPeaks = filterFun,
        msLevel = as.integer(msLevel))
}

#' @exportMethod neutralLoss
#'
#' @rdname neutralLoss
setMethod(
    "neutralLoss", c(object = "Spectra", param = "PrecursorMzParam"),
    function(object, param, ...) {
        .nl <- function(x, spectrumMsLevel, precursorMz, msl,
                        filterPeaks, ...) {
            if (spectrumMsLevel %in% msl) {
                x <- filterPeaks(x, precursorMz, ...)
                x[, "mz"] <- precursorMz - x[, "mz"]
                x <- x[!is.na(x[, "mz"]), , drop = FALSE]
                if (is.unsorted(x[, "mz"]))
                    x <- x[order(x[, "mz"]), , drop = FALSE]
            }
            x
        }
        .fpeaks <- param@filterPeaks
        if (!length(.fpeaks))
            .fpeaks <- function(x, ...) x
        object <- addProcessing(
            object, .nl, spectraVariables = c("msLevel", "precursorMz"),
            msl = param@msLevel, filterPeaks = .fpeaks)
        object@processing <- .logging(
            object@processing, "Neutral loss calculation using ",
            "'PrecursorMzParam'")
        object@metadata <- c(object@metadata, list(param = param))
        object
    })

## .hypothetical_loss <- function() {
##     ## pairwise differences of m/z values. intensity is mean.
##     ## start with largest m/z, difference to all others with smaller m/z.
##     ## peaks matrix: nrow: (n-1) * n/2
##     for (i in seq_len(nrow(pks))) {
##         for(j in (i+1):nrow(pks)) {
##         }
##     }
## }
