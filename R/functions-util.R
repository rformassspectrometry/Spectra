
#' Simple helper to convert a variable to integer and throw an error if it is
#' not a numeric.
#'
#' @param x `numeric`
#'
#' @author Johannes Rainer
#'
#' @noRd
.as_integer <- function(x) {
    if (is.numeric(x))
        as.integer(x)
    else {
        arg <- deparse(substitute(x))
        stop("Argument ", arg, " should be a numeric or integer", call. = FALSE)
    }
}

#' Helper to convert i to integer index for subsetting.
#'
#' @param i `character` `logical` or `integer` used in `[i]` for subsetting
#'
#' @param length `integer` representing the `length` of the object to be
#'     subsetted
#'
#' @param names `character` with the names (rownames or similar) of the object
#'
#' @return `integer` with the indices
#'
#' @author Johannes Rainer
#'
#' @noRd
.i_to_index <- function(i, length, names = NULL) {
    if (is.character(i)) {
        if (!length(names))
            stop("can not subset by name, object does not have names")
        if (!all(i %in% names))
            stop("not all names present in object")
        i <- match(i, names)
    }
    if (is.logical(i)) {
        if (length(i) != length)
            stop("if i is logical it has to match the length of the object (",
                 length, ")")
        i <- which(i)
    }
    if (is.numeric(i))
        i <- as.integer(i)
    if (is.integer(i) && !all(abs(i) %in% seq_len(length)))
        stop("index out of bounds: index has to be between 1 and ", length)
    i
}

#' Helper function to check the data type even if `x` is an `Rle`.
#'
#' @param x
#'
#' @param class2 `character(1)` specifying the class for which should be tested.
#'
#' @author Johannes Rainer
#'
#' @noRd
.is_class <- function(x, class2) {
    if (is(x, "Rle"))
        is(x@values, class2)
    else is(x, class2)
}

.class_rle <- function(x) {
    if (is(x, "Rle"))
        class(x@values)
    else class(x)
}

#' Function to convert a `numeric`, `logical` or `character` vector into an
#' `Rle` if it will use less memory than the bare vector.
#'
#' @param x vector to convert to Rle.
#'
#' @return `x` or an `Rle` version of `x`
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @importFrom S4Vectors nrun
#'
#' @noRd
.as_rle <- function(x) {
    if (length(x) > 2L && (is.numeric(x) || is.character(x) || is.logical(x))) {
        r <- Rle(x)
        if (nrun(r) < length(x) / 2L)
            return(r)
    }
    x
}

#' Removes zeros from input except the ones that in the direct neighbourhood of
#' non-zero values.
#'
#' @note
#'
#' Copied from MSnbase
#'
#' @param x \code{numeric}, vector to be cleaned
#' @param all \code{logical}, should all zeros be removed?
#' @param na.rm \code{logical}, should NAs removed before looking for zeros?
#' @return logical vector, \code{TRUE} for keeping the value
#' @note The return value for \code{NA} is always \code{FALSE}.
#' @examples
#' x <- c(1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0)
#' #      T, T, F, T, T, T, T, T, T, T, T, F, F
#' r <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
#'        FALSE, FALSE)
#' stopifnot(utils.clean(x) == r)
#' @noRd
utils.clean <- function(x, all=FALSE, na.rm=FALSE) {
  notNA <- !is.na(x)
  notZero <- x != 0 & notNA

  if (all) {
    notZero
  } else if (na.rm) {
    notNA[notNA] <- utils.enableNeighbours(notZero[notNA])
    notNA
  } else {
    utils.enableNeighbours(notZero)
  }
}

#' Switch FALSE to TRUE in the direct neighborhod of TRUE.
#' (used in utils.clean)
#'
#' @param x logical
#' @return logical
#' @examples
#' x <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
#'        FALSE, FALSE)
#' r <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
#'        FALSE, FALSE)
#' stopifnot(utils.enableNeighbours(x) == r)
#' @noRd
utils.enableNeighbours <- function(x) {
  stopifnot(is.logical(x))
  x | c(x[-1], FALSE) | c(FALSE, x[-length(x)])
}

#' @param acquisitionNum `integer` with the acquisition numbers of all scans.
#'
#' @param precursorScanNum `integer` with the precursor scan numbers.
#'
#' @param an integer, acquisitionNum of spectrum of interest (parent and
#' children will be selected)
#'
#' @author Sebastian Gibb
#'
#' @noRd
.filterSpectraHierarchy <- function(acquisitionNum = integer(),
                                    precursorScanNum = integer(), an) {
    if (length(acquisitionNum) != length(precursorScanNum))
        stop("length of 'acquisitionNum' and 'precursorScanNum' have to be ",
             "the same")
    ## we could use recursion which is slow in R
    ## or reformat the adjacency list into a nested tree
    ## list model but most ms data are limited to at most 3 levels and the
    ## filtering isn't done very often, so we use for loops here

    parents <- logical(length(acquisitionNum))

    ## find current scan
    parents[acquisitionNum %in% an] <- TRUE
    children <- parents

    ## find parent scan
    nLastParents <- 0L
    nParents <- 1L
    while (nLastParents < nParents) {
        parents[acquisitionNum %in% precursorScanNum[parents]] <- TRUE
        nLastParents <- nParents
        nParents <- sum(parents)
    }

    ## find children scans
    nLastChildren <- 0L
    nChildren <- 1L
    while (nLastChildren < nChildren) {
        children[precursorScanNum %in% acquisitionNum[children]] <- TRUE
        nLastChildren <- nChildren
        nChildren <- sum(children)
    }
    parents | children
}

#' @title Combine two two-dimensional arrays
#'
#' @importMethodsFrom S4Vectors cbind nrow rownames colnames
#'
#' @description
#'
#' Combine instances of `matrix`, `data.frame` or `DataFrame` objects into a
#' single instance adding eventually missing columns filling them with `NA`s.
#'
#' @note
#'
#' This function might not work if one of the columns contains S4classes.
#'
#' @param ... 2 or more: `matrix`, `data.frame` or `DataFrame`.
#'
#' @return The merged object.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @noRd
.rbind_fill <- function(...) {
    l <- list(...)

    if (length(l) == 1L)
        l <- l[[1L]]

    cl <- vapply(l, class, character(1L))

    stopifnot(all(cl %in% c("matrix", "data.frame", "DataFrame")))

    ## convert matrix to data.frame for easier and equal subsetting and class
    ## determination
    isMatrix <- cl == "matrix"
    l[isMatrix] <- lapply(l[isMatrix], as.data.frame)

    allcl <- unlist(lapply(l, vapply, class, character(1L)))
    allnms <- unique(names(allcl))
    allcl <- allcl[allnms]

    for (i in seq(along=l)) {
        diffcn <- setdiff(allnms, names(l[[i]]))
        if (length(diffcn))
            l[[i]][, diffcn] <- lapply(allcl[diffcn], as, object = NA)
    }
    r <- do.call(rbind, l)

    ## if we had just matrices as input we need to convert our temporary
    ## data.frame back to a matrix
    if (all(isMatrix))
        r <- as.matrix(r)
    r
}

.logging <- function(x, ...) {
    c(x, paste0(..., " [", date(), "]"))
}

setAs("logical", "factor", function(from, to) factor(from))

#' The function aggregates `x` for `toBin` falling into bins defined
#' by `breaks` using the `fun` function.
#'
#' @details
#'
#' This is a combination of the code from the former bin_Spectrum.
#'
#' @param x `numeric` with the values that should be binned.
#'
#' @param toBin `numeric`, same length than `x`, with values to be used for the
#'     binning.
#'
#' @param binSize `numeric(1)` with the size of the bins.
#'
#' @param breaks `numeric` defining the breaks/bins.
#'
#' @param fun `function` to be used to aggregate values of `x` falling into the
#'     bins defined by `breaks`.
#'
#' @return `list` with elements `x` and `mids` being the aggregated values
#'     of `x` for values in `toBin` falling within each bin and the bin mid
#'     points.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @noRd
.bin_values <- function(x, toBin, binSize = 1, breaks = seq(floor(min(toBin)),
                                                            ceiling(max(toBin)),
                                                            by = binSize),
                        fun = max) {
    if (length(x) != length(toBin))
        stop("lengths of 'x' and 'toBin' have to match.")
    fun <- match.fun(fun)
    breaks <- .fix_breaks(breaks, range(toBin))
    nbrks <- length(breaks)
    idx <- findInterval(toBin, breaks)
    ## Ensure that indices are within breaks.
    idx[idx < 1L] <- 1L
    idx[idx >= nbrks] <- nbrks - 1L

    ints <- double(nbrks - 1L)
    ints[unique(idx)] <- unlist(lapply(base::split(x, idx), fun),
                                use.names = FALSE)
    list(x = ints, mids = (breaks[-nbrks] + breaks[-1L]) / 2L)
}

#' Simple function to ensure that breaks (for binning) are spaning at least the
#' expected range.
#'
#' @param brks `numeric` with *breaks* such as calculated by `seq`.
#'
#' @param rng `numeric(2)` with the range of original numeric values on which
#'     the breaks were calculated.
#'
#' @noRd
.fix_breaks <- function(brks, rng) {
    ## Assuming breaks being sorted.
    if (brks[length(brks)] <= rng[2])
        brks <- c(brks, max((rng[2] + 1e-6),
                            brks[length(brks)] + mean(diff(brks))))
    brks
}

#' @title Match values given a user defined tolerance
#'
#' @description
#'
#' `matchApprox` allows to match numeric values/vectors accepting differences
#' defined by parameter `tolerance`. The function supports `tolerance` of
#' length 1 to use the same tolerance value for all comparison, or of length
#' equal `length(x)` to use value-specific tolerance (e.g. for matching with
#' parts-per-million *ppm*). See examples for more details.
#'
#' @note
#'
#' Similar to the `all.equal.numeric` function a value of
#' `sqrt(.Machine$double.eps)` is added to each `tolerance`.
#'
#' @param x `numeric` with values to match against `table`.
#'
#' @param table `numeric` with the values to be matched against.
#'
#' @param nomatch `integer(1)` with the value to be returned in case when no
#'     match is found.
#'
#' @param tolerance `numeric` defining the maximal acceptable difference between
#'     values in `x` and `table` to consider values to be equal (i.e. match).
#'     Can be either of length `1` or `length(x)`. See description and examples
#'     for details.
#'
#' @return `integer` of length equal to `length(x)` with the indices in `table`
#'     where values in `x` match.
#'
#' @importFrom Rcpp evalCpp
#'
#' @useDynLib Spectra
#'
#' @author Johannes Rainer
#'
#' @export
#'
#' @examples
#'
#' ## Define two vectors to match
#' x <- c(1.11, 45.02, 556.45)
#' y <- c(3.01, 34.12, 45.021, 46.1, 556.449)
#'
#' ## Elements do not match exactly
#' match(x, y)
#' matchApprox(x, y)
#'
#' ## Using a single tolerance value
#' matchApprox(x, y, tolerance = 0.01)
#'
#' ## Using a value-specific tolerance accepting differences of 20 ppm
#' matchApprox(x, y, tolerance = x * 20 / 1e6)
#'
#' ## Same using 50ppm
#' matchApprox(x, y, tolerance = x * 50 / 1e6)
#'
#' ## Compare to MSnbase::relaxedMatch and MSnbase::matchPeaks
#' MSnbase:::relaxedMatch(x, y, tolerance = x * 20 / 1e6, relative = FALSE)
#' ## relaxedMatch adds the ppm to the table, not x - which seems unintuitive
#' ## to me.
#'
#' ## Speed:
#' library(microbenchmark)
#' microbenchmark(
#' matchApprox(x, y, tolerance = x * 50 / 1e6),
#' MSnbase:::relaxedMatch(x, y, tolerance = y * 20 / 1e6, relative = FALSE),
#' MALDIquant::match.closest(x, y, tolerance = y * 20 / 1e6)
#' )
#'
#' ## Real data:
#' library(msdata)
#' fl <- system.file("TripleTOF-SWATH/PestMix1_SWATH.mzML", package ="msdata")
#' sps <- Spectra(fl, source = MsBackendMzR())
#'
#' x <- sps$mz[[198]]
#' y <- sps$mz[[199]]
#' microbenchmark(
#' matchApprox(x, y, tolerance = x * 40 / 1e6),
#' MSnbase:::relaxedMatch(x, y, tolerance = y * 40 / 1e6, relative = FALSE),
#' MALDIquant::match.closest(x, y, tolerance = y * 40 / 1e6)
#' )
#' ## Compare to MALDIquant match.closest
matchApprox <- function(x, table = numeric(), nomatch = NA_integer_,
                        tolerance = 0.0) {
    if (!is.integer(nomatch))
        stop("'nomatch' must be an integer")
    if (!is.numeric(tolerance) || tolerance < 0.0)
        stop("'tolerance' must be a positive number")
    x_len <- length(x)
    tolerance <- tolerance + sqrt(.Machine$double.eps)
    if (length(tolerance) != x_len)
        tolerance <- rep(tolerance[1], x_len)
    match_approx_cpp(x, table, nomatch, tolerance)
}


#' Group a sorted `numeric` of m/z values from consecutive scans by ion assuming
#' that the variation between m/z values for the same ion in consecutive scan
#' is much lower than the difference between m/z values within one scan.
#'
#' @param x `numeric` with ordered and combined m/z values from consecutive
#'     scans.
#'
#' @param mzd `numeric(1)` with the m/z difference below which m/z values are
#'     grouped together. If not provided the `.estimate_mz_scattering` function
#'     is used to estimate it.
#'
#' @param ppm `numeric(1)` defining an optional ppm. `mzd` will be increased
#'     by the ppm of the m/z to allow m/z dependent peak groups.
#'
#' @return `integer` of same length than `x` grouping m/z values.
#'
#' @author Sebastian Gibb, Johannes Rainer
#'
#' @noRd
.group_mz_values <- function(x, mzd, ppm = 0) {
    mzdiff <- diff(x)
    if (ppm > 0)
        cumsum(c(0L, mzdiff >= (mzd + x[-length(x)] * ppm / 1e6))) + 1L
    else
        cumsum(c(0L, mzdiff >= mzd)) + 1L
}

#' @title Align/collocate sorted numeric vectors
#'
#' @description
#'
#' Collocate (align) two (increasingly ordered!) numeric vectors re-sizing them
#' to equal lengths filling non-matching elements with `NA`. Parameters
#' `tolerance` and `ppm` allow to loosen the stringency of the value matching:
#' values with an absolute difference `<=` the value of `tolerance` plus value
#' times `ppm` divided by 1,000,000 are considered to be matching.
#'
#' @note
#'
#' Similar to the `match` call only the **first** matching values are
#' considered for multiple matches.
#'
#' @param x `numeric` of **increasingly ordered** values.
#'
#' @param y `numeric` of **increasingly ordered** values.
#'
#' @param tolerance `numeric(1)`: accepted maximal difference between
#'     matching elements in `x` and `y`.
#'
#' @param ppm `numeric(1)` (parts-per-million) allowing to define a
#'     value-dependent tolerance.
#'
#' @param index `logical(1)`: if `TRUE` the indices are returned rather then
#'     the actual values.
#'
#' @return
#'
#' `list` with two elements `x` and `y` representing the aligned vectors `x`
#' and `y` (or their indices if `index = TRUE`).
#'
#' @author Johannes Rainer
#'
#' @seealso
#'
#' [matchApprox] for approximate matching.
#'
#' @noRd
#'
#' @examples
#'
#' x <- c(2, 5, 6)
#' y <- c(1, 2, 3, 5, 7)
#'
#' ## Align/collocate the two vectors. Their length will be the length of
#' ## unique values
#' .collocate_ordered_vectors(x, y)
#'
#' x <- c(1.11, 45.02, 556.45)
#' y <- c(3.01, 34.12, 45.021, 46.1, 556.449)
#'
#' ## Align/collocate the two vectors - no matching elements.
#' .collocate_ordered_vectors(x, y)
#' .collocate_ordered_vectors(x, y, index = TRUE)
#'
#' ## With tolerance
#' .collocate_ordered_vectors(x, y, tolerance = 0.01)
#' .collocate_ordered_vectors(x, y, tolerance = 0.01, index = TRUE)
#'
#' ## With ppm
#' .collocate_ordered_vectors(x, y, ppm = 20)
#' .collocate_ordered_vectors(x, y, ppm = 20, index = TRUE)
#'
#' .collocate_ordered_vectors(x, y, ppm = 40)
#' .collocate_ordered_vectors(x, y, ppm = 40, index = TRUE)
.collocate_ordered_vectors <- function(x, y, tolerance = 0, ppm = 0,
                                    index = FALSE) {
    ordr <- sort(c(x, y))
    grp <- .group_mz_values(ordr, mzd = tolerance + sqrt(.Machine$double.eps),
                            ppm = ppm)
    x_idx <- y_idx <- rep(NA_integer_, grp[length(grp)])
    if (index) {
        x_idx[grp[match(x, ordr)]] <- seq_along(x)
        y_idx[grp[match(y, ordr)]] <- seq_along(y)
    } else {
        x_idx[grp[match(x, ordr)]] <- x
        y_idx[grp[match(y, ordr)]] <- y
    }
    list(x = x_idx, y = y_idx)
}
