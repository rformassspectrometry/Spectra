
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

#' @importFrom stats quantile
#'
#' @noRd
.isCentroided <- function(pk, k = 0.025, qtl = 0.9) {
    .qtl <- quantile(pk[, 2], qtl)
    x <- pk[pk[, 2] > .qtl, 1]
    quantile(diff(x), 0.25) > k
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
    idx[which(idx >= nbrks)] <- nbrks - 1L

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
