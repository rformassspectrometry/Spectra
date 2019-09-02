
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

.logging <- function(x, ...) {
    c(x, paste0(..., " [", date(), "]"))
}

setAs("logical", "factor", function(from, to) factor(from))
