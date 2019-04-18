
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
