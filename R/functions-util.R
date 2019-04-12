
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
