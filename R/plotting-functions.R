#' @description
#'
#' Plot a single spectrum (m/z on x against intensity on y) with the optional
#' possibility to label the individual peaks.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @examples
#'
#' ints <- c(4.3412, 12, 8, 34, 23.4)
#' mzs <- c(13.453421, 43.433122, 46.6653553, 129.111212, 322.24432)
#'
#' df <- DataFrame(msLevel = 1L, rtime = 123.12)
#' df$mz <- list(mzs)
#' df$intensity <- list(ints)
#' sp <- Spectra(df)
#'
#' .plot_single_spectrum(sp)
#'
#' .plot_single_spectrum(sp, col = 1:5)
#'
#' .plot_single_spectrum(sp, labels = 1:5, col = 1:5)
#'
#' .plot_single_spectrum(sp, labels = format(mz(sp)[[1]], digits = 5),
#'     labelPos = 2, labelOffset = 0.1, labelSrt = -30)
#' grid()
#'
#' .plot_single_spectrum(sp, labels = function(z) format(mz(z)[[1]], digits = 5),
#'     labelPos = 2, labelOffset = 0.1, labelSrt = -30)
#' grid()
#'
#' @noRd
.plot_single_spectrum <- function(x, xlab = "m/z", ylab = "intensity",
                                  type = "h", xlim = numeric(),
                                  ylim = numeric(), main = paste("RT", rtime(x)),
                                  col = "#00000080", labels = character(),
                                  labelCol = col, labelCex = 1, labelSrt = 0,
                                  labelAdj = NULL, labelPos = NULL,
                                  labelOffset = 0.5, add = FALSE, ...) {
    v <- as.list(x)[[1L]]
    mzs <- v[, "mz"]
    ints <- v[, "intensity"]
    if (!length(xlim))
        xlim <- range(mzs, na.rm = TRUE)
    if (!length(ylim))
        ylim <- c(0, max(ints, na.rm = TRUE))
    if (length(labels)) {
        if (!add)
            plot(x = NA, xlim = xlim, ylim = ylim, axes = FALSE)
        if (is.function(labels))
            labels <- labels(x)
        wdths <- max(strwidth(labels, cex = labelCex)) / 2
        ylim[2L] <- ylim[2L] + wdths * diff(ylim) / diff(xlim)
        xlim[1L] <- xlim[1L] - wdths
        xlim[2L] <- xlim[2L] + wdths
    }
    if (!add)
        plot(x = NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
             main = main)
    points(mzs, ints, type = type, col = col, ...)
    if (length(labels))
        text(mzs, ints, labels = labels, adj = labelAdj, pos = labelPos,
             col = labelCol, cex = labelCex, srt = labelSrt,
             offset = labelOffset)
}
