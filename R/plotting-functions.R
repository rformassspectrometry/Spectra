#' @title Plotting Spectra
#'
#' @description
#'
#' [Spectra()] can be plotted with one of the following functions
#'
#' - `plotSpectra`: plots each spectrum in its separate plot by splitting
#'   the plot area into as many panels as there are spectra.
#'
#' - `plotSpectraOverlay`: plots all spectra **into the same** plot (as an
#'   overlay).
#'
#' - `plotSpectraMirror`: plots a pair of spectra as a *mirror plot*.
#'
#' @param x a [Spectra()] object. For `plotSpectraMirror` it has to be an
#'     object of length 2.
#'
#' @param xlab `character(1)` with the label for the x-axis (by default
#'     `xlab = "m/z"`).
#'
#' @param ylab `character(1)` with the label for the y-axis (by default
#'     `ylab = "intensity"`).
#'
#' @param type `character(1)` specifying the type of plot. See [plot.default()]
#'     for details. Defaults to `type = "h"` which draws each peak as a line.
#'
#' @param xlim `numeric(2)` defining the x-axis limits. The range of m/z values
#'     are used by default.
#'
#' @param ylim `numeric(2)` defining the y-axis limits. The range of intensity
#'     values are used by default.
#'
#' @param main `character(1)` with the title for the plot. By default the
#'     spectrum's retention time (in seconds) is used.
#'
#' @param col color to be used to draw the peaks. Should be either of length 1,
#'     or equal to the number of spectra (to plot each spectrum in a different
#'     color) or be a `list` with colors for each individual peak in each
#'     spectrum.
#'
#' @param labels allows to specify a label for each peak. Can be a `character`
#'     with length equal to the number of peaks, or, ideally, a `function` that
#'     uses one of the `Spectra`'s variables (see examples below).
#'
#' @param labelCex `numeric(1)` giving the amount by which the text should be
#'     magnified relative to the default. See parameter `cex` in [par()].
#'
#' @param labelSrt `numeric(1)` defining the rotation of the label. See
#'     parameter `srt` in [text()].
#'
#' @param labelAdj see parameter `adj` in [text()].
#'
#' @param labelPos see parameter `pos` in [text()].
#'
#' @param labelOffset see parameter `offset` in [text()].
#'
#' @param ... additional parameters to be passed to the [plot.default()]
#'     function.
#'
#' @author Johannes Rainer
#'
#' @name spectra-plotting
#'
#' @examples
#'
#' ints <- list(c(4.3412, 12, 8, 34, 23.4),
#'     c(8, 25, 16, 32))
#' mzs <- list(c(13.453421, 43.433122, 46.6653553, 129.111212, 322.24432),
#'     c(13.452, 43.5122, 129.112, 322.245))
#'
#' df <- DataFrame(msLevel = c(1L, 1L), rtime = c(123.12, 124))
#' df$mz <- mzs
#' df$intensity <- ints
#' sp <- Spectra(df)
#'
#' #### --------------------------------------------- ####
#' ##                   plotSpectra                     ##
#'
#' ## Plot one spectrum
#' plotSpectra(sp[1])
#'
#' ## Plot both spectra
#' plotSpectra(sp)
#'
#' ## Define a color for each peak in each spectrum
#' plotSpectra(sp, col = list(c(1, 2, 3, 4, 5), 1:4))
#'
#' ## Color peaks from each spectrum in different colors.
#' plotSpectra(sp, col = c("green", "blue"))
#'
#' ## Label each peak with its m/z
#' plotSpectra(sp, labels = function(z) format(unlist(mz(z)), digits = 4))
#'
#' ## Rotate the labels
#' plotSpectra(sp, labels = function(z) format(unlist(mz(z)), digits = 4),
#'     labelPos = 2, labelOffset = 0.1, labelSrt = -30)
#'
#' ## Add a custom annotation for each peak
#' sp$label <- list(c("", "A", "B", "C", "D"),
#'     c("Frodo", "Bilbo", "Peregrin", "Samwise"))
#' ## Plot each peak in a different color
#' plotSpectra(sp, labels = function(z) unlist(z$label),
#'     col = list(1:5, 1:4))
#'
#' ## Plot a single spectrum specifying the label
#' plotSpectra(sp[2], labels = c("A", "B", "C", "D"))
NULL

#' @rdname spectra-plotting
#'
#' @importFrom graphics par
#'
#' @export plotSpectra
plotSpectra <- function(x, xlab = "m/z", ylab = "intensity",
                        type = "h", xlim = numeric(),
                        ylim = numeric(), main = paste("RT", rtime(x)),
                        col = "#00000080", labels = character(),
                        labelCex = 1, labelSrt = 0,
                        labelAdj = NULL, labelPos = NULL,
                        labelOffset = 0.5, ...) {
    nsp <- length(x)
    if (nsp == 1)
        col <- list(col)
    if (length(col) != nsp)
        col <- rep(col[1], nsp)
    if (length(main) != nsp)
        main <- rep(main[1], nsp)
    par(mfrow = c(round(sqrt(nsp)), ceiling(sqrt(nsp))))
    for (i in seq_len(nsp))
        .plot_single_spectrum(x[i], xlab = xlab, ylab = ylab, type = type,
                              xlim = xlim, ylim = ylim, main = main[i],
                              col = col[[i]], labels = labels,
                              labelCex = labelCex, labelSrt = labelSrt,
                              labelAdj = labelAdj, labelPos = labelPos,
                              labelOffset = labelOffset, ...)
}

#' @description
#'
#' Plot a single spectrum (m/z on x against intensity on y) with the optional
#' possibility to label the individual peaks.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @importFrom graphics axis box plot.new plot.window plot.xy strwidth text title
#'
#' @importFrom grDevices dev.flush dev.hold xy.coords
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
#' .plot_single_spectrum(sp, main = "hello")
#' .plot_single_spectrum(sp, bty = "n")
#' .plot_single_spectrum(sp, frame.plot = FALSE)
#'
#' .plot_single_spectrum(sp, col = 1:5)
#'
#' .plot_single_spectrum(sp, labels = 1:5, col = 1:5)
#'
#' .plot_single_spectrum(sp, labels = format(mz(sp)[[1]], digits = 5),
#'     labelPos = 2, labelOffset = 0.1, labelSrt = -30)
#' grid()
#' .plot_single_spectrum(sp, col = "red", type = "p", add = TRUE)
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
                                  labelOffset = 0.5, add = FALSE,
                                  axes = TRUE, frame.plot = axes, ...) {
    v <- as.list(x)[[1L]]
    mzs <- v[, "mz"]
    ints <- v[, "intensity"]
    if (!length(xlim))
        xlim <- range(mzs, na.rm = TRUE)
    if (!length(ylim))
        ylim <- c(0, max(ints, na.rm = TRUE))
    if (!add) {
        dev.hold()
        on.exit(dev.flush())
        plot.new()
        plot.window(xlim = xlim, ylim = ylim)
    }
    if (length(labels)) {
        if (is.function(labels))
            labels <- labels(x)
        wdths <- max(strwidth(labels, cex = labelCex)) / 2
        ylim[2L] <- ylim[2L] + wdths * diff(ylim) / diff(xlim)
        xlim[1L] <- xlim[1L] - wdths
        xlim[2L] <- xlim[2L] + wdths
        if (!add)
            plot.window(xlim = xlim, ylim = ylim, ...)
    }
    if (!add) {
        if (axes) {
            axis(side = 1, ...)
            axis(side = 2, ...)
        }
        if (frame.plot)
            box(...)
        title(main = main, xlab = xlab, ylab = ylab, ...)
    }
    plot.xy(xy.coords(mzs, ints), type = type, col = col, ...)

    if (length(labels))
        text(mzs, ints, labels = labels, adj = labelAdj, pos = labelPos,
             col = labelCol, cex = labelCex, srt = labelSrt,
             offset = labelOffset)
}
