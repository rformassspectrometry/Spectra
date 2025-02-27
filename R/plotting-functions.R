#' @title Plotting Spectra
#'
#' @aliases plotSpectraMirror
#'
#' @description
#'
#' [Spectra()] can be plotted with one of the following functions
#'
#' - `plotSpectra()`: plots each spectrum in its separate plot by splitting
#'   the plot area into as many panels as there are spectra.
#'
#' - `plotSpectraOverlay()`: plots all spectra in `x` **into the same** plot (as
#'   an overlay).
#'
#' - `plotSpectraMirror()`: plots a pair of spectra as a *mirror plot*.
#'   Parameters `x` and `y` both have to be a `Spectra` of length 1. Matching
#'   peaks (considering `ppm` and `tolerance`) are highlighted. See
#'   [MsCoreUtils::common()] for details on peak matching. Parameters
#'   `matchCol`, `matchLty`, `matchLwd` and `matchPch` allow to customize
#'   how matching peaks are indicated.
#'
#' @param x a [Spectra()] object. For `plotSpectraMirror()` it has to be an
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
#'     spectrum's MS level and retention time (in seconds) is used.
#'
#' @param col color to be used to draw the peaks. Should be either of length 1,
#'     or equal to the number of spectra (to plot each spectrum in a different
#'     color) or be a `list` with colors for each individual peak in each
#'     spectrum.
#'
#' @param labels allows to specify a label for each peak. Can be a `character`
#'     with length equal to the number of peaks, or, ideally, a `function` that
#'     uses one of the `Spectra`'s variables (see examples below).
#'     `plotSpectraMirror()` supports only `labels` of type *function*.
#'
#' @param labelCex `numeric(1)` giving the amount by which the text should be
#'     magnified relative to the default. See parameter `cex` in [par()].
#'
#' @param labelSrt `numeric(1)` defining the rotation of the label. See
#'     parameter `srt` in [text()].
#'
#' @param labelCol color for the label(s).
#'
#' @param labelAdj see parameter `adj` in [text()].
#'
#' @param labelPos see parameter `pos` in [text()].
#'
#' @param labelOffset see parameter `offset` in [text()].
#'
#' @param axes `logical(1)` whether (x and y) axes should be drawn.
#'
#' @param frame.plot `logical(1)` whether a box should be drawn around the
#'     plotting area.
#'
#' @param ppm for `plotSpectraMirror()`: m/z relative acceptable difference (in
#'     ppm) for peaks to be considered matching (see [MsCoreUtils::common()]
#'     for more details).
#'
#' @param tolerance for `plotSpectraMirror()`: absolute acceptable difference of
#'     m/z values for peaks to be considered matching (see
#'     [MsCoreUtils::common()] for more details).
#'
#' @param matchCol for `plotSpectraMirror()`: color for matching peaks.
#'
#' @param matchLwd for `plotSpectraMirror()`: line width (`lwd`) to draw
#'     matching peaks. See [par()] for more details.
#'
#' @param matchLty for `plotSpectraMirror()`: line type (`lty`) to draw matching
#'     peaks. See [par()] for more details.
#'
#' @param matchPch for `plotSpectraMirror()`: point character (`pch`) to label
#'     matching peaks. Defaults to `matchPch = 16`, set to `matchPch = NA` to
#'     disable. See [par()] for more details.
#'
#' @param y for `plotSpectraMirror()`: `Spectra` object of length 1 against
#'     which `x` should be plotted against.
#'
#' @param asp for `plotSpectra()`: the target ratio (columns / rows) when
#'     plotting mutliple spectra (e.g. for 20 spectra use `asp = 4/5` for 4
#'     columns and 5 rows or `asp = 5/4` for 5 columns and 4 rows; see
#'     [grDevices::n2mfrow()] for details).
#'
#' @param ... additional parameters to be passed to the [plot.default()]
#'     function.
#'
#' @return These functions create a plot.
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
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
#' ## Plot one spectrum.
#' plotSpectra(sp[1])
#'
#' ## Plot both spectra.
#' plotSpectra(sp)
#'
#' ## Define a color for each peak in each spectrum.
#' plotSpectra(sp, col = list(c(1, 2, 3, 4, 5), 1:4))
#'
#' ## Color peaks from each spectrum in different colors.
#' plotSpectra(sp, col = c("green", "blue"))
#'
#' ## Label each peak with its m/z.
#' plotSpectra(sp, labels = function(z) format(unlist(mz(z)), digits = 4))
#'
#' ## Rotate the labels.
#' plotSpectra(sp, labels = function(z) format(unlist(mz(z)), digits = 4),
#'     labelPos = 2, labelOffset = 0.1, labelSrt = -30)
#'
#' ## Add a custom annotation for each peak.
#' sp$label <- list(c("", "A", "B", "C", "D"),
#'     c("Frodo", "Bilbo", "Peregrin", "Samwise"))
#' ## Plot each peak in a different color
#' plotSpectra(sp, labels = function(z) unlist(z$label),
#'     col = list(1:5, 1:4))
#'
#' ## Plot a single spectrum specifying the label.
#' plotSpectra(sp[2], labels = c("A", "B", "C", "D"))
#'
#'
#' #### --------------------------------------------- ####
#' ##                plotSpectraOverlay                 ##
#'
#' ## Plot both spectra overlaying.
#' plotSpectraOverlay(sp)
#'
#' ## Use a different color for each spectrum.
#' plotSpectraOverlay(sp, col = c("#ff000080", "#0000ff80"))
#'
#' ## Label also the peaks with their m/z if their intensity is above 15.
#' plotSpectraOverlay(sp, col = c("#ff000080", "#0000ff80"),
#'     labels = function(z) {
#'         lbls <- format(mz(z)[[1L]], digits = 4)
#'         lbls[intensity(z)[[1L]] <= 15] <- ""
#'         lbls
#'     })
#' abline(h = 15, lty = 2)
#'
#' ## Use different asp values
#' plotSpectra(sp, asp = 1/2)
#' plotSpectra(sp, asp = 2/1)
#'
#' #### --------------------------------------------- ####
#' ##                plotSpectraMirror                  ##
#'
#' ## Plot two spectra against each other.
#' plotSpectraMirror(sp[1], sp[2])
#'
#' ## Label the peaks with their m/z
#' plotSpectraMirror(sp[1], sp[2],
#'     labels = function(z) format(mz(z)[[1L]], digits = 3),
#'     labelSrt = -30, labelPos = 2, labelOffset = 0.2)
#' grid()
#'
#' ## The same plot with a tolerance of 0.1 and using a different color to
#' ## highlight matching peaks
#' plotSpectraMirror(sp[1], sp[2],
#'     labels = function(z) format(mz(z)[[1L]], digits = 3),
#'     labelSrt = -30, labelPos = 2, labelOffset = 0.2, tolerance = 0.1,
#'     matchCol = "#ff000080", matchLwd = 2)
#' grid()
NULL

#' @rdname spectra-plotting
#'
#' @importFrom graphics par
#' @importFrom grDevices n2mfrow
#'
#' @export plotSpectra
plotSpectra <- function(x, xlab = "m/z", ylab = "intensity", type = "h",
                        xlim = numeric(), ylim = numeric(),
                        main = character(), col = "#00000080",
                        labels = character(), labelCex = 1, labelSrt = 0,
                        labelAdj = NULL, labelPos = NULL, labelOffset = 0.5,
                        labelCol = "#00000080", asp = 1, ...) {
    if (!length(main))
        main <- paste0("MS", msLevel(x), " RT: ", round(rtime(x), 1))
    nsp <- length(x)
    if (nsp == 1)
        col <- list(col)
    if (length(col) != nsp)
        col <- rep(col[1], nsp)
    if (length(main) != nsp)
        main <- rep(main[1], nsp)
    if (nsp > 1)
        par(mfrow = n2mfrow(nsp, asp = asp))
    if (length(labels)) {
        if (is.function(labels)) {
            labels <- labels(x)
        }
        if (length(labels) != length(x))
            stop("Please provide a list of annotations of 'length == length(x)'.")
    } else {labels <- NULL}
    
    for (i in seq_len(nsp))
        .plot_single_spectrum(x[i], xlab = xlab, ylab = ylab, type = type,
                              xlim = xlim, ylim = ylim, main = main[i],
                              col = col[[i]], labels = labels[[i]],
                              labelCex = labelCex, labelSrt = labelSrt,
                              labelAdj = labelAdj, labelPos = labelPos,
                              labelOffset = labelOffset, labelCol = labelCol,
                              ...)
}

#' @rdname spectra-plotting
#'
#' @export plotSpectraOverlay
plotSpectraOverlay <- function(x, xlab = "m/z", ylab = "intensity",
                               type = "h", xlim = numeric(),
                               ylim = numeric(),
                               main = paste(length(x), "spectra"),
                               col = "#00000080", labels = character(),
                               labelCex = 1, labelSrt = 0,
                               labelAdj = NULL, labelPos = NULL,
                               labelOffset = 0.5, labelCol = "#00000080",
                               axes = TRUE, frame.plot = axes, ...) {
    nsp <- length(x)
    if (nsp == 1)
        col <- list(col)
    if (length(col) != nsp)
        col <- rep(col[1], nsp)
    if (!length(xlim))
        xlim <- range(unlist(mz(x)), na.rm = TRUE)
    if (!length(ylim))
        ylim <- c(0, max(unlist(intensity(x)), na.rm = TRUE))
    dev.hold()
    on.exit(dev.flush())
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    if (axes) {
        axis(side = 1, ...)
        axis(side = 2, ...)
    }
    if (frame.plot)
        box(...)
    title(main = main, xlab = xlab, ylab = ylab, ...)
    if (length(labels)) {
        if (is.function(labels)) {
            labels <- labels(x)
        }
        if (length(labels) != length(x))
            stop("Please provide a list of annotations of 'length == length(x)'.")
    } else {labels <- NULL}
    for (i in seq_len(nsp))
        .plot_single_spectrum(x[i], add = TRUE, type = type, col = col[[i]],
                              labels = labels[[i]], labelCex = labelCex,
                              labelSrt = labelSrt, labelAdj = labelAdj,
                              labelPos = labelPos, labelOffset = labelOffset,
                              labelCol = labelCol, ...)
}

#' @rdname spectra-plotting
#'
#' @importFrom MsCoreUtils common
#'
#' @importFrom graphics abline
#'
#' @exportMethod plotSpectraMirror
setMethod(
    "plotSpectraMirror", "Spectra",
    function(x, y, xlab = "m/z", ylab = "intensity",
             type = "h", xlim = numeric(),
             ylim = numeric(), main = character(),
             col = "#00000080", labels = character(),
             labelCex = 1, labelSrt = 0,
             labelAdj = NULL, labelPos = NULL,
             labelOffset = 0.5, labelCol = "#00000080",
             axes = TRUE, frame.plot = axes, ppm = 20,
             tolerance = 0, matchCol = "#80B1D3", matchLwd = 1,
             matchLty = 1, matchPch = 16, ...) {
        if (length(x) != 1 || length(y) != 1)
            stop("'x' and 'y' have to be of length 1")
        if (length(col) != 2)
            col <- rep(col[1], 2)
        if (!length(xlim))
            suppressWarnings(
                xlim <- range(unlist(mz(x)), unlist(mz(y)), na.rm = TRUE))
        if (!length(ylim))
            suppressWarnings(
                ylim <- c(-1, 1) * max(unlist(intensity(x)),
                                       unlist(intensity(y)),
                                       na.rm = TRUE))
        if (any(is.infinite(xlim)))
            xlim <- c(0, 0)
        if (any(is.infinite(ylim)))
            ylim <- c(0, 0)
        dev.hold()
        on.exit(dev.flush())
        plot.new()
        plot.window(xlim = xlim, ylim = ylim)
        ## Stop if variable modifications are used
        ## Will need to be removed once plotSpectra accepts variable modifications
        ## See issue: https://github.com/rformassspectrometry/Spectra/issues/346
        if (length(labels)) {
            if (is.function(labels)) {
                labels <- c(labels(x), labels(y))
            } else {
                if (length(labels) != length(x))
                    stop("This Error occurs either because\n1) Annotations are not of length 2\n2) Variable modifications are not yet supported.")
            }
            l <- c(labels[[1]], labels[[2]])
            wdths <- max(strwidth(l, cex = labelCex)) / 2
            usr_lim <- par("usr")
            ylim[1L] <- ylim[1L] - wdths *
                diff(usr_lim[3:4]) / diff(usr_lim[1:2])
            ylim[2L] <- -ylim[1L]
            xlim[1L] <- xlim[1L] - wdths
            xlim[2L] <- xlim[2L] + wdths
            plot.window(xlim = xlim, ylim = ylim, ...)
        } else {labels <- NULL}
        if (axes) {
            axis(side = 1, ...)
            axis(side = 2, ...)
        }
        if (frame.plot)
            box(...)
        title(main = main, xlab = xlab, ylab = ylab, ...)
        ## Find common peaks
        x_data <- peaksData(x)[[1L]]
        y_data <- peaksData(y)[[1L]]
        .plot_single_spectrum(x, add = TRUE, type = type, col = col[[1L]],
                              labels = labels[[1L]], labelCex = labelCex,
                              labelSrt = labelSrt, labelAdj = labelAdj,
                              labelPos = labelPos, labelOffset = labelOffset,
                              labelCol = labelCol, ...)
        idx <- which(common(x_data[, "mz"], y_data[, "mz"],
                            tolerance = tolerance, ppm = ppm))
        if (length(idx)) {
            plot.xy(xy.coords(x_data[idx, "mz"], x_data[idx, "intensity"]),
                    type = "h", col = matchCol, lwd = matchLwd, ...)
            plot.xy(xy.coords(x_data[idx, "mz"], x_data[idx, "intensity"]),
                    type = "p", col = matchCol, pch = matchPch, ...)
        }
        if (length(labelPos) && labelPos == 1)
            labelPos <- 3
        if (length(labelPos) && labelPos == 3)
            labelPos <- 1
        labelSrt <- -1 * labelSrt
        .plot_single_spectrum(y, add = TRUE, type = type, col = col[[1L]],
                              labels = labels[[2L]], labelCex = labelCex,
                              labelSrt = labelSrt, labelAdj = labelAdj,
                              labelPos = labelPos, labelOffset = labelOffset,
                              orientation = -1, labelCol = labelCol, ...)
        idx <- which(common(y_data[, "mz"], x_data[, "mz"],
                            tolerance = tolerance, ppm = ppm))
        if (length(idx)) {
            plot.xy(xy.coords(y_data[idx, "mz"], -y_data[idx, "intensity"]),
                    type = "h", col = matchCol, lwd = matchLwd, ...)
            plot.xy(xy.coords(y_data[idx, "mz"], -y_data[idx, "intensity"]),
                    type = "p", col = matchCol, pch = matchPch, ...)
        }
        abline(h = 0)
    })

#' @description
#'
#' Plot a single spectrum (m/z on x against intensity on y) with the optional
#' possibility to label the individual peaks.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @importFrom graphics axis box plot.new plot.window plot.xy strwidth
#'
#' @importFrom graphics text title
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
#' .plot_single_spectrum(sp,
#'     labels = function(z) format(mz(z)[[1]], digits = 5),
#'     labelPos = 2, labelOffset = 0.1, labelSrt = -30)
#' grid()
#'
#' @noRd
.plot_single_spectrum <- function(x, xlab = "m/z", ylab = "intensity",
                                  type = "h", xlim = numeric(),
                                  ylim = numeric(),
                                  main = paste("RT", round(rtime(x), 1)),
                                  col = "#00000080", labels = character(),
                                  labelCol = col, labelCex = 1, labelSrt = 0,
                                  labelAdj = NULL, labelPos = NULL,
                                  labelOffset = 0.5, add = FALSE,
                                  axes = TRUE, frame.plot = axes,
                                  orientation = 1, ...) {
    v <- peaksData(x)[[1L]]
    mzs <- v[, "mz"]
    ints <- orientation * v[, "intensity"]
    if (!length(xlim))
        suppressWarnings(xlim <- range(mzs, na.rm = TRUE))
    if (!length(ylim))
        suppressWarnings(
            ylim <- range(orientation * c(0, max(abs(ints), na.rm = TRUE))))
    if (any(is.infinite(xlim)))
        xlim <- c(0, 0)
    if (any(is.infinite(ylim)))
        ylim <- c(0, 0)
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
        usr_lim <- par("usr")
        ylim[2L] <- ylim[2L] + wdths * diff(usr_lim[3:4]) / diff(usr_lim[1:2])
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
