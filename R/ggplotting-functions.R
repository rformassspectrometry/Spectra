#' @title Plotting Spectra with ggplot
#'
#' @aliases ggplotSpectraMirror
#'
#' @description
#'
#' [Spectra()] can be plotted with one of the following functions
#'
#' - `ggplotSpectra()`: plots each spectrum in its separate plot by splitting
#'   the plot area into as many panels as there are spectra.
#'
#' - `ggplotSpectraOverlay()`: plots all spectra in `x` **into the same** plot
#'   (as an overlay).
#'
#' - `ggplotSpectraMirror()`: plots a pair of spectra as a *mirror plot*.
#'   Parameters `x` and `y` both have to be a `Spectra` of length 1. Matching
#'   peaks (considering `ppm` and `tolerance`) are highlighted. See
#'   [MsCoreUtils::common()] for details on peak matching. Parameters
#'   `matchCol`, `matchLty`, `matchLwd` and `matchPch` allow to customize
#'   how matching peaks are indicated.
#'
#' @param x a [Spectra()] object. For `ggplotSpectraMirror()` it has to be an
#'     object of length 2.
#'
#' @param xlab `character(1)` with the label for the x-axis (by default
#'     `xlab = "m/z"`).
#'
#' @param ylab `character(1)` with the label for the y-axis (by default
#'     `ylab = "intensity"`).
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
#' @param labels allows to specify a label for each peak. Needs to be a `list()`
#'     with length equal to the number of spectra (each element of the list
#'     being a `character()` with length equal to the number of peaks for that
#'     spectrum), or, ideally, a `function` that uses one of the `Spectra`'s
#'     variables (see examples below). `ggplotSpectraMirror()` supports only
#'     `labels` of type *function*.
#'
#' @param labelAngle `numeric(1)` defining the rotation of the label. See
#'     parameter `angle` in *Aesthetics* of [ggplot2::geom_text()].
#'
#' @param labelCol color for the label(s).
#'
#' @param labelSize size of the label(s).
#'
#' @param labelVjust vertical justification of the label(s).
#'
#' @param labelHjust horizontal justification of the label(s).
#'
#' @param axes `logical(1)` whether (x and y) axes should be drawn.
#'
#' @param frame.plot `logical(1)` whether a box should be drawn around the
#'     plotting area.
#'
#' @param ppm for `ggplotSpectraMirror()`: m/z relative acceptable difference
#'     (in ppm) for peaks to be considered matching (see [MsCoreUtils::common()]
#'     for more details).
#'
#' @param tolerance for `ggplotSpectraMirror()`: absolute acceptable difference
#'     of m/z values for peaks to be considered matching (see
#'     [MsCoreUtils::common()] for more details).
#'
#' @param matchCol for `ggplotSpectraMirror()`: color for matching peaks.
#'
#' @param matchLwd for `ggplotSpectraMirror()`: line width to draw matching
#'     peaks.
#'
#' @param matchPch for `ggplotSpectraMirror()`: point character to label
#'     matching peaks. Defaults to `matchPch = 16`, set to `matchPch = NA` to
#'     disable.
#'
#' @param matchCex for `ggplotSpectraMirror()`: point size to draw matching
#'     peaks.
#'
#' @param y for `ggplotSpectraMirror()`: `Spectra` object of length 1 against
#'     which `x` should be plotted against.
#'
#' @param asp aspect ratio of the plot(s) Default: `1/2`.
#'
#' @param interactive `logical(1)` return the interactive ggplot based on
#'     ggiraph.
#'
#' @return These functions create a ggplot.
#'
#' @author Gabriele Tomè, Johannes Rainer
#'
#' @name spectra-ggplotting
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
#' ##                   ggplotSpectra                   ##
#'
#' ## Plot one spectrum.
#' ggplotSpectra(sp[1])
#'
#' ## Plot both spectra.
#' ggplotSpectra(sp)
#'
#' ## Define a color for each peak in each spectrum.
#' ggplotSpectra(sp, col = list(c(1, 2, 3, 4, 5), 1:4))
#'
#' ## Color peaks from each spectrum in different colors.
#' ggplotSpectra(sp, col = c("green", "blue"))
#'
#' ## Label each peak with its m/z.
#' ggplotSpectra(sp, labels = function(z) lapply(mz(z), format, digits = 4))
#'
#' ## Rotate the labels.
#' ggplotSpectra(sp, labels = function(z) lapply(mz(z), format, digits = 4),
#'     labelAngle = -30, labelHjust = 0)
#'
#' ## Add a custom annotation for each peak.
#' sp$label <- list(c("", "A", "B", "C", "D"),
#'     c("Frodo", "Bilbo", "Peregrin", "Samwise"))
#'
#' ## Plot each peak in a different color
#' ggplotSpectra(sp, labels = sp$label,
#'     col = list(1:5, 1:4))
#'
#' ## Plot a single spectrum specifying the label.
#' ggplotSpectra(sp[2], labels = list(c("A", "B", "C", "D")))
#'
#'
#' #### --------------------------------------------- ####
#' ##                ggplotSpectraOverlay               ##
#'
#' ## Plot both spectra overlaying.
#' ggplotSpectraOverlay(sp)
#'
#' ## Use a different color for each spectrum.
#' ggplotSpectraOverlay(sp, col = c("#ff000080", "#0000ff80"))
#'
#' ## Label also the peaks with their m/z if their intensity is above 15.
#' ggplotSpectraOverlay(sp, col = c("#ff000080", "#0000ff80"),
#' labels = function(z) {
#'     lapply(seq_along(mz(z)), function(i) {
#'             lbls <- format(mz(z)[[i]], digits = 4)
#'             lbls[intensity(z)[[i]] <= 15] <- ""
#'             lbls
#'      })
#'  }) + ggplot2::geom_hline(yintercept = 15, linetype = 2)
#'
#' ## Use different asp values
#' ggplotSpectra(sp, asp = 1/3)
#' ggplotSpectra(sp, asp = 2/1)
#'
#' #### --------------------------------------------- ####
#' ##                ggplotSpectraMirror                ##
#'
#' ## Plot two spectra against each other.
#' ggplotSpectraMirror(sp[1], sp[2])
#'
#' ## Label the peaks with their m/z
#' ggplotSpectraMirror(sp[1], sp[2],
#'     labels = function(z) list(format(mz(z)[[1L]], digits = 3)),
#'     labelAngle = -30, labelHjust = 0)
#'
#' ## The same ggplot with a tolerance of 0.1 and using a different color to
#' ## highlight matching peaks
#' ggplotSpectraMirror(sp[1], sp[2],
#'     labels = function(z) list(format(mz(z)[[1L]], digits = 3)),
#'     labelAngle = -30, labelHjust = 0, tolerance = 0.1,
#'     matchCol = "#ff000080", matchLwd = 2, matchCex = 8)
NULL

#' @rdname spectra-ggplotting
#'
#' @importFrom ggiraph girafe
#'
#' @export ggplotSpectra
ggplotSpectra <- function(x, xlab = "m/z", ylab = "intensity",
                        xlim = numeric(), ylim = numeric(),
                        main = character(), col = "#00000080",
                        labels = character(), labelCol = col, labelSize = 5,
                        labelAngle = 0, labelVjust = -0.2, labelHjust = 0.5, asp = 0.5, axes = TRUE, frame.plot = axes,
                        interactive = FALSE){
    if (!length(main))
        main <- paste0("MS", msLevel(x), " RT: ", round(rtime(x), 1))
    nsp <- length(x)
    if (length(col) != nsp)
        col <- rep(col[1], nsp)
    if (length(main) != nsp)
        main <- rep(main[1], nsp)

    if (length(labels)) {
        if (is.function(labels))
            labels <- labels(x)
        if (is.character(labels))
            labels <- list(labels)
        if (length(labels) != length(x))
            stop("Please provide a list of annotations of length equal to 'x'.")
    } else {labels <- NULL}

    if (!interactive){
        .plot_single_spectrum_ggplot(x, xlab = xlab, ylab = ylab,
                                xlim = xlim, ylim = ylim, main = main,
                                col = col, labels = labels,
                                labelCol = labelCol, labelSize = labelSize,
                                labelAngle = labelAngle,
                                labelVjust = labelVjust,
                                labelHjust = labelHjust,
                                asp = asp, axes = axes, frame.plot = frame.plot)
    } else {
        gg <- .plot_single_spectrum_ggplot_interactive(x, xlab = xlab,
                            ylab = ylab, xlim = xlim, ylim = ylim, main = main,
                            col = col, labels = labels, labelCol = labelCol,
                            labelSize = labelSize, labelAngle = labelAngle,
                            labelVjust = labelVjust, labelHjust = labelHjust,
                            asp = asp, axes = axes, frame.plot = frame.plot)
        girafe(gg)
    }
}

#' @rdname spectra-ggplotting
#'
#' @export ggplotSpectraOverlay
ggplotSpectraOverlay <- function(x, xlab = "m/z", ylab = "intensity",
                               xlim = numeric(), ylim = numeric(),
                               main = paste(length(x), "spectra"),
                               col = "#00000080", labels = character(),
                               labelCol = col, labelSize = 5,
                               labelAngle = 0, labelVjust = -0.2,
                               labelHjust = 0.5, asp = 0.5, axes = TRUE,
                               frame.plot = axes, interactive = FALSE) {
    nsp <- length(x)
    if (length(col) != nsp)
        col <- rep(col[1], nsp)
    if (!length(xlim))
        xlim <- c(min(unlist(mz(x)))*0.9, max(unlist(mz(x)))*1.1)
    if (!length(ylim))
        ylim <- c(0, max(unlist(intensity(x))*1.1, na.rm = TRUE))

    if (length(labels)) {
        if (is.function(labels))
            labels <- labels(x)
        if (is.character(labels))
            labels <- list(labels)
        if (length(labels) != length(x))
            stop("Please provide a list of annotations of length equal to 'x'.")
    } else {labels <- NULL}

    if (!interactive){
        .plot_single_spectrum_ggplot(x, add = TRUE, xlab = xlab, ylab = ylab,
                            xlim = xlim, ylim = ylim, main = main, col = col,
                            labels = labels, labelCol = labelCol,
                            labelSize = labelSize, labelAngle = labelAngle,
                            labelVjust = labelVjust, labelHjust = labelHjust,
                            asp = asp, axes = axes, frame.plot = frame.plot)
    } else {
        gg <- .plot_single_spectrum_ggplot_interactive(x, add = TRUE,
                            xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
                            main = main, col = col, labels = labels,
                            labelCol = labelCol, labelSize = labelSize,
                            labelAngle = labelAngle, labelVjust = labelVjust,
                            labelHjust = labelHjust, asp = asp, axes = axes,
                            frame.plot = frame.plot)
        girafe(gg)
    }
}

#' @rdname spectra-ggplotting
#'
#' @exportMethod ggplotSpectraMirror
setMethod(
    "ggplotSpectraMirror", "Spectra",
    function(x, y, xlab = "m/z", ylab = "intensity",
            xlim = numeric(), ylim = numeric(), main = character(),
            col = "#00000080", labels = character(), labelCol = col,
            labelSize = 5, labelAngle = 0, labelVjust = -0.2, labelHjust = 0.5,
            axes = TRUE, frame.plot = axes, ppm = 20, tolerance = 0,
            matchCol = "#80B1D3", matchCex = 5, matchPch = 16,
            matchLwd = 0.5, asp = 0.5, interactive = FALSE) {
        if (length(x) != 1 || length(y) != 1)
            stop("'x' and 'y' have to be of length 1")
        if (length(col) != 2)
            col <- rep(col[1], 2)
        if (!length(main))
            main <- list(paste0("MS", msLevel(x), " RT: ", round(rtime(x), 1)),
                         paste0("MS", msLevel(y), " RT: ", round(rtime(y), 1)))

        ## Stop if variable modifications are used
        ## Will need to be removed once plotSpectra accepts variable modifications
        ## See issue: https://github.com/rformassspectrometry/Spectra/issues/346
        if (length(labels)) {
            if (is.function(labels)) {
                x_labels <- labels(x)
                y_labels <- labels(y)
                if (is.character(x_labels)) {
                    x_labels <- list(x_labels)
                    y_labels <- list(y_labels)
                }
                labels <- c(x_labels, y_labels)
            } else {
                if (length(labels) != 2)
                    stop("This Error occurs either because\n1) Annotations are not of length 2\n2) For 'labelFragments', variable modifications are not yet supported.")
            }
            l <- c(labels[[1]], labels[[2]])
        } else {labels <- NULL}

        if (!interactive){
            .plot_single_spectrum_ggplot(x, y = y, add = TRUE, xlab = xlab,
                            ylab = ylab, xlim = xlim, ylim = ylim,
                            col = col, main = main, ppm = ppm,
                            tolerance = tolerance, labels = labels,
                            labelCol = labelCol, labelSize = labelSize,
                            labelAngle = labelAngle, labelVjust = labelVjust,
                            labelHjust = labelHjust, matchCol = matchCol,
                            matchCex = matchCex, matchPch = matchPch,
                            matchLwd = matchLwd, axes = axes,
                            frame.plot = frame.plot)
        } else {
            gg <- .plot_single_spectrum_ggplot_interactive(x, y = y, add = TRUE,
                                xlab = xlab, ylab = ylab, xlim = xlim,
                                ylim = ylim, col = col, main = main, ppm = ppm,
                                tolerance = tolerance, labels = labels,
                                labelCol = labelCol, labelSize = labelSize,
                                labelAngle = labelAngle,
                                labelVjust = labelVjust,
                                labelHjust = labelHjust, matchCol = matchCol,
                                matchCex = matchCex, matchPch = matchPch,
                                matchLwd = matchLwd, axes = axes,
                                frame.plot = frame.plot)
            girafe(gg)
        }
    })


#' @description
#'
#' Plot a single spectrum (m/z on x against intensity on y) with the optional
#' possibility to label the individual peaks.
#'
#' @author Gabriele Tomè, Johannes Rainer
#'
#' @importFrom ggplot2 ggplot geom_bar scale_color_identity labs theme_bw theme
#' @importFrom ggplot2 geom_text facet_wrap scale_y_continuous xlim ylim ggtitle
#' @importFrom ggplot2 aes element_blank geom_hline geom_point
#'
#' @importFrom MsCoreUtils common
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
#' .plot_single_spectrum_ggplot(sp, main = "hello")
#' .plot_single_spectrum_ggplot(sp, frame.plot = FALSE)
#'
#' .plot_single_spectrum_ggplot(sp, col = 1:5)
#' .plot_single_spectrum_ggplot(sp, col = "red")
#'
#' .plot_single_spectrum_ggplot(sp, labels = 1:5, col = 1:5)
#'
#' .plot_single_spectrum_ggplot(sp, labels = format(mz(sp)[[1]], digits = 5),
#'     labelHjust = 0, labelAngle = -30)
#'
#' @noRd
.plot_single_spectrum_ggplot <- function(x, y = NULL,
                                    xlab = "m/z", ylab = "intensity",
                                    xlim = numeric(), ylim = numeric(),
                                    main = paste("RT", round(rtime(x), 1)),
                                    col = "#00000080", labels = character(),
                                    labelCol = col, labelSize = 5,
                                    labelAngle = 0, labelVjust = -0.2,
                                    labelHjust = 0.5, add = FALSE,
                                    orientation = 1, bs = 16, asp = 0.5,
                                    tolerance = 0, ppm = 20,
                                    matchCol = "#80B1D3", matchPch = 16,
                                    matchCex = 5, matchLwd = 0.5,
                                    axes = TRUE, frame.plot = axes) {
    v <- asDataFrame(x)
    v$intensity_orient <- orientation * v[, "intensity"]
    v$rtime <- factor(v$rtime, levels = sort(unique(v$rtime)),
                      labels = main[[1]])

    if (!is.null(y)) {
        v_y <- asDataFrame(y)
        v_y$intensity_orient <- -orientation * v_y[, "intensity"]
        v_y$rtime <- factor(v_y$rtime, levels = sort(unique(v_y$rtime)))
        ## Find common peaks
        v$common <- common(v[, "mz"], v_y[, "mz"],
                            tolerance = tolerance, ppm = ppm)
        v_y$common <- common(v_y[, "mz"], v[, "mz"],
                            tolerance = tolerance, ppm = ppm)

        v <- rbind(v, v_y)
    }
    if (length(labels))
        v$label <- unlist(labels)

    names(col) <- unique(v$rtime)
    col_df <- stack(col)
    if(nrow(col_df) != length(unique(v$rtime))){
        col_df$values <- factor(col_df$values, levels = unique(col_df$values))
        v$color <- col_df$values
    } else {
        names(col_df) <- c("color", "rtime")
        v <- merge(v, col_df)
    }

    gg <- ggplot(v, aes(x = mz, y = intensity_orient)) +
        labs(x = xlab, y = ylab) +
        theme_bw(base_size = bs) +
        theme(legend.position = "none", panel.grid = element_blank(),
              aspect.ratio = asp)

    if(!is.numeric(v$color)) {
        gg <- gg + scale_color_identity()
    }

    if(!add) {
        gg <- gg +
            geom_bar(aes(group = rtime, color = color), width = 0.5,
                    na.rm = TRUE, stat = "identity") +
            theme(strip.background = element_blank()) +
            facet_wrap(.~rtime, scales = "free") +
            scale_y_continuous(limits = c(0, max(v$intensity_orient)*1.15))
    } else if (!is.null(y)) {
        gg <- gg +
            geom_bar(data = v[!v$common, ], aes(group = rtime),
                width = 0.5, na.rm = TRUE,
                stat = "identity") +
            geom_bar(data = v[v$common, ], aes(group = rtime), width = matchLwd,
                color = matchCol, fill = matchCol, na.rm = TRUE,
                stat = "identity") +
            geom_hline(yintercept = 0) +
            geom_point(data = v[v$common, ],
                    aes(x = mz, y = intensity_orient, group = rtime),
                    color = matchCol, size = matchCex, shape = matchPch)
    } else {
        gg <- gg +
            geom_bar(aes(group = rtime, color = color), width = 0.5,
                    na.rm = TRUE, stat = "identity") +
            xlim(xlim) + ylim(ylim) + ggtitle(main)
    }

    if (length(labels)) {
        gg <- gg +
            geom_text(aes(label = label, color = color),
                        angle = labelAngle, vjust = labelVjust,
                        hjust = labelHjust, size = labelSize)
    }

    if (!axes) {
        gg <- gg +
            theme(axis.text = element_blank(), axis.ticks = element_blank(),
                axis.line = element_blank())
    }
    if (!frame.plot) {
        gg <- gg +
            theme(panel.border = element_blank())
    }
    gg
}

#' @description
#'
#' Generate an interactive plot a single spectrum (m/z on x against intensity
#' on y) with the optional possibility to label the individual peaks.
#'
#' @author Gabriele Tomè, Johannes Rainer
#'
#' @importFrom ggplot2 ggplot scale_color_identity labs theme_bw theme aes
#' @importFrom ggplot2 scale_y_continuous xlim ylim ggtitle element_blank
#'
#' @importFrom ggiraph geom_bar_interactive geom_text_interactive
#' @importFrom ggiraph geom_point_interactive geom_hline_interactive
#' @importFrom ggiraph facet_wrap_interactive set_girafe_defaults
#' @importFrom ggiraph opts_zoom opts_tooltip opts_sizing opts_toolbar
#'
#' @importFrom MsCoreUtils common
#'
#' @noRd
.plot_single_spectrum_ggplot_interactive <- function(x, y = NULL,
                                  xlab = "m/z", ylab = "intensity",
                                  xlim = numeric(), ylim = numeric(),
                                  main = paste("RT", round(rtime(x), 1)),
                                  col = "#00000080", labels = character(),
                                  labelCol = col, labelSize = 5,
                                  labelAngle = 0, labelVjust = 0,
                                  labelHjust = -0.1, add = FALSE,
                                  orientation = 1, bs = 16, asp = 0.5,
                                  tolerance = 0, ppm = 20,
                                  matchCol = "#80B1D3", matchPch = 16,
                                  matchCex = 5, matchLwd = 0.5,
                                  axes = TRUE, frame.plot = axes, ...) {
    v <- asDataFrame(x)
    v$intensity_orient <- orientation * v[, "intensity"]
    v$rtime <- factor(v$rtime, levels = sort(unique(v$rtime)),
                      labels = main[[1]])

    if (!is.null(y)) {
        v_y <- asDataFrame(y)
        v_y$intensity_orient <- -orientation * v_y[, "intensity"]
        v_y$rtime <- factor(v_y$rtime, levels = sort(unique(v_y$rtime)),
                        labels = main[[2]])

        v$common <- common(v[, "mz"], v_y[, "mz"],
                            tolerance = tolerance, ppm = ppm)
        v_y$common <- common(v_y[, "mz"], v[, "mz"],
                            tolerance = tolerance, ppm = ppm)

        v <- rbind(v, v_y)
    }
    v$index <- 1:nrow(v)
    if (length(labels))
        v$label <- unlist(labels)

    names(col) <- unique(v$rtime)
    col_df <- stack(col)
    if(nrow(col_df) != length(unique(v$rtime))){
        v$color <- col_df$values
    } else {
        names(col_df) <- c("color", "rtime")
        v <- merge(v, col_df)
    }

    set_girafe_defaults(
        opts_zoom = opts_zoom(min = 1, max = 4),
        opts_tooltip = opts_tooltip(
            css = "padding:3px;background-color:#333333;color:white;"),
        opts_sizing = opts_sizing(rescale = TRUE),
        opts_toolbar = opts_toolbar(saveaspng = TRUE, position = "topright",
                                    delay_mouseout = 5000, fixed = TRUE)
    )

    gg <- ggplot(v, aes(x = mz, y = intensity_orient)) +
        labs(x = xlab, y = ylab) +
        theme_bw(base_size = bs) +
        theme(legend.position = "none", panel.grid = element_blank(),
              aspect.ratio = asp)

    if(!is.numeric(v$color)) {
        gg <- gg + scale_color_identity()
    }

    if(!add) {
        gg <- gg +
            geom_bar_interactive(
                aes(group = rtime, color = color, data_id = index,
                    tooltip = do.call(paste,
                        c(lapply(names(v),
                                 function(x){paste(x, ": ", v[, x])}),
                        sep = "\n"))),
                width = 0.5, na.rm = TRUE, stat = "identity",
                hover_nearest = TRUE) +
            theme(strip.background = element_blank()) +
            facet_wrap_interactive(.~rtime, scales = "free") +
            scale_y_continuous(limits = c(0, max(v$intensity_orient)*1.15))
    } else if (!is.null(y)) {
        gg <- gg +
            geom_bar_interactive(data = v[!v$common, ],
                aes(group = rtime, data_id = index), width = 0.5, na.rm = TRUE,
                hover_nearest = TRUE, stat = "identity") +
            geom_bar_interactive(data = v[v$common, ],
                aes(group = rtime, data_id = index), width = matchLwd,
                color = matchCol, fill = matchCol, na.rm = TRUE,
                stat = "identity", hover_nearest = TRUE) +
            geom_hline_interactive(yintercept = 0) +
            geom_point_interactive(data = v[v$common, ],
                    aes(x = mz, y = intensity_orient, group = rtime,
                        data_id = index,
                        tooltip = do.call(paste,
                            c(lapply(names(v),
                                    function(x){paste(x, ": ", v[, x])}),
                            sep = "\n"))),
                    color = matchCol, size = matchCex, shape = matchPch,
                    hover_nearest = TRUE)
    } else {
        gg <- gg +
            geom_bar_interactive(
                    aes(group = rtime, color = color, data_id = index,
                    tooltip = do.call(paste,
                        c(lapply(names(v),
                                 function(x){paste(x, ": ", v[, x])}),
                        sep = "\n"))), width = 0.5,
                    na.rm = TRUE, stat = "identity", hover_nearest = TRUE) +
            xlim(xlim) + ylim(ylim) + ggtitle(main)
    }

    if (length(labels)) {
        gg <- gg +
            geom_text_interactive(aes(label = label, color = color),
                                    angle = labelAngle, vjust = labelVjust,
                                    hjust = labelHjust, size = labelSize)
    }

    if (!axes) {
        gg <- gg +
            theme(axis.text = element_blank(), axis.ticks = element_blank(),
                  axis.line = element_blank())
    }
    if (!frame.plot) {
        gg <- gg + theme(panel.border = element_blank())
    }
    gg
}
