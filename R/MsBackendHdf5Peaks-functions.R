#' @include hidden_aliases.R
NULL

.valid_ms_backend_mod_count <- function(x, y) {
    if (length(x) != length(y))
        "Different number of source files and modification counters."
    else
        NULL
}

#' @rdname MsBackend
#'
#' @export MsBackendHdf5Peaks
MsBackendHdf5Peaks <- function() {
    if (!requireNamespace("rhdf5", quietly = TRUE))
        stop("The use of 'MsBackendHdf5Peaks' requires package 'rhdf5'. ",
             "Please install with 'Biobase::install(\"rhdf5\")'")
    new("MsBackendHdf5Peaks")
}

#' @description
#'
#' - check that lengths of `x` and `files` matches.
#' - check that `x` contains unique names.
#' - check that `x` exist.
#' - check that `x` are in the correct format.
#'
#' @param x `character` with the file names of the Hdf5 files.
#'
#' @author Johannes Rainer
#'
#' @noRd
.valid_h5files <- function(x) {
    msg <- NULL
    if (length(x) != length(unique(x)))
        msg <- c(msg, paste0("No duplicated HDF5 files allowed"))
    if (!all(file.exists(x)))
        msg <- c(msg, paste0("File(s) ", paste(x[!file.exists(x)],
                                               collapse = ", "),
                             " do not exist"))
    msg <- c(msg, unlist(lapply(x, .valid_h5peaks_file), use.names = FALSE))
    msg
}

#' @description
#'
#' Initializes a **new** hdf5 file.
#'
#' @noRd
.initialize_h5peaks_file <- function(x, modCount = 0L) {
    h5 <- rhdf5::H5Fcreate(x)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    comp_level <- .hdf5_compression_level()
    rhdf5::h5createGroup(h5, "header")
    rhdf5::h5write("Spectra::MsBackendHdf5Peaks", h5, "/header/class",
                   level = comp_level)
    rhdf5::h5write(modCount, h5, "/header/modcount", level = comp_level)
    rhdf5::h5createGroup(h5, "peaks")
}

.valid_h5peaks_file <- function(x) {
    msg <- NULL
    fid <- try(.Call("_H5Fopen", x, 0L, PACKAGE = "rhdf5"))
    if (is(fid, "try-error"))
        return(paste0("File ", x, " is not a Hdf5 file"))
    on.exit(invisible(.Call("_H5Fclose", fid, PACKAGE = "rhdf5")))
    cls <- try(.h5_read_bare(fid, "/header/class"))
    if (is(cls, "try-error"))
        msg <- "Wrong Hdf5 file format: element /header/class not found"
    else {
        if (cls != "Spectra::MsBackendHdf5Peaks")
            msg <- "Wrong Hdf5 file format: unexpected class"
    }
    modc <- try(.h5_read_bare(fid, "/header/modcount"))
    if (is(modc, "try-error")) {
        msg <- c(msg, paste0("Wrong Hdf5 file format: element /header/modcount",
                             " not found"))
    } else {
        if (!is.integer(modc))
            msg <- c(msg, paste0("Wrong Hdf5 file format: modification counter is",
                                 " not an integer"))
    }
    if (length(msg)) msg
    else NULL
}

#' This is based on Mike Smith's function (issue #395) that directly accesses
#' the data without validation and file checking.
#'
#' @param file the ID as returned by `_H5Fopen`.
#'
#' @param name `character` defining the data set to be read.
#'
#' @return the imported data set (in most cases a `matrix`).
#'
#' @author Mike Smith, Johannes Rainer
#'
#' @md
#'
#' @noRd
.h5_read_bare <- function(file, name = "") {
    did <- .Call("_H5Dopen", file, name, NULL, PACKAGE = "rhdf5")
    res <- .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, FALSE,
                 PACKAGE = "rhdf5")
    invisible(.Call("_H5Dclose", did, PACKAGE = "rhdf5"))
    res
}

.hdf5_compression_level <- function() {
    getOption("HDF5_COMPRESSION_LEVEL", 3L)
}

#' Read peak data from an hdf5 file and return a list of `matrix` objects with
#' columns `"mz"` and `"intensity"`.
#'
#' @param x `character`: the name of the HDF5 file.
#'
#' @param scanIndex `integer` with the scan indices/ids.
#'
#' @param modCount `integer`: the modification counter stored in the R object
#'     This is checked against one stored in the hdf5 file and an error is thrown
#'     if they differ, i.e. if the object was copied and the hdf5 files modified
#'     by the original/source object.
#'
#' @return list of `matrix` objects with the peak data.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.h5_read_peaks <- function(x = character(), scanIndex = integer(),
                           modCount = 0L) {
    if (length(x) != 1)
        stop("'x' should have length 1")
    if (!length(scanIndex))
        return(list(matrix(ncol = 2, nrow = 0,
                           dimnames = list(character(), c("mz", "intensity")))))
    requireNamespace("rhdf5", quietly = TRUE)
    fid <- .Call("_H5Fopen", x, 0L, PACKAGE = "rhdf5")
    on.exit(invisible(.Call("_H5Fclose", fid, PACKAGE = "rhdf5")))
    h5modCount <- .h5_read_bare(fid, "/header/modcount")
    if (h5modCount != modCount)
        stop("The data in the hdf5 files associated with this object appear ",
             "to have changed! Please see the Notes section in ?MsBackend ",
             "for more information.")
    base::lapply(paste0("/spectra/", scanIndex), function(z, file) {
        res <- .h5_read_bare(name = z, file = file)
        colnames(res) <- c("mz", "intensity")
        res
    }, file = fid)
}

#' Write spectra of an mzML file to a group within a hdf5 file. Depending on
#' the value of `prune` the hdf5 group will be deleted before writing the data.
#' This is required as datasets in a hdf5 file can not be *updated* if their
#' dimensions don't match.
#'
#' @param x `list` of `matrix` objects with columns `"mz"` and `"intensity"`.
#'
#' @param scanIndex `integer` with the scan index for the peaks of each
#'     spectrum.
#'
#' @param h5file `character(1)` with the name of the hdf5 file into which the
#'     data should be written.
#'
#' @param modCount `integer`, modification counter, incremented by one for each
#'     write
#'
#' @param prune `logical(1)` whether the group should be deleted before writing
#'     the data (see description above for more details).
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.h5_write_peaks <- function(x = list(), scanIndex = integer(), h5file,
                              modCount = 0L, prune = TRUE) {
    if (length(x) != length(scanIndex))
        stop("lengths of 'x' and 'scanIndex' have to match")
    if (any(duplicated(scanIndex)))
        stop("no duplicated values in 'scanIndex' allowed")
    requireNamespace("rhdf5", quitely = TRUE)
    h5 <- rhdf5::H5Fopen(h5file)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    comp_level <- .hdf5_compression_level()
    x <- force(x)
    if (prune && rhdf5::H5Lexists(h5, "/spectra"))
        rhdf5::h5delete(h5, "/spectra")
    if (!rhdf5::H5Lexists(h5, "/spectra"))
        rhdf5::h5createGroup(h5, "/spectra")
    spids <- paste0("/spectra/", scanIndex)
    for (i in seq_along(x)) {
        xi <- x[[i]]
        if (!is.matrix(xi))
            stop("Peak data has to be provided as a matrix but I got ",
                 class(xi))
        if (ncol(xi) != 2 || !all(c("mz", "intensity") %in% colnames(xi)))
            stop("A peak matrix should have two columns named \"mz\" and ",
                 "\"intensity\"")
        if (is.unsorted(xi[, 1]))
            stop("m/z values have to be ordered")
        rhdf5::h5write(xi, h5, name = spids[i], level = comp_level)
    }
    rhdf5::h5write(modCount, h5, "/header/modcount", level = comp_level)
}
