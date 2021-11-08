.onLoad <- function(libname, pkgname) {
    options(HDF5_COMPRESSION_LEVEL = 3L)
    ## options(READ_HEADER = Sys.info()["sysname"] == "Darwin")
    options(READ_HEADER = FALSE)        # Also on macOS it's now better to disable
}
