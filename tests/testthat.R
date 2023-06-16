library("testthat")
library("patrick")
library("Spectra")

register(SerialParam())

sciex_file <- normalizePath(
    dir(system.file("sciex", package = "msdata"), full.names = TRUE))
cdf_file <- normalizePath(
    dir(system.file("cdf", package = "msdata"), full.names = TRUE))

sciex_mzr <- backendInitialize(MsBackendMzR(), files = sciex_file)
sciex_pks <- peaksData(sciex_mzr)
fl <- normalizePath(
    dir(system.file("proteomics", package = "msdata"), full.names = TRUE))
tmt_mzr <- backendInitialize(MsBackendMzR(), files = fl[5])

fl <- system.file("TripleTOF-SWATH", "PestMix1_SWATH.mzML",
                  package = "msdata")
sps_dia <- Spectra(fl)

sciex_hd5 <- backendInitialize(MsBackendHdf5Peaks(),
                               data = spectraData(sciex_mzr),
                               hdf5path = tempdir())

test_suite <- system.file("test_backends", "test_MsBackend",
                          package = "Spectra")
be <- sciex_mzr[1:10]
test_dir(test_suite, stop_on_failure = TRUE)

be <- backendInitialize(MsBackendDataFrame(), spectraData(be))
test_dir(test_suite, stop_on_failure = TRUE)

be <- backendInitialize(MsBackendMemory(), spectraData(be))
test_dir(test_suite, stop_on_failure = TRUE)

## be <- sciex_hd5[1:10]
## test_dir(test_suite, stop_on_failure = TRUE)

test_check("Spectra")
