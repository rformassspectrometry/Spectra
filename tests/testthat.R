library("testthat")
library("Spectra")

sciex_file <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
cdf_file <- dir(system.file("cdf", package = "msdata"), full.names = TRUE)

sciex_mzr <- backendInitialize(MsBackendMzR(), files = sciex_file)
sciex_pks <- peaks(sciex_mzr)
fl <- dir(system.file("proteomics", package = "msdata"), full.names = TRUE)
tmt_mzr <- backendInitialize(MsBackendMzR(), files = fl[5])

sciex_hd5 <- backendInitialize(MsBackendHdf5Peaks(),
                               spectraData = spectraData(sciex_mzr),
                               hdf5path = tempdir())

test_check("Spectra")
