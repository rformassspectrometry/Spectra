library("testthat")
library("patrick")
library("Spectra")
library("MsDataHub")

register(SerialParam())

## MsDataHub
cdf_file <- unname(MsDataHub::ko15.CDF())

sciex_file <- unname(c(MsDataHub::X20171016_POOL_POS_1_105.134.mzML(),
                       MsDataHub::X20171016_POOL_POS_3_105.134.mzML()))
sciex_mzr <- backendInitialize(MsBackendMzR(), files = sciex_file)
sciex_pks <- peaksData(sciex_mzr)

fl <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzML.gz()
tmt_mzr <- backendInitialize(MsBackendMzR(), files = fl)
sps_tmt <- setBackend(Spectra(tmt_mzr), MsBackendMemory())

fl <- MsDataHub::PestMix1_SWATH.mzML()
sps_dia <- Spectra(fl)

fl <- MsDataHub::PestMix1_DDA.mzML()
sps_dda <- Spectra(fl)

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
