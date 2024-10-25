# Low level infrastructure to handle MS spectra

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/Spectra/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/Spectra/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov](https://codecov.io/gh/rformassspectrometry/Spectra/branch/main/graph/badge.svg?token=jy0Mid9gKn)](https://codecov.io/gh/rformassspectrometry/Spectra)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)
[![years in bioc](http://bioconductor.org/shields/years-in-bioc/Spectra.svg)](https://bioconductor.org/packages/release/bioc/html/Spectra.html)
[![Ranking by downloads](http://bioconductor.org/shields/downloads/release/Spectra.svg)](https://bioconductor.org/packages/stats/bioc/Spectra/)
[![build release](http://bioconductor.org/shields/build/release/bioc/Spectra.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/Spectra/)
[![build devel](http://bioconductor.org/shields/build/devel/bioc/Spectra.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/Spectra/)

The Spectra package defines an efficient infrastructure for storing and handling
mass spectrometry spectra and functionality to subset, process, visualize and
compare spectra data. It provides different implementations (backends) to store
mass spectrometry data. These comprise backends tuned for fast data access and
processing and backends for very large data sets ensuring a small memory
footprint.

A (possibly incomplete) list of available backends (along with a link to the R
package providing it) is shown below:

- `MsBackendCompDb` (package
  [*CompoundDb*](https://github.com/rformassspectrometry/CompoundDb): provides
  access to spectra data (spectra and peaks variables) from a *CompDb*
  database. Has a small memory footprint because all data (except precursor m/z
  values) are retrieved on-the-fly from the database.

- `MsBackendDataFrame` (package: *Spectra*): alternative to the
  `MsBackendMemory` also keeping all data in memory, but supporting `S4` objects
  as spectra variables because the data is stored internally in a `DataFrame`.

- `MsBackendHdf5Peaks` (package: *Spectra*): on-disk backend similar to
  `MsBackendMzR`, but the peaks data is stored in HDF5 files (general spectra
  variables are kept in memory).

- `MsBackendHmdbXml` (package
  [*MsbackendHmdb*](https://github.com/rformassspectrometry/MsBackendHmdb)):
  allows import of MS data from xml files of the Human Metabolome Database
  (HMDB). Extends the `MsBackendDataFrame` and keeps thus all data, after
  import, in memory.

- `MsBackendMassbank` (package
  [*MsBackendMassbank*](https://github.com/rformassspectrometry/MsBackendMassbank)):
  allows to import/export data in MassBank text file format. Extends the
  `MsBackendDataFrame` and keeps thus all data, after import, in memory.

- `MsBackendMassbankSql` (package
  [*MsBackendMassbank*](https://github.com/rformassspectrometry/MsBackendMassbank)):
  allows to directly connect to a MassBank SQL database to retrieve all MS data
  and variables. Has a minimal memory footprint because all data is retrieved
  on-the-fly from the SQL database.

- `MsBackendMemory` (package: *Spectra*): *default* backend which keeps all data
  in memory. Optimized for fast processing.

- `MsBackendMetaboLights` (package
  [*MsBackendMetaboLights*](https://github.com/rformassspectrometry/MsBackendMetaboLights)):
  retrieves and caches MS data files from MetaboLights.

- `MsBackendMgf` (package
  [*MsBackendMgf*](https://github.com/rformassspectrometry/MsBackendMgf)): allows
  to import/export data in mascot generic format (MGF). Extends the
  `MsBackendDataFrame` and keeps thus all data, after import, in memory.

- `MsBackendMsp` (package
  [*MsbackendMsp*](https://github.com/rformassspectrometry/MsBackendMsp)): allows
  to import/export data in NIST MSP format. Extends the `MsBackendDataFrame` and
  keeps thus all data, after import, in memory.

- `MsBackendMzR` (package: *Spectra*): by using the `mzR` package it supports
  import of MS data from mzML, mzXML and CDF files. This backend keeps only
  general spectra variables in memory and retrieves the peaks data (m/z and
  intensity values) on-the-fly from the original data files. The backend has
  thus a smaller memory footprint compared to in-memory backends.

- `MsBackendOfflineSql` (package
  [*MsBackendSql*](https://github.com/rformassspectrometry/MsBackendSql)):
  stores all MS data in a SQL database and has thus a minimal memory footprint.
  Does, in contrast to `MsBackendSql`, not keep an active SQL database
  connection and can thus support parallel processing.

- `MsBackendRawFileReader` (package
  [*MsBackendRawFileReader*](https://github.com/fgcz/MsBackendRawFileReader)):
  implements a backend for reading MS data from Thermo Fisher Scientific's raw
  data files using the manufacturer's NewRawFileReader .Net libraries. The
  package generalizes the functionality introduced by the `rawrr` package.

- `MsBackendSql` (package
  [*MsBackendSql*](https://github.com/rformassspectrometry/MsBackendSql)):
  stores all MS data in a SQL database and has thus a minimal memory footprint.

- `MsBackendTimsTof` (package
  [*MsBackendTimsTof*](https://github.com/rformassspectrometry/MsBackendTimsTof):
  allows import of data from Bruker TimsTOF raw data files (using the
  `opentimsr` R package).

- `MsBackendWeizMass` (package
  [*MsBackendWeizMass*](https://github.com/rformassspectrometry/MsBackendWeizMass):
  allows to access MS data from WeizMass MS/MS spectral databases.


For more information see the package
[homepage](https://rformassspectrometry.github.io/Spectra).


# Installation

The package can be installed with

```r
install.packages("BiocManager")
BiocManager::install("Spectra")
```


# Contributions

Contributions are highly welcome and should follow the [contribution
guidelines](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#contributions).
Also, please check the coding style guidelines in the [RforMassSpectrometry
vignette](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html)
and importantly, follow our [code of
conduct](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#code-of-conduct).
