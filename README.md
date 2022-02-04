# Low level infrastructure to handle MS spectra

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/Spectra/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/Spectra/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](http://codecov.io/github/rformassspectrometry/Spectra/coverage.svg?branch=master)](http://codecov.io/github/rformassspectrometry/Spectra?branch=master)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)
[![years in bioc](http://bioconductor.org/shields/years-in-bioc/Spectra.svg)](https://bioconductor.org/packages/release/bioc/html/Spectra.html)
[![Ranking by downloads](http://bioconductor.org/shields/downloads/release/Spectra.svg)](https://bioconductor.org/packages/stats/bioc/Spectra/)
Release: [![build release](http://bioconductor.org/shields/build/release/bioc/Spectra.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/Spectra/)
Devel: [![build devel](http://bioconductor.org/shields/build/devel/bioc/Spectra.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/Spectra/)

The Spectra package defines an efficient infrastructure
for storing and handling mass spectrometry spectra and functionality to
subset, process, visualize and compare spectra data. It provides different
implementations (backends) to store mass spectrometry data. These comprise
backends tuned for fast data access and processing and backends for very
large data sets ensuring a small memory footprint.
