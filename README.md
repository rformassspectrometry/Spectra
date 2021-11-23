# Low level infrastructure to handle MS spectra

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/Spectra/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/Spectra/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](http://codecov.io/github/rformassspectrometry/Spectra/coverage.svg?branch=master)](http://codecov.io/github/rformassspectrometry/Spectra?branch=master)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)

<img
src="https://raw.githubusercontent.com/rformassspectrometry/stickers/master/Spectra/Spectra.png"
height="150">

The Spectra package defines an efficient infrastructure
for storing and handling mass spectrometry spectra and functionality to
subset, process, visualize and compare spectra data. It provides different
implementations (backends) to store mass spectrometry data. These comprise
backends tuned for fast data access and processing and backends for very
large data sets ensuring a small memory footprint.
