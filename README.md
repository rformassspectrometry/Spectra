# Low level infrastructure to handle MS spectra

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/Spectra/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/Spectra/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](http://codecov.io/github/rformassspectrometry/Spectra/coverage.svg?branch=master)](http://codecov.io/github/rformassspectrometry/Spectra?branch=master)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)

<img
src="https://raw.githubusercontent.com/rformassspectrometry/stickers/master/Spectra/Spectra.png"
height="150">

Externalising the new MS spectra backend-supporting infrastructure
from `MSnbase`. Could either survive as its own package or be
phagocyted back into by `MSnbase`.

# Basic ideas and concepts

- `Spectra` should become the main object to represent MS data.

- No explicit `Spectrum` objects, *spectrum* data (m/z - intensity
  pairs) will be represented as a `matrix` stored/provided by the
  `Backend` class.
