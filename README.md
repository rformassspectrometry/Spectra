# Low level infrastructure to handle MS spectra

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![build status](https://travis-ci.org/rformassspectrometry/Spectra.svg?branch=master)](https://travis-ci.org/rformassspectrometry/Spectra)
[![codecov.io](http://codecov.io/github/rformassspectrometry/Spectra/coverage.svg?branch=master)](http://codecov.io/github/rformassspectrometry/Spectra?branch=master)
[![license](http://img.shields.io/badge/license-GPL%20%28%3E=%203%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)

Externalising the new MS spectra backend-supporting infrastructure
from `MSnbase`. Could either survive as its own package or be
phagocyted back into by `MSnbase`.

# Basic ideas and concepts

- `Spectra` should become the main object to represent MS data.

- No explicit `Spectrum` objects, *spectrum* data (m/z - intensity
  pairs) will be represented as a `matrix` stored/provided by the
  `Backend` class.
