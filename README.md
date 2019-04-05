# Low level infrastructure to handle MS spectra

Externalising the new MS spectra backend-supporting infrastructure from
`MSnbase`. Could either survive as its own package or be phagocyted back into by
`MSnbase`.

# Basic ideas and concepts

- `Spectra` should become the main object to represent MS data.
- No `Spectrum` objects, *spectrum* data (m/z - intensity pairs) will be
  represented as a `data.frame` stored/provided by the `Backend` class.
