# Spectra 0.99

## Changes in 0.99.1

- Change parameter `source` in `Spectra,character` to `MsBackendMzR` and set
  parameter `backend = source`. Thus by default, the import backend will also
  be used to store the data.

## Changes in 0.99.0

- Add `reset` method.
- Add processing by chunk to `compareSpectra`.
