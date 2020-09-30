# Spectra 0.99

## Changes in 0.99.7

- Fix problem in `dropNaSpectraVariables` that would also drop m/z and
  intensity values for most backends (issue
  [#138](https://github.com/rformassspectrometry/Spectra/issues/138).

## Changes in 0.99.6

- Support `intensity` in `filterIntensity` method to be a function to enable
  peak intensity-based filtering of spectra (issue
  [#126](https://github.com/rformassspectrometry/Spectra/issues/126)).

## Changes in 0.99.5

- Add `filterMzRange` and `filterMzValues` to filter spectra based on an m/z
  range or a list of target m/z values, respectively.

## Changes in 0.99.4

- Add `export,MsBackendMzR` to export spectra data to mzML or mzXML file(s).
- Add an `export,MsBackend` method to allow backends to take care of data
  export.
- Refactor `export,Spectra` to use the `MsBackend` class to export the data.
- Change parameter `source` in `Spectra,character` to `MsBackendMzR` and set
  parameter `backend = source`. Thus by default, the import backend will also
  be used to store the data.

## Changes in 0.99.3

- Replace `lapply,Spectra` with `spectrapply,Spectra`.

## Changes in 0.99.2

- Replace `asDataFrame,MsBackend` with `spectraData,MsBackend`.
- Replace `asDataFrame<-,MsBackend` with `spectraData<-,MsBackend`.
- Replace `as.list,MsBackend` with `peaksData,MsBackend`.
- Replace `replaceList<-,MsBackend` with `peaksData<-,MsBackend`.
- Replace `as.list,Spectra` with `peaksData,Spectra` and add methods to coerce a
  `Spectra` to a `list` or `SimpleList`.

## Changes in 0.99.0

- Add `reset` method.
- Add processing by chunk to `compareSpectra`.
