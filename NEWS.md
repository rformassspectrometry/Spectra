# Spectra 1.1

## Changes in 1.1.3

- Implement a generic `Spectra,ANY` constructor replacing `Spectra,DataFrame`
  and `Spectra,character`.

## Changes in 1.1.2

- Fix problem in export to mzML files that failed for empty spectra (issue
  [#145](https://github.com/rformassspectrometry/Spectra/issues/145))

## Changes in 1.1.1

- Round retention time in figure titles.
- Document differences between `spectrumId` (`spectrumID`),
  `acquisitionNum` and `scanIndex`.

## Changes in 1.1.0

- New Bioc devel version

# Spectra 0.99

## Changes in 0.99.11

- Re-add `mz` and `intensity` as core spectra variables.

## Changes in 0.99.10

- Fix in `spectraData<-,Spectra` to avoid removing m/z and intensity values
  (issue [#146](https://github.com/rformassspectrometry/Spectra/issues/146)).
- Add default implementations of filter functions for `MsBackend`.

## Changes in 0.99.9

- Fix in `Spectra,character` constructor to ensure the backend is changed even
  if `source` inherits from `backend` (issue
  [#143](https://github.com/rformassspectrometry/Spectra/issues/143)).

## Changes in 0.99.8

- `combineSpectra` applies data processing steps in the processing queue prior to
  combination (issue
  [#140](https://github.com/rformassspectrometry/Spectra/issues/140)).

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
