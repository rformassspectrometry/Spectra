# Spectra 0.99

## Changes in 0.99.3

- Add `export,MsBackendMzR` to export spectra data to mzML or mzXML file(s).
- Add an `export,MsBackend` method to allow backends to take care of data
  export.
- Refactor `export,Spectra` to use the `MsBackend` class to export the data.
- Change parameter `source` in `Spectra,character` to `MsBackendMzR` and set
  parameter `backend = source`. Thus by default, the import backend will also
  be used to store the data.

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
