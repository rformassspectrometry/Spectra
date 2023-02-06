# Spectra 1.9

## Changes in 1.9.8

- Address comments from Michele Stravs regarding the `MsBackend` vignette.
- Add additional tests checking for `MsBackend` compliance.

## Changes in 1.9.7

- Add a vignette describing how to build a `MsBackend` from scratch (issue
  #262).
- Extend unit test suite to evaluate validity of `MsBackend` implementations.

## Changes in 1.9.6

- Replace `<=` with `between` calls.

## Changes in 1.9.5

- Fix bug in `containsMz()` when `mz` isn't ordered (see #258).

## Changes in 1.9.4

- Fix error when extracting spectra variables from a `MsBackendMzR` of
  length 0.

## Changes in 1.9.3

- Add `chunkapply` function to split a `Spectra` into chunks and
  stepwise apply a function `FUN` to each.

## Changes in 1.9.2

- `combineSpectra` on `Spectra` with read-only backends change backend
  to an `MsBackendMemory` instead of an `MsBackendDataFrame`.

## Changes in 1.9.1

- Expand documentation on `compareSpectra` for GNPS-like similarity
  scoring.

## Changes in 1.9.0

- Bioconductor 3.17 developmental version.

# Spectra 1.7

## Changes in 1.7.5

- Force serial processing in some unit tests to avoid potential
  failures on some Bioconductor build and check servers (under some
  circumstances).

## Changes in 1.7.4

- Add `MsBackendMemory` backend class providing a more efficient
  in-memory data representation than `MsBackendDataFrame`.

## Changes in 1.7.3

- Import `spectrapply` from `ProtGenerics`.

## Changes in 1.7.2

- Fix `setBackend` if provided `Spectra` is empty.
- `backendInitialize,Spectra,MsBackendDataFrame` returns a `Spectra`
  object with the full provided spectra data.

## Changes in 1.7.1

- Add `uniqueMsLevels` function to allow more efficient,
  backend-specific, implementations for retrieving unique MS levels
  from a data set.

# Spectra 1.5

## Changes in 1.5.20

- Add parameters `ppm` and `tolerance` to `PrecursorMzParam` (for
  neutral loss calculation) and add option `filterPeaks =
  "removePrecursor"`.

## Changes in 1.5.19

- Improved the `bin` method.

## Changes in 1.5.18

- Set default for parameter `columns` in `peaksData,Spectra` and
  `peaksData,MsBackend` to `c("mz", "intensity")`.

## Changes in 1.5.17

- Add `peaksVariables` method and add parameter `columns` (or `...`) to
  `peaksData`.
- Add `columns` parameter to the `peaksData` method of
  `MsBackendDataFrame`, `MsBackendMzR` and `MsBackendHdf5peaks`.

## Changes in 1.5.16

- Fix issue in `neutralLoss` that would prevent calculation of neutral
  loss spectra if

## Changes in 1.5.15

- Fix typo in MZ delta plot title.

## Changes in 1.5.14

- Add `coreSpectraVariables` function to export the *core* spectra
  variables and their expected data types.

## Changes in 1.5.13

- Fix figure sizes in vignette.

## Changes in 1.5.12

- Add `neutralLoss` method and first algorithm to calculate neutral loss
  spectra.

## Changes in 1.5.11

- Fix neutral loss example in the vignette.

## Changes in 1.5.10

- Add citation.

## Changes in 1.5.9

- Add examples for `combineSpectra` to the vignette.

## Changes in 1.5.8

- Add `spectraVariableMapping` generic.

## Changes in 1.5.7

- Add missing export of the `filterPrecursorMz` method.

## Changes in 1.5.6

- Add `filterPrecursorMzValue` method which allows to filter using
  multiple precursor m/z values (issue
  [#230](https://github.com/rformassspectrometry/Spectra/issues/230)).
- Fix unit test suite.

## Changes in 1.5.5

- Add a testing framework allowing to run standardized unit tests for
  new `MsBackend` implementations (issue
  [#186](https://github.com/rformassspectrometry/Spectra/issues/186)).

## Changes in 1.5.4

- Add the `MsBackendCached` backend.

## Changes in 1.5.3

- Only calculate number of peaks per spectra if the processing queue
  of the `Spectra` is not empty. Otherwise call the backend's
  implementation (issue [MsBackendSql
  #31](https://github.com/rformassspectrometry/MsBackendSql/issues/31)).

## Changes in 1.5.2

- Small documentation update (related to `MsCoreUtils` issue
  [#87](https://github.com/rformassspectrometry/MsCoreUtils/issues/87)).
- New `countIdentifications()` function.
- Add `filterFourierTransformArtefacts` function to remove fast
  fourier artefact peaks seen on e.g. Orbitrap instruments (issue
  #223).

## Changes in 1.5.1

- Don't read header information when importing peaks matrix on macOS.

# Spectra 1.3

## Changes in 1.3.11

- Fix error message in `setBackend` (issue
  [#217](https://github.com/rformassspectrometry/Spectra/issues/217)).

## Changes in 1.3.10

- Fix bug in `plotSpectra` and `plotSpectraMirror` that would cause an error if
  the number of peaks in a spectrum was 1 and labels were provided.

## Changes in 1.3.9

- New features: `joinSpectraData()` now check for duplicated keys in
  `x` (throws an error) and `y` (thows a warning).

## Changes in 1.3.8

- New features: `plotMzDelta()` function to M/Z delta QC (ported from
  MSnbase).

## Changes in 1.3.7

- Add fix from `MSnbase` (issue
  [#170](https://github.com/lgatto/MSnbase/issues/170)) to `Spectra`: on macOS
  require reading also the spectrum header before reading the peaks data.

## Changes in 1.3.6

- Documentation updates for `combineSpectra` and `combinePeaks`.

## Changes in 1.3.5

- `filterMzValues` supports also removing peaks matching specified m/z values
  (issue [#209](https://github.com/rformassspectrometry/Spectra/issues/209)).

## Changes in 1.3.4

- Add list of additional R packages and repositories providing `MsBackend`
  backends to the vignette.

## Changes in 1.3.3

- Move generics for `bin` and `compareSpectra` to `ProtGenerics`.

## Changes in 1.3.2

- Add parameter `f` to `filterPrecursorScan` to fix issue
  [#194](https://github.com/rformassspectrometry/Spectra/issues/194).

## Changes in 1.3.1

- Add `estimatePrecursorIntensity` function (issue
  [#202](https://github.com/rformassspectrometry/Spectra/issues/202)).


# Spectra 1.1

## Changes in 1.1.20

- Fix concatenating empty spectra (issue
  [#200](https://github.com/rformassspectrometry/Spectra/issues/200)).

## Changes in 1.1.19

- New `filterPrecursorCharge()` method.

## Changes in 1.1.18

- Define `plotSpectraMirror` as a method.

## Changes in 1.1.17

- Fix issue [#187](https://github.com/rformassspectrometry/Spectra/issues/187).
- Add function `concatenateSpectra` to allow concatenating `Spectra` objects and
  list of `Spectra` objects.

## Changes in 1.1.16

- Support arbitrary spectra variables to be passed to the functions
  provided/added with `addProcessing`; issue
  [#182](https://github.com/rformassspectrometry/Spectra/issues/182).

## Changes in 1.1.15

- Pass spectras' precursor m/z to the `MAPFUN` in `compareSpectra`; issue
  [#171](https://github.com/rformassspectrometry/Spectra/issues/171).
- Add `joinPeaksGnps` to perform a peak matching between spectra similar to the
  one performed in GNPS (issue
  [#171](https://github.com/rformassspectrometry/Spectra/issues/171)).

## Changes in 1.1.14

- Support plotting of empty spectra (issue
  [175](https://github.com/rformassspectrometry/Spectra/issues/175)).

## Changes in 1.1.13

- Move `ProcessingStep` to `ProtGenerics`.

## Changes in 1.1.12

- Fix `show` method for `Spectra` to list only the 3 most recent processing
  steps (issue
  [173](https://github.com/rformassspectrometry/Spectra/issues/173)).
- Add `processingLog` function to display the log messages of all processing
  steps of a `Spectra` object.

## Changes in 1.1.11

- Add support for `...` to `pickPeaks` and `smooth` (issue
  [168](https://github.com/rformassspectrometry/Spectra/issues/168)).

## Changes in 1.1.10

- Import `filterIntensity` from `ProtGenerics`.

## Changes in 1.1.9

- Fix label in `plotSpectra`.

## Changes in 1.1.8

- `filterIntensity` supports passing of additional parameters to the used
  filter function (issue
  [164](https://github.com/rformassspectrometry/Spectra/issues/164)).

## Changes in 1.1.7

- Fix bug in `show,ProcessingStep` (issue
  [162](https://github.com/rformassspectrometry/Spectra/issues/162)).

## Changes in 1.1.6

- New `joinSpectraData()` function.

## Changes in 1.1.5

- Add `[[,Msbackend` and `[[<-,MsBackend` methods (issue
  [149](https://github.com/rformassspectrometry/Spectra/issues/149)).
- Add `[[,Spectra` and `[[<-,Spectra` methods.

## Changes in 1.1.4

- Fix issue with `labelCol` in `plotSpectra` (issue
  [#157](https://github.com/rformassspectrometry/Spectra/issues/157)).

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
