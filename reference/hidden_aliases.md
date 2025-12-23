# Internal page for hidden aliases

For S4 methods that require a documentation entry but only clutter the
index.

## Usage

``` r
# S4 method for class 'numeric'
bin(
  x,
  y,
  size = 1,
  breaks = seq(floor(min(y)), ceiling(max(y)), by = size),
  FUN = max,
  returnMids = TRUE,
  .check = TRUE
)

# S4 method for class 'MsBackendCached'
lengths(x, use.names = FALSE)

# S4 method for class 'MsBackendDataFrame'
show(object)

# S4 method for class 'MsBackendDataFrame'
backendMerge(object, ...)

# S4 method for class 'MsBackendDataFrame'
backendRequiredSpectraVariables(object, ...)

# S4 method for class 'MsBackendDataFrame'
acquisitionNum(object)

# S4 method for class 'MsBackendDataFrame'
peaksData(object, columns = c("mz", "intensity"))

# S4 method for class 'MsBackendDataFrame'
centroided(object)

# S4 method for class 'MsBackendDataFrame'
centroided(object) <- value

# S4 method for class 'MsBackendDataFrame'
collisionEnergy(object)

# S4 method for class 'MsBackendDataFrame'
collisionEnergy(object) <- value

# S4 method for class 'MsBackendDataFrame'
dataOrigin(object)

# S4 method for class 'MsBackendDataFrame'
dataOrigin(object) <- value

# S4 method for class 'MsBackendDataFrame'
dataStorage(object)

# S4 method for class 'MsBackendDataFrame'
dataStorage(object) <- value

# S4 method for class 'MsBackendDataFrame,ANY'
extractByIndex(object, i)

# S4 method for class 'MsBackendDataFrame'
intensity(object)

# S4 method for class 'MsBackendDataFrame'
intensity(object) <- value

# S4 method for class 'MsBackendDataFrame'
isEmpty(x)

# S4 method for class 'MsBackendDataFrame'
isolationWindowLowerMz(object)

# S4 method for class 'MsBackendDataFrame'
isolationWindowLowerMz(object) <- value

# S4 method for class 'MsBackendDataFrame'
isolationWindowTargetMz(object)

# S4 method for class 'MsBackendDataFrame'
isolationWindowTargetMz(object) <- value

# S4 method for class 'MsBackendDataFrame'
isolationWindowUpperMz(object)

# S4 method for class 'MsBackendDataFrame'
isolationWindowUpperMz(object) <- value

# S4 method for class 'MsBackendDataFrame'
length(x)

# S4 method for class 'MsBackendDataFrame'
lengths(x, use.names = FALSE)

# S4 method for class 'MsBackendDataFrame'
msLevel(object, ...)

# S4 method for class 'MsBackendDataFrame'
msLevel(object) <- value

# S4 method for class 'MsBackendDataFrame'
mz(object)

# S4 method for class 'MsBackendDataFrame'
mz(object) <- value

# S4 method for class 'MsBackendDataFrame'
polarity(object)

# S4 method for class 'MsBackendDataFrame'
polarity(object) <- value

# S4 method for class 'MsBackendDataFrame'
precScanNum(object)

# S4 method for class 'MsBackendDataFrame'
precursorCharge(object)

# S4 method for class 'MsBackendDataFrame'
precursorIntensity(object)

# S4 method for class 'MsBackendDataFrame'
precursorMz(object)

# S4 method for class 'MsBackendDataFrame'
peaksData(object) <- value

# S4 method for class 'MsBackendDataFrame'
peaksVariables(object)

# S4 method for class 'MsBackendDataFrame'
rtime(object)

# S4 method for class 'MsBackendDataFrame'
rtime(object) <- value

# S4 method for class 'MsBackendDataFrame'
scanIndex(object)

# S4 method for class 'MsBackendDataFrame'
selectSpectraVariables(object, spectraVariables = spectraVariables(object))

# S4 method for class 'MsBackendDataFrame'
smoothed(object)

# S4 method for class 'MsBackendDataFrame'
smoothed(object) <- value

# S4 method for class 'MsBackendDataFrame'
spectraData(object, columns = spectraVariables(object))

# S4 method for class 'MsBackendDataFrame'
spectraData(object) <- value

# S4 method for class 'MsBackendDataFrame'
spectraNames(object)

# S4 method for class 'MsBackendDataFrame'
spectraNames(object) <- value

# S4 method for class 'MsBackendDataFrame'
spectraVariables(object)

# S4 method for class 'MsBackendDataFrame'
tic(object, initial = TRUE)

# S4 method for class 'MsBackendDataFrame'
x$name

# S4 method for class 'MsBackendDataFrame'
x$name <- value

# S4 method for class 'MsBackendDataFrame'
x[i, j, ..., drop = FALSE]

# S4 method for class 'MsBackendDataFrame,dataframeOrDataFrameOrmatrix'
cbind2(x, y = data.frame(), ...)

# S4 method for class 'MsBackendDataFrame,ANY'
split(x, f, drop = FALSE, ...)

# S4 method for class 'MsBackendDataFrame'
filterAcquisitionNum(
  object,
  n = integer(),
  dataStorage = character(),
  dataOrigin = character()
)

# S4 method for class 'MsBackendHdf5Peaks'
backendRequiredSpectraVariables(object, ...)

# S4 method for class 'MsBackendHdf5Peaks'
backendInitialize(
  object,
  files = character(),
  data = DataFrame(),
  hdf5path = character(),
  ...,
  BPPARAM = bpparam()
)

# S4 method for class 'MsBackendHdf5Peaks'
show(object)

# S4 method for class 'MsBackendHdf5Peaks'
peaksData(object, columns = peaksVariables(object))

# S4 method for class 'MsBackendHdf5Peaks'
intensity(object)

# S4 method for class 'MsBackendHdf5Peaks'
intensity(object) <- value

# S4 method for class 'MsBackendHdf5Peaks'
ionCount(object)

# S4 method for class 'MsBackendHdf5Peaks'
isCentroided(object, ...)

# S4 method for class 'MsBackendHdf5Peaks'
mz(object)

# S4 method for class 'MsBackendHdf5Peaks'
mz(object) <- value

# S4 method for class 'MsBackendHdf5Peaks'
peaksData(object) <- value

# S4 method for class 'MsBackendHdf5Peaks'
spectraData(object, columns = spectraVariables(object))

# S4 method for class 'MsBackendHdf5Peaks'
spectraData(object) <- value

# S4 method for class 'MsBackendHdf5Peaks'
x$name <- value

# S4 method for class 'MsBackendHdf5Peaks'
x[i, j, ..., drop = FALSE]

# S4 method for class 'MsBackendHdf5Peaks,ANY'
extractByIndex(object, i)

# S4 method for class 'MsBackendHdf5Peaks'
backendMerge(object, ...)

# S4 method for class 'MsBackendMemory'
show(object)

# S4 method for class 'MsBackendMemory'
backendMerge(object, ...)

# S4 method for class 'MsBackendMemory'
backendRequiredSpectraVariables(object, ...)

# S4 method for class 'MsBackendMemory'
acquisitionNum(object)

# S4 method for class 'MsBackendMemory'
centroided(object)

# S4 method for class 'MsBackendMemory'
centroided(object) <- value

# S4 method for class 'MsBackendMemory'
collisionEnergy(object)

# S4 method for class 'MsBackendMemory'
collisionEnergy(object) <- value

# S4 method for class 'MsBackendMemory'
dataOrigin(object)

# S4 method for class 'MsBackendMemory'
dataOrigin(object) <- value

# S4 method for class 'MsBackendMemory'
dataStorage(object)

# S4 method for class 'MsBackendMemory'
dataStorage(object) <- value

# S4 method for class 'MsBackendMemory,ANY'
extractByIndex(object, i)

# S4 method for class 'MsBackendMemory'
intensity(object)

# S4 method for class 'MsBackendMemory'
intensity(object) <- value

# S4 method for class 'MsBackendMemory'
ionCount(object)

# S4 method for class 'MsBackendMemory'
isolationWindowLowerMz(object)

# S4 method for class 'MsBackendMemory'
isolationWindowLowerMz(object) <- value

# S4 method for class 'MsBackendMemory'
isolationWindowTargetMz(object)

# S4 method for class 'MsBackendMemory'
isolationWindowTargetMz(object) <- value

# S4 method for class 'MsBackendMemory'
isolationWindowUpperMz(object)

# S4 method for class 'MsBackendMemory'
isolationWindowUpperMz(object) <- value

# S4 method for class 'MsBackendMemory'
length(x)

# S4 method for class 'MsBackendMemory'
msLevel(object, ...)

# S4 method for class 'MsBackendMemory'
msLevel(object) <- value

# S4 method for class 'MsBackendMemory'
mz(object)

# S4 method for class 'MsBackendMemory'
mz(object) <- value

# S4 method for class 'MsBackendMemory'
peaksData(object, columns = c("mz", "intensity"))

# S4 method for class 'MsBackendMemory'
peaksData(object) <- value

# S4 method for class 'MsBackendMemory'
polarity(object)

# S4 method for class 'MsBackendMemory'
polarity(object) <- value

# S4 method for class 'MsBackendMemory'
precScanNum(object)

# S4 method for class 'MsBackendMemory'
precursorCharge(object)

# S4 method for class 'MsBackendMemory'
precursorIntensity(object)

# S4 method for class 'MsBackendMemory'
precursorMz(object)

# S4 method for class 'MsBackendMemory'
rtime(object)

# S4 method for class 'MsBackendMemory'
rtime(object) <- value

# S4 method for class 'MsBackendMemory'
scanIndex(object)

# S4 method for class 'MsBackendMemory'
selectSpectraVariables(object, spectraVariables = spectraVariables(object))

# S4 method for class 'MsBackendMemory'
smoothed(object)

# S4 method for class 'MsBackendMemory'
smoothed(object) <- value

# S4 method for class 'MsBackendMemory'
spectraData(object, columns = spectraVariables(object))

# S4 method for class 'MsBackendMemory'
spectraData(object) <- value

# S4 method for class 'MsBackendMemory'
spectraNames(object)

# S4 method for class 'MsBackendMemory'
spectraNames(object) <- value

# S4 method for class 'MsBackendMemory'
spectraVariables(object)

# S4 method for class 'MsBackendMemory'
peaksVariables(object)

# S4 method for class 'MsBackendMemory'
tic(object, initial = TRUE)

# S4 method for class 'MsBackendMemory'
x$name

# S4 method for class 'MsBackendMemory'
x$name <- value

# S4 method for class 'MsBackendMemory'
x[i, j, ..., drop = FALSE]

# S4 method for class 'MsBackendMemory,dataframeOrDataFrameOrmatrix'
cbind2(x, y = data.frame(), ...)

# S4 method for class 'MsBackendMemory,ANY'
split(x, f, drop = FALSE, ...)

# S4 method for class 'MsBackendMemory'
filterAcquisitionNum(
  object,
  n = integer(),
  dataStorage = character(),
  dataOrigin = character()
)

# S4 method for class 'MsBackendMemory'
longForm(object, columns = spectraVariables(object))

# S4 method for class 'MsBackendMzR'
backendRequiredSpectraVariables(object, ...)

# S4 method for class 'MsBackendMzR'
backendInitialize(object, files, ..., BPPARAM = bpparam())

# S4 method for class 'MsBackendMzR'
show(object)

# S4 method for class 'MsBackendMzR'
peaksData(object, columns = peaksVariables(object))

# S4 method for class 'MsBackendMzR'
intensity(object)

# S4 method for class 'MsBackendMzR'
intensity(object) <- value

# S4 method for class 'MsBackendMzR'
ionCount(object)

# S4 method for class 'MsBackendMzR'
isCentroided(object, ...)

# S4 method for class 'MsBackendMzR'
mz(object)

# S4 method for class 'MsBackendMzR'
mz(object) <- value

# S4 method for class 'MsBackendMzR'
spectraData(object, columns = spectraVariables(object))

# S4 method for class 'MsBackendMzR'
spectraData(object) <- value

# S4 method for class 'MsBackendMzR'
spectraNames(object)

# S4 method for class 'MsBackendMzR'
spectraNames(object) <- value

# S4 method for class 'MsBackendMzR'
spectraVariables(object)

# S4 method for class 'MsBackendMzR'
x$name <- value

# S4 method for class 'MsBackendMzR'
export(
  object,
  x,
  file = tempfile(),
  format = c("mzML", "mzXML"),
  copy = FALSE,
  BPPARAM = bpparam()
)

# S4 method for class 'Spectra'
show(object)

# S4 method for class 'list'
combinePeaks(object, ...)
```

## Value

Not applicable

## Note

: this replaces **all** the data in the backend.
