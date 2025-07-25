% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spectra.R
\name{Spectra}
\alias{Spectra}
\alias{Spectra-class}
\alias{setBackend}
\alias{export}
\alias{Spectra,missing-method}
\alias{Spectra,MsBackend-method}
\alias{Spectra,character-method}
\alias{Spectra,ANY-method}
\alias{setBackend,Spectra,MsBackend-method}
\alias{export,Spectra-method}
\alias{dataStorageBasePath,Spectra-method}
\alias{dataStorageBasePath<-,Spectra-method}
\title{The Spectra class to manage and access MS data}
\usage{
\S4method{Spectra}{missing}(
  object,
  processingQueue = list(),
  metadata = list(),
  ...,
  backend = MsBackendMemory(),
  BPPARAM = bpparam()
)

\S4method{Spectra}{MsBackend}(object, processingQueue = list(), metadata = list(), ...)

\S4method{Spectra}{character}(
  object,
  processingQueue = list(),
  metadata = list(),
  source = MsBackendMzR(),
  backend = source,
  ...,
  BPPARAM = bpparam()
)

\S4method{Spectra}{ANY}(
  object,
  processingQueue = list(),
  metadata = list(),
  source = MsBackendMemory(),
  backend = source,
  ...,
  BPPARAM = bpparam()
)

\S4method{setBackend}{Spectra,MsBackend}(
  object,
  backend,
  f = processingChunkFactor(object),
  ...,
  BPPARAM = bpparam()
)

\S4method{export}{Spectra}(object, backend, ...)

\S4method{dataStorageBasePath}{Spectra}(object)

\S4method{dataStorageBasePath}{Spectra}(object) <- value
}
\arguments{
\item{object}{For \code{Spectra()}: an object to instantiate the \code{Spectra}
object and initialize the with data.. See section on creation of
\code{Spectra} objects for details. For all other methods a \code{Spectra} object.}

\item{processingQueue}{For \code{Spectra()}: optional \code{list} of
\link[ProtGenerics:ProcessingStep]{ProtGenerics::ProcessingStep} objects.}

\item{metadata}{For \code{Spectra()}: optional \code{list} with metadata information.}

\item{...}{Additional arguments.}

\item{backend}{For \code{Spectra()}: \linkS4class{MsBackend} to be used as backend. See
section on creation of \code{Spectra} objects for details. For \code{setBackend()}:
instance of \linkS4class{MsBackend} that supports \code{setBackend()} (i.e. for
which \code{supportsSetBackend()} returns \code{TRUE}). Such backends have a
parameter \code{data} in their \code{backendInitialize()} function that support
passing the full spectra data to the initialize method. See section on
creation of \code{Spectra} objects for details.
For \code{export()}: \linkS4class{MsBackend} to be used to export the data.}

\item{BPPARAM}{Parallel setup configuration. See \code{\link[BiocParallel:register]{BiocParallel::bpparam()}}
for more information. This is passed directly to the
\code{\link[=backendInitialize]{backendInitialize()}} method of the \linkS4class{MsBackend}.}

\item{source}{For \code{Spectra()}: instance of \linkS4class{MsBackend} that can be
used to import spectrum data from the provided files. See section
\emph{Creation of objects} for more details.}

\item{f}{For \code{setBackend()}: factor defining how to split the data
for parallelized copying of the spectra data to the new backend. For
some backends changing this parameter can lead to errors. Defaults to
\code{\link[=processingChunkFactor]{processingChunkFactor()}}.}

\item{value}{For \code{dataStorageBasePath()}: A \code{character} vector that defines
the base directory where the data storage files can be found.}
}
\description{
The \code{Spectra} class encapsules spectral mass spectrometry (MS) data and
related metadata. The MS data is represented by a \emph{backend} extending the
virual \link{MsBackend} class which provides the data to the \code{Spectra} object.
The \code{Spectra} class implements only data accessor, filtering and analysis
methods for the MS data and relies on its \emph{backend} to provide the MS data.
This allows to change data representations of a \code{Spectra} object depending
on the user's needs and properties of the data. Different backends and
their properties are explained in the \link{MsBackend} documentation.

Documentation on other topics and functionality of \code{Spectra}can be found in:
\itemize{
\item \code{\link[=spectraData]{spectraData()}} for accessing and using MS data through \code{Spectra} objects.
\item \code{\link[=filterMsLevel]{filterMsLevel()}} to subset and filter \code{Spectra} objects.
\item \code{\link[=plotSpectra]{plotSpectra()}} for visualization of \code{Spectra} objects.
\item \code{\link[=processingChunkSize]{processingChunkSize()}} for information on parallel and chunk-wise data
processing.
\item \code{\link[=combineSpectra]{combineSpectra()}} for merging, aggregating and splitting of \code{Spectra}
objects.
\item \code{\link[=combinePeaks]{combinePeaks()}} for merging and aggregating \code{Spectra}'s mass peaks data.
\item \code{\link[=addProcessing]{addProcessing()}} for data analysis functions.
\item \code{\link[=compareSpectra]{compareSpectra()}} for spectra similarity calculations.
}
}
\details{
The \code{Spectra} class uses by default a lazy data manipulation strategy,
i.e. data manipulations such as performed with \code{replaceIntensitiesBelow()}
are not applied immediately to the data, but applied on-the-fly to the
spectrum data once it is retrieved. This enables data manipulation
operations also for \emph{read only} data representations. For some backends that
allow to write data back to the data storage (such as the
\code{\link[=MsBackendMemory]{MsBackendMemory()}}, \code{\link[=MsBackendDataFrame]{MsBackendDataFrame()}} and \code{\link[=MsBackendHdf5Peaks]{MsBackendHdf5Peaks()}}) it
is possible to apply to queue with the \code{\link[=applyProcessing]{applyProcessing()}} function (see
the \code{\link[=applyProcessing]{applyProcessing()}} function for details).

Clarifications regarding scan/acquisition numbers and indices:
\itemize{
\item A \code{spectrumId} (or \code{spectrumID}) is a vendor specific field in
the mzML file that contains some information about the
run/spectrum, e.g.: \verb{controllerType=0 controllerNumber=1 scan=5281 file=2}
\item \code{acquisitionNum} is a more a less sanitize spectrum id generated
from the \code{spectrumId} field by \code{mzR} (see
\href{https://github.com/sneumann/mzR/blob/master/src/pwiz/data/msdata/MSData.cpp#L552-L580}{here}).
\item \code{scanIndex} is the \code{mzR} generated sequence number of the
spectrum in the raw file (which doesn't have to be the same as
the \code{acquisitionNum})
}

See also \href{https://github.com/lgatto/MSnbase/issues/525}{this issue}.
}
\section{Data stored in a \code{Spectra} object}{


The \code{Spectra} object is a container for MS data that includes mass peak
data (\emph{m/z} and related intensity values, also referred to as \emph{peaks data}
in the context of \code{Spectra}) and metadata of individual spectra (so called
\emph{spectra variables}). While a core set of spectra variables (the
\code{coreSpectraVariables()}) are guaranteed to be provided by a
\code{Spectra}, it is possible to add arbitrary additional spectra variables to
a \code{Spectra} object.

The \code{Spectra} object is designed to contain MS data of a (large) set of mass
spectra. The data is organized \emph{linearly} and can be thought of a list of
mass spectra, i.e. each element in the \code{Spectra} is one spectrum.
}

\section{Creation of objects}{


\code{Spectra} classes can be created with the \code{Spectra()} constructor function
which supports the following formats:
\itemize{
\item parameter \code{object} is a \code{data.frame} or \code{DataFrame} containing the
full spectrum data (spectra variables in columns as well as columns
with the individual MS peak data, \emph{m/z} and intensity). The provided
\code{backend} (by default a \linkS4class{MsBackendMemory}) will be initialized
with that data.
\item parameter \code{object} is a \linkS4class{MsBackend} (assumed to be already
initialized).
\item parameter \code{object} is missing, in which case it is supposed that the data
is provided by the \linkS4class{MsBackend} class passed along with the \code{backend}
argument.
\item parameter \code{object} is of type \code{character} and is expected to be the file
names(s) from which spectra should be imported. Parameter \code{source} allows
to define a \linkS4class{MsBackend} that is able to import the data from the
provided source files. The default value for \code{source} is \code{\link[=MsBackendMzR]{MsBackendMzR()}}
which allows to import spectra data from mzML, mzXML or CDF files.
}

With \code{...} additional arguments can be passed to the backend's
\code{\link[=backendInitialize]{backendInitialize()}} method. Parameter \code{backend} allows to specify which
\linkS4class{MsBackend} should be used for data representation and storage.
}

\section{Data representation of a \code{Spectra}}{


The MS data which can be accessed through the \code{Spectra} object is
\emph{represented} by its backend, which means that this backend defines how
and where the data is stored (e.g. in memory or on disk). The \code{Specrta}
object relies on the backend to provide the MS data whenever it needs it
for data processing.
Different backends with different properties, such as minimal memory
requirement or fast data access, are defined in the \emph{Spectra} package or
one of the MsBackend* packages. More information on backends and their
properties is provided in the documentation of \link{MsBackend}.

On-disk backends keep only a limited amount of data in memory retrieving
most of the data (usually the MS peak data) upon request on-the-fly from
their on-disk data representations. Moving the on-disk data storage of such
a backend or a serialized object to a different location in the file
system will cause data corruption. The \code{dataStorageBasePath()} and
\verb{dataStorageBasePath<-} functions allow in such cases (and if thebackend
classes support this operation), to get or change the \emph{base}
path to the directory of the backend's data storage. In-memory backends
such as \link{MsBackendMemory} or \link{MsBackendDataFrame} keeping all MS data in
memory don't support, and need, this function, but for \link{MsBackendMzR} this
function can be used to update/adapt the path to the directory containing
the original data files. Thus, for \code{Spectra} objects (using this backend)
that were moved to another file system or computer, these functions allow to
adjust/adapt the base file path.
}

\section{Changing data representation of a \code{Spectra}}{


The data representation, i.e. the backend of a \code{Spectra} object can be
changed with the \code{setBackend()} method that takes an instance of the new
backend as second parameter \code{backend}. A call to
\code{setBackend(sps, backend = MsBackendDataFrame())}
would for example change the backend of \code{sps} to the \emph{in-memory}
\code{MsBackendDataFrame}. Changing to a backend is only supported if that
backend has a \code{data} parameter in its \code{backendInitialize()} method and if
\code{supportsSetBackend()} returns \code{TRUE} for that backend. \code{setBackend()} will
transfer the full spectra data from the originating backend as a \code{DataFrame}
to the new backend.

Generally, it is not possible to change \strong{to} a read-only backend such as
the \code{\link[=MsBackendMzR]{MsBackendMzR()}} backend.

The definition of the function is:
\code{setBackend(object, backend, ..., f = dataStorage(object), BPPARAM = bpparam())} and its parameters are:
\itemize{
\item \code{object}: the \code{Spectra} object.
\item \code{backend}: an instance of the new backend, e.g. \verb{[MsBackendMemory()]}.
\item \code{f}: factor allowing to parallelize the change of the backends. By
default the process of copying the spectra data from the original to the
new backend is performed separately (and in parallel) for each file. Users
are advised to use the default setting.
\item \code{...}: optional additional arguments passed to the \code{\link[=backendInitialize]{backendInitialize()}}
method of the new \code{backend}.
\item \code{BPPARAM}: setup for the parallel processing. See
\code{\link[BiocParallel:register]{BiocParallel::bpparam()}} for details.
}
}

\section{Exporting data from a \code{Spectra} object}{


Data from a \code{Spectra} object can be \strong{exported} to a file with the
\code{export()} function. The actual export of the data is performed by
the \code{export} method of the \link{MsBackend} class defined with the mandatory
parameter \code{backend} which defines also the format in which the data
is exported. Note however that not all backend classes support
export of data. From the \code{MsBackend} classes in the \code{Spectra} package
currently only the \code{MsBackendMzR} backend supports data export (to
mzML/mzXML file(s)); see the help page of the \linkS4class{MsBackend} for
information on its arguments or the examples below or the vignette
for examples.

The definition of the function is
\code{export(object, backend,  ...)} and its
parameters are:
\itemize{
\item \code{object}: the \code{Spectra} object to be exported.
\item \code{backend}: instance of a class extending \link{MsBackend} which supports export
of the data (i.e. which has a defined \code{export} method).
\item \code{...}: additional parameters specific for the \code{MsBackend} passed with
parameter \code{backend}.
}
}

\examples{

##  --------  CREATION OF SPECTRA OBJECTS  --------

## Create a Spectra providing a `DataFrame` containing the spectrum data.

spd <- DataFrame(msLevel = c(1L, 2L), rtime = c(1.1, 1.2))
spd$mz <- list(c(100, 103.2, 104.3, 106.5), c(45.6, 120.4, 190.2))
spd$intensity <- list(c(200, 400, 34.2, 17), c(12.3, 15.2, 6.8))

data <- Spectra(spd)
data

## Create a Spectra from mzML files and use the `MsBackendMzR` on-disk
## backend.
sciex_file <- dir(system.file("sciex", package = "msdata"),
    full.names = TRUE)
sciex <- Spectra(sciex_file, backend = MsBackendMzR())
sciex


##  --------  CHANGING DATA REPRESENTATIONS  --------

## The MS data is on disk and will be read into memory on-demand. We can
## however change the backend to a MsBackendMemory backend which will
## keep all of the data in memory.
sciex_im <- setBackend(sciex, MsBackendMemory())
sciex_im

## The `MsBackendMemory()` supports the `setBackend()` method:
supportsSetBackend(MsBackendMemory())

## Thus, it is possible to change to that backend with `setBackend()`. Most
## read-only backends however don't support that, such as the
## `MsBackendMzR` and `setBackend()` would fail to change to that backend.
supportsSetBackend(MsBackendMzR())

## The on-disk object `sciex` is light-weight, because it does not keep the
## MS peak data in memory. The `sciex_im` object in contrast keeps all the
## data in memory and its size is thus much larger.
object.size(sciex)
object.size(sciex_im)

## The spectra variable `dataStorage` returns for each spectrum the location
## where the data is stored. For in-memory objects:
head(dataStorage(sciex_im))

## While objects that use an on-disk backend will list the files where the
## data is stored.
head(dataStorage(sciex))

## The spectra variable `dataOrigin` returns for each spectrum the *origin*
## of the data. If the data is read from e.g. mzML files, this will be the
## original mzML file name:
head(dataOrigin(sciex))
head(dataOrigin(sciex_im))


##  -------- DATA EXPORT  --------

## Some `MsBackend` classes provide an `export()` method to export the data
## to the file format supported by the backend.
## The `MsBackendMzR` for example allows to export MS data to mzML or
## mzXML file(s), the `MsBackendMgf` (defined in the MsBackendMgf R package)
## would allow to export the data in mgf file format.
## Below we export the MS data in `data`. We call the `export()` method on
## this object, specify the backend that should be used to export the data
## (and which also defines the output format) and provide a file name.
fl <- tempfile()
export(data, MsBackendMzR(), file = fl)

## This exported our data in mzML format. Below we read the first 6 lines
## from that file.
readLines(fl, n = 6)

## If only a single file name is provided, all spectra are exported to that
## file. To export data with the `MsBackendMzR` backend to different files, a
## file name for each individual spectrum has to be provided.
## Below we export each spectrum to its own file.
fls <- c(tempfile(), tempfile())
export(data, MsBackendMzR(), file = fls)

## Reading the data from the first file
res <- Spectra(backendInitialize(MsBackendMzR(), fls[1]))

mz(res)
mz(data)
}
\author{
Sebastian Gibb, Johannes Rainer, Laurent Gatto, Philippine Louail
}
