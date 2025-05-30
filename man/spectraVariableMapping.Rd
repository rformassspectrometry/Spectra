% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{spectraVariableMapping}
\alias{spectraVariableMapping}
\alias{spectraVariableMapping<-}
\title{Mapping between spectra variables and data file fields}
\usage{
spectraVariableMapping(object, ...)

spectraVariableMapping(object, ...) <- value
}
\arguments{
\item{object}{An instance of an object extending \code{\link[=MsBackend]{MsBackend()}}.}

\item{...}{Optional parameters.}

\item{value}{For \verb{spectraVariableMapping<-}: a named \code{character} vector.}
}
\value{
A named \code{character} with names being spectra variable names (use
\code{\link[=spectraVariables]{spectraVariables()}} for a list of supported names) and values being the
data field names.
}
\description{
The \code{spectraVariableMapping} function provides the mapping
between \emph{spectra variables} of a \code{\link[=Spectra]{Spectra()}} object with data fields from a
data file. Such name mapping is expected to enable an easier import of data
files with specific \emph{dialects}, e.g. files in MGF format that use a
different naming convention for core spectra variables.

\code{\link[=MsBackend]{MsBackend()}} implementations are expected to implement this function
(if needed) to enable import of data from file formats with non-standardized
data fields.
}
\author{
Johannes Rainer
}
