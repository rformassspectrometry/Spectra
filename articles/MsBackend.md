# Creating new \`MsBackend\` classes

**Package**:
*[Spectra](https://bioconductor.org/packages/3.23/Spectra)*  
**Authors**: RforMassSpectrometry Package Maintainer \[cre\], Laurent
Gatto \[aut\] (ORCID: <https://orcid.org/0000-0002-1520-2268>), Johannes
Rainer \[aut\] (ORCID: <https://orcid.org/0000-0002-6977-7147>),
Sebastian Gibb \[aut\] (ORCID: <https://orcid.org/0000-0001-7406-4443>),
Philippine Louail \[aut\] (ORCID:
<https://orcid.org/0009-0007-5429-6846>), Jan Stanstrup \[ctb\] (ORCID:
<https://orcid.org/0000-0003-0541-7369>), Nir Shahaf \[ctb\], Mar
Garcia-Aloy \[ctb\] (ORCID: <https://orcid.org/0000-0002-1330-6610>),
Guillaume Deflandre \[ctb\] (ORCID:
<https://orcid.org/0009-0008-1257-2416>), Ahlam Mentag \[ctb\] (ORCID:
<https://orcid.org/0009-0008-5438-7067>)  
**Last modified:** 2025-12-23 13:01:48.121851  
**Compiled**: Tue Dec 23 13:21:43 2025

## Introduction

This vignette briefly describes the `MsBackend` class which is used by
the `Spectra` package to *represent* and provide Mass Spectrometry (MS)
data and illustrates how a new such *backend* class can be created and
tested for validity.

Contributions to this vignette (content or correction of typos) or
requests for additional details and information are highly welcome
(ideally *via* pull requests or github issues).

## What is a `MsBackend`?

The `Spectra` package separates the code for the analysis of MS data
from the code needed to import, represent and provide the data. The
former is implemented for the `Spectra` class which is the main object
users will use for their analyses. The `Spectra` object relies on a
so-called *backend* to provide the MS data. The `MsBackend` virtual
class defines the API that new *backend* classes need to implement in
order to be used with the `Spectra` object. Each `Spectra` object
contains an implementation of such a `MsBackend` within its `@backend`
slot which provides the MS data to the `Spectra` object. All data
management is thus hidden from the user. In addition this separation
allows to define new, alternative, data representations and integrate
them seamlessly into a `Spectra`-based data analysis workflow.

This concept is an extension of the of *in-memory* and *on-disk* data
representations from the
*[MSnbase](https://bioconductor.org/packages/3.23/MSnbase)* package
(Gatto et al. 2020).

### Conventions and definitions

General conventions for MS data of a `Spectra` are:

- One `Spectra` object is supposed to contain MS (spectral) data of
  multiple MS spectra.
- m/z values within each spectrum are expected to be sorted
  increasingly.
- Missing values (`NA`) for m/z values are not supported.
- Properties of a spectrum are called *spectra variables*. While
  backends can define their own properties, a minimum required set of
  spectra variables **must** be provided by each backend (even if their
  values are empty). These *core spectra variables* are listed (along
  with their expected data type) by the
  [`coreSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
  function.
- `dataStorage` and `dataOrigin` are two special spectra variables that
  define for each spectrum where the data is stored and from where the
  data derived (or was loaded, such as the data origin). Both are
  expected to be of type`character` and need to be defined by the
  backend (i.e., they can not be empty or missing).
- `MsBackend` implementations can also represent purely *read-only* data
  resources. In this case only data accessor methods need to be
  implemented but not data replacement methods. Whether a backend is
  read-only can be set with the `@readonly` slot of the virtual
  `MsBackend` class (the `isReadOnly()` function can be used to retrieve
  the value for this slot). The default is `@readonly = FALSE` and thus
  all data replacement method listed in section *Data replacement
  methods* have to be implemented. For read-only backends
  (`@readonly = TRUE`) only the methods in section *Required methods*
  need to be implemented. Backends can also be *partially* read-only,
  such as the `MsBackendMzR`. This backend allows for example to change
  spectra variables, but not the peaks data (i.e. the m/z and intensity
  values). Also, backends for purely read-only resources could extend
  the `MsBackendCached` from the
  *[Spectra](https://bioconductor.org/packages/3.23/Spectra)* package to
  enable support for modifying (or adding) spectra variables. Any
  changes to spectra variables will be internally cached by the
  `MsBackendCached` without the need of them being propagating to the
  underlying data resource (see for example the `MsBackendMassbankSql`
  from the
  *[MsBackendMassbank](https://bioconductor.org/packages/3.23/MsBackendMassbank)*
  package).

### Notes on parallel processing

For parallel processing, `Spectra` splits the backend based on a defined
`factor` and processes each in parallel (or *in serial* if a
`SerialParam` is used). The splitting `factor` can be defined for
`Spectra` by setting the parameter `processingChunkSize`. Alternatively,
through the
[`backendParallelFactor()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
function the backend can also *suggest* a `factor` that should/could be
used for splitting and parallel processing. The default implementation
for
[`backendParallelFactor()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
is to return an empty `factor`
([`factor()`](https://rdrr.io/r/base/factor.html)) hence not suggesting
any preferred splitting.
[`backendParallelFactor()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
for `MsBackendMzR` on the other hand returns a `factor` based on the
data files the data is stored in (i.e. based on the `dataStorage` of the
MS data).

Besides parallel processing, this chunk-wise processing can also reduce
the memory demand for operations, because only the peak data of the
current chunk needs to be realized in memory.

## API

The `MsBackend` class defines core methods that have to be implemented
by a MS *backend* as well as *optional* methods with default
implementations that might be implemented for a new backend but don’t
necessarily have to. These functions are described in sections *Required
methods* and *Optional methods*, respectively.

To create a new backend a class extending the virtual `MsBackend` needs
to be implemented. In the example below we create thus a simple class
with a `data.frame` to contain general spectral properties (*spectra
variables*) and two slots for m/z and intensity values. These are stored
as `NumericList` objects since both m/z and intensity values are
expected to be of type `numeric` and to allow to store data from
multiple spectra into a single backend object. We also define a simple
constructor function that returns an empty instance of our new class.

``` r

library(Spectra)
library(IRanges)

setClass("MsBackendTest",
         contains = "MsBackend",
         slots = c(
             spectraVars = "data.frame",
             mz = "NumericList",
             intensity = "NumericList"
         ),
         prototype = prototype(
             spectraVars = data.frame(),
             mz = NumericList(compress = FALSE),
             intensity = NumericList(compress = FALSE)
         ))

MsBackendTest <- function() {
    new("MsBackendTest")
}
```

The 3 slots `spectraVars`, `mz` and `intensity` will be used to store
our MS data, each row in `spectraVars` being data for one spectrum with
the columns being the different *spectra variables* (i.e. additional
properties of a spectrum such as its retention time or MS level) and
each element in `mz` and `intensity` being a `numeric` with the m/z and
intensity values of the respective spectrum.

We should ideally also add some basic validity function that ensures the
data to be OK. The function below simply checks that the number of rows
of the `spectraVars` slot matches the length of the `mz` and `intensity`
slot.

``` r

setValidity("MsBackendTest", function(object) {
    if (length(object@mz) != length(object@intensity) ||
        length(object@mz) != nrow(object@spectraVars))
        return("length of 'mz' and 'intensity' has to match the number of ",
               "rows of 'spectraVars'")
    NULL
})
```

    ## Class "MsBackendTest" [in ".GlobalEnv"]
    ## 
    ## Slots:
    ##                                                                   
    ## Name:  spectraVars          mz   intensity    readonly     version
    ## Class:  data.frame NumericList NumericList     logical   character
    ## 
    ## Extends: "MsBackend"

We can now create an instance of our new class with the
`MsBackendTest()` function.

``` r

MsBackendTest()
```

    ## An object of class "MsBackendTest"
    ## Slot "spectraVars":
    ## data frame with 0 columns and 0 rows
    ## 
    ## Slot "mz":
    ## NumericList of length 0
    ## 
    ## Slot "intensity":
    ## NumericList of length 0
    ## 
    ## Slot "readonly":
    ## [1] FALSE
    ## 
    ## Slot "version":
    ## [1] "0.1"

Note that a *backend* class does not necessarily need to contain all the
data like the one from our example. Backends such as the `MsBackendMzR`
for example retrieve the data on the fly from the raw MS data files or
the `MsBackendSql` from the
*[MsBackendSql](https://bioconductor.org/packages/3.23/MsBackendSql)* a
SQL database.

### Required methods

Methods listed in this section must be implemented for a new class
extending `MsBackend`. Methods should ideally also implemented in the
order they are listed here. Also, it is strongly advised to write
dedicated unit tests for each newly implemented method or function
already **during** the development.

#### `spectraData()`

The
[`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method should return the **full** spectra data within a backend as a
`DataFrame` object (defined in the
*[S4Vectors](https://bioconductor.org/packages/3.23/S4Vectors)*
package). The second parameter `columns` allows to define the names of
the spectra variables that should be returned in the `DataFrame`. Each
row in this data frame should represent one spectrum, each column a
spectra variable. Columns `"mz"` and `"intensity"` (if requested) have
to contain each a `NumericList` with the m/z and intensity values of the
spectra. The `DataFrame` **must** provide values (even if they are `NA`)
for **all** requested spectra variables of the backend (**including**
the core spectra variables).

In our toy backend class we keep the spectra variable data in a simple
`data.frame` without any constraints such as required columns etc. To
ensure that
[`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
always returns all required *core* spectra variables (of the correct
data type) we can use however the
[`fillCoreSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/fillCoreSpectraVariables.md)
function. This function adds eventually missing core spectra variables
to a `data.frame` thus we use it in our implementation of the
[`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method below:

``` r

setMethod(
    "spectraData", "MsBackendTest",
    function(object, columns = spectraVariables(object)) {
        if (!all(columns %in% spectraVariables(object)))
            stop("Some of the requested spectra variables are not available")
        res <- fillCoreSpectraVariables(
            object@spectraVars,
            columns = columns[!columns %in% c("mz", "intensity")])
        ## Add m/z and intensity values to the result
        res$mz <- object@mz
        res$intensity <- object@intensity
        as(res[, columns, drop = FALSE], "DataFrame")
})
```

As an alternative, we could also initialize the `@spectraVars` data
frame within the
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method adding columns for spectra variables that are not provided by the
user and require that this data frame always contains all core spectra
variables. Extracting spectra data (single spectra variables or the full
data) might thus be more efficient then the on-the-fly initialization
with eventual missing spectra variables, but the backend class would
also have a larger memory footprint because even spectra variables with
only missing values for all spectra need to be stored within the object.

#### `spectraVariables()`

The
[`spectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method should return a `character` vector with the names of all
available spectra variables of the backend. While a backend class should
support defining and providing their own spectra variables, each
`MsBackend` class **must** provide also the *core spectra variables* (in
the correct data type). Since not all data file formats provide values
for all these spectra variables they can however also be `NA` (with the
exception of the spectra variable `"dataStorage"`).

The
[`coreSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
function returns the full list of mandatory spectra variables along with
their expected data type.

``` r

coreSpectraVariables()
```

    ##                 msLevel                   rtime          acquisitionNum 
    ##               "integer"               "numeric"               "integer" 
    ##               scanIndex                      mz               intensity 
    ##               "integer"           "NumericList"           "NumericList" 
    ##             dataStorage              dataOrigin              centroided 
    ##             "character"             "character"               "logical" 
    ##                smoothed                polarity             precScanNum 
    ##               "logical"               "integer"               "integer" 
    ##             precursorMz      precursorIntensity         precursorCharge 
    ##               "numeric"               "numeric"               "integer" 
    ##         collisionEnergy  isolationWindowLowerMz isolationWindowTargetMz 
    ##               "numeric"               "numeric"               "numeric" 
    ##  isolationWindowUpperMz 
    ##               "numeric"

A typical
[`spectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method for a `MsBackend` class will thus be implemented similarly to the
one for our `MsBackendTest` test backend: it will return the union of
the core spectra variables and the names for all available spectra
variables within the backend object.

``` r

setMethod("spectraVariables", "MsBackendTest", function(object) {
    union(names(coreSpectraVariables()), colnames(object@spectraVars))
})
```

Note that we avoided calling
[`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
in the definition of the
[`spectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method, because in our implementation of the
[`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method above we were already using
[`spectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
to get the available spectra variables. Calling
[`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
also in the implementation of the
[`spectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method for `MsBackendTest` would have resulted in a cyclic function
call.

#### `backendInitialize()`

The
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method is expected to be called after creating an instance of the
backend class and should prepare (initialize) the backend which in most
cases means that MS data is loaded. This method can take any parameters
needed by the backend to get loaded/initialized with data (which can be
file names from which to load the data, a database connection or
object(s) containing the data). During
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
usually also the special spectra variables `dataStorage` and
`dataOrigin` are set.

Below we define a
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method that takes as arguments a `data.frame` with spectra variables and
two `list`s with the m/z and intensity values for each spectrum.

``` r

setMethod(
    "backendInitialize", "MsBackendTest",
    function(object, svars, mz, intensity) {
        if (!is.data.frame(svars))
            stop("'svars' needs to be a 'data.frame' with spectra variables")
        if (is.null(svars$dataStorage))
            svars$dataStorage <- "<memory>"
        if (is.null(svars$dataOrigin))
            svars$dataOrigin <- "<user provided>"
        object@spectraVars <- svars
        object@mz <- NumericList(mz, compress = FALSE)
        object@intensity <- NumericList(intensity, compress = FALSE)
        validObject(object)
        object
    })
```

In addition to adding the data to object, the function also defined the
`dataStorage` and `dataOrigin` spectra variables. The purpose of these
two variables is to provide some information on where the data is stored
(*in memory* as in our example) and from where the data is originating.
The `dataOrigin` would for example allow to specify from which original
data files individual spectra derive.

We can now create an instance of our backend class and fill it with
data. We thus first define our MS data and pass this to the
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method.

``` r

## A data.frame with spectra variables.
svars <- data.frame(msLevel = c(1L, 2L, 2L),
                    rtime = c(1.2, 1.3, 1.4))
## m/z values for each spectrum.
mzs <- list(c(12.3, 13.5, 16.5, 17.5),
            c(45.1, 45.2),
            c(64.4, 123.1, 124.1))
## intensity values for each spectrum.
ints <- list(c(123.3, 153.6, 2354.3, 243.4),
             c(100, 80.1),
             c(12.3, 35.2, 100))

## Create and initialize the backend
be <- backendInitialize(MsBackendTest(),
                        svars = svars, mz = mzs, intensity = ints)
be
```

    ## An object of class "MsBackendTest"
    ## Slot "spectraVars":
    ##   msLevel rtime dataStorage      dataOrigin
    ## 1       1   1.2    <memory> <user provided>
    ## 2       2   1.3    <memory> <user provided>
    ## 3       2   1.4    <memory> <user provided>
    ## 
    ## Slot "mz":
    ## NumericList of length 3
    ## [[1]] 12.3 13.5 16.5 17.5
    ## [[2]] 45.1 45.2
    ## [[3]] 64.4 123.1 124.1
    ## 
    ## Slot "intensity":
    ## NumericList of length 3
    ## [[1]] 123.3 153.6 2354.3 243.4
    ## [[2]] 100 80.1
    ## [[3]] 12.3 35.2 100
    ## 
    ## Slot "readonly":
    ## [1] FALSE
    ## 
    ## Slot "version":
    ## [1] "0.1"

While this method works and is compliant with the `MsBackend` API
(because there is no requirement on the input parameters for the
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method), it would be good practice for backends that are supposed to
support replacing data, to add an optional additional parameter `data`
that would allow passing the *complete* MS data (including m/z and
intensity values) to the function as a `DataFrame`. This would simplify
the implementation of some replacement methods and would in addition
also allow to change the backend of a `Spectra` using the
[`setBackend()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
to our new backend. We thus re-implement the
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method supporting also to initialize the backend with such a data frame
and we also implement a helper function that checks spectra variables
for the correct data type.

``` r

#' Helper function to check if core spectra variables have the correct
#' data type.
#'
#' @param x `data.frame` with the data for spectra variables.
#'
#' @param name `character` defining the column names (spectra variables) of `x`
#'     for which the correct data type should be evaluated.
.sv_valid_data_type <- function(x, name = colnames(x)) {
    sv <- coreSpectraVariables()[names(coreSpectraVariables()) %in% name]
    for (i in seq_along(sv)) {
        if (!is(x[, names(sv[i])], sv[i]))
            stop("Spectra variabe \"", names(sv[i]), "\" is not of type ",
                 sv[i], call. = FALSE)
    }
    TRUE
}
```

This function is then used to check the input data in our new
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method.

``` r

setMethod(
    "backendInitialize", "MsBackendTest",
    function(object, svars, mz, intensity, data) {
        if (!missing(data)) {
            svars <- as.data.frame(
                data[, !colnames(data) %in% c("mz", "intensity")])
            if (any(colnames(data) == "mz"))
                mz <- data$mz
            if (any(colnames(data) == "intensity"))
                intensity <- data$intensity
        }
        if (!is.data.frame(svars))
            stop("'svars' needs to be a 'data.frame' with spectra variables")
        if (is.null(svars$dataStorage))
            svars$dataStorage <- "<memory>"
        if (is.null(svars$dataOrigin))
            svars$dataOrigin <- "<user provided>"
        .sv_valid_data_type(svars)
        object@spectraVars <- svars
        object@mz <- NumericList(mz, compress = FALSE)
        object@intensity <- NumericList(intensity, compress = FALSE)
        validObject(object)
        object
    })
```

We below create the backend again with the updated
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md).

``` r

## Create and initialize the backend
be <- backendInitialize(MsBackendTest(),
                        svars = svars, mz = mzs, intensity = ints)
be
```

    ## An object of class "MsBackendTest"
    ## Slot "spectraVars":
    ##   msLevel rtime dataStorage      dataOrigin
    ## 1       1   1.2    <memory> <user provided>
    ## 2       2   1.3    <memory> <user provided>
    ## 3       2   1.4    <memory> <user provided>
    ## 
    ## Slot "mz":
    ## NumericList of length 3
    ## [[1]] 12.3 13.5 16.5 17.5
    ## [[2]] 45.1 45.2
    ## [[3]] 64.4 123.1 124.1
    ## 
    ## Slot "intensity":
    ## NumericList of length 3
    ## [[1]] 123.3 153.6 2354.3 243.4
    ## [[2]] 100 80.1
    ## [[3]] 12.3 35.2 100
    ## 
    ## Slot "readonly":
    ## [1] FALSE
    ## 
    ## Slot "version":
    ## [1] "0.1"

The
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method that we implemented for our backend class expects the user to
provide the full MS data. This does however not always have to be the
case. The
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method of the `MsBackendMzR` backend takes for example the file names of
the raw mzML, mzXML or CDF files as input and initializes the backend by
importing part of the data from these. Also the backends defined by the
*[MsBackendMgf](https://bioconductor.org/packages/3.23/MsBackendMgf)* or
`r Biocpkg("MsBackendMsp")` packages work in the same way and thus allow
to import MS data from these specific file formats. The
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method of the backend defined in the
*[MsBackendSql](https://bioconductor.org/packages/3.23/MsBackendSql)* on
the other hand takes only the connection to a database containing the
data as input and performs some sanity checks on the data but does not
load the data into the backend. Any subsequent data access is handled by
the methods of the backend class through SQL calls to the database.

The purpose of the
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method is to *initialize* and prepare the data in a way that it can be
accessed by a `Spectra` object (through the initialized backend class).
Whether the data is loaded by the
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method into memory or simply referenced to within the backend class does
not matter as long as the backend is able to provide the data with its
accessor methods.

Note also that a
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
function should ideally also perform some data sanity checks
(e.g. whether spectra variables have the correct data type etc).

After creating our backend and initializing it with data we can also
test and use the
[`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method:

``` r

## Full data
spectraData(be)
```

    ## DataFrame with 3 rows and 19 columns
    ##     msLevel     rtime acquisitionNum scanIndex                 mz
    ##   <integer> <numeric>      <integer> <integer>      <NumericList>
    ## 1         1       1.2             NA        NA 12.3,13.5,16.5,...
    ## 2         2       1.3             NA        NA          45.1,45.2
    ## 3         2       1.4             NA        NA   64.4,123.1,124.1
    ##                  intensity dataStorage      dataOrigin centroided  smoothed
    ##              <NumericList> <character>     <character>  <logical> <logical>
    ## 1  123.3, 153.6,2354.3,...    <memory> <user provided>         NA        NA
    ## 2              100.0, 80.1    <memory> <user provided>         NA        NA
    ## 3         12.3, 35.2,100.0    <memory> <user provided>         NA        NA
    ##    polarity precScanNum precursorMz precursorIntensity precursorCharge
    ##   <integer>   <integer>   <numeric>          <numeric>       <integer>
    ## 1        NA          NA          NA                 NA              NA
    ## 2        NA          NA          NA                 NA              NA
    ## 3        NA          NA          NA                 NA              NA
    ##   collisionEnergy isolationWindowLowerMz isolationWindowTargetMz
    ##         <numeric>              <numeric>               <numeric>
    ## 1              NA                     NA                      NA
    ## 2              NA                     NA                      NA
    ## 3              NA                     NA                      NA
    ##   isolationWindowUpperMz
    ##                <numeric>
    ## 1                     NA
    ## 2                     NA
    ## 3                     NA

``` r

## Selected variables
spectraData(be, c("rtime", "mz", "centroided"))
```

    ## DataFrame with 3 rows and 3 columns
    ##       rtime                 mz centroided
    ##   <numeric>      <NumericList>  <logical>
    ## 1       1.2 12.3,13.5,16.5,...         NA
    ## 2       1.3          45.1,45.2         NA
    ## 3       1.4   64.4,123.1,124.1         NA

``` r

## Only missing core spectra variables
spectraData(be, c("centroided", "polarity"))
```

    ## DataFrame with 3 rows and 2 columns
    ##   centroided  polarity
    ##    <logical> <integer>
    ## 1         NA        NA
    ## 2         NA        NA
    ## 3         NA        NA

And the
[`spectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method:

``` r

spectraVariables(be)
```

    ##  [1] "msLevel"                 "rtime"                  
    ##  [3] "acquisitionNum"          "scanIndex"              
    ##  [5] "mz"                      "intensity"              
    ##  [7] "dataStorage"             "dataOrigin"             
    ##  [9] "centroided"              "smoothed"               
    ## [11] "polarity"                "precScanNum"            
    ## [13] "precursorMz"             "precursorIntensity"     
    ## [15] "precursorCharge"         "collisionEnergy"        
    ## [17] "isolationWindowLowerMz"  "isolationWindowTargetMz"
    ## [19] "isolationWindowUpperMz"

#### `peaksData()`

The
[`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method extracts the MS peaks data from a backend, which includes the m/z
and intensity values of each MS peak of a spectrum. These are expected
to be returned as a `List` of numerical matrices with columns in each
`matrix` being the requested *peaks variables* (with the default being
`"mz"` and `"intensity"`) of one spectrum. Backends must provide at
least these two peaks variables.

Below we implement the
[`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method for our backend. We need to loop over the `@mz` and `@intensity`
slots to merge the m/z and intensity of each spectrum into a `matrix`.
Also, for simplicity reasons, we accept only `c("mz", "intensity")` for
the `columns` parameter. This is the expected default behavior for a
`MsBackend`, but in general the `columns` parameter is thought to allow
the user to specify which peaks variables should be returned in each
`matrix`.

``` r

setMethod(
    "peaksData", "MsBackendTest",
    function(object, columns = c("mz", "intensity")) {
        if (length(columns) != 2 && columns != c("mz", "intensity"))
            stop("'columns' supports only \"mz\" and \"intensity\"")
        mapply(mz = object@mz, intensity = object@intensity,
               FUN = cbind, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    })
```

And with this method we can now extract the peaks data from our backend.

``` r

peaksData(be)
```

    ## [[1]]
    ##        mz intensity
    ## [1,] 12.3     123.3
    ## [2,] 13.5     153.6
    ## [3,] 16.5    2354.3
    ## [4,] 17.5     243.4
    ## 
    ## [[2]]
    ##        mz intensity
    ## [1,] 45.1     100.0
    ## [2,] 45.2      80.1
    ## 
    ## [[3]]
    ##         mz intensity
    ## [1,]  64.4      12.3
    ## [2,] 123.1      35.2
    ## [3,] 124.1     100.0

The
[`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method is used in many data analysis functions of the `Spectra` object
to extract the MS data, thus ideally this method should be implemented
in an efficient way. For our backend we need to loop over the lists of
m/z and intensity values which is obviously not ideal. Thus, storing the
m/z and intensity values in separate slots as done in this backend might
not be ideal. The `MsBackendMemory` backend for example stores the MS
data already as a `list` of matrices which results in a more efficient
[`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method (but comes also with a larger overhead when adding, replacing or
checking MS data).

Note also that while a backend needs to provide m/z and intensity
values, additional peak variables would also be supported. The
`MsBackendMemory` class for example allows to store and provide
additional peak variables that can then be added as additional columns
to each returned `matrix`. In this case the default
[`peaksVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method should also be overwritten to list the additionally available
variables and the `columns` parameter of the
[`peaksData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method should allow selection of these additional peaks variables (in
addition to the required `"mz"` and `"intensity"` variables).

#### `extractByIndex()` and `[`

The
[`extractByIndex()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
and `[` methods allows to subset `MsBackend` objects. This operation is
expected to reduce a `MsBackend` object to the selected spectra. These
methods must also support duplication (e.g. `[c(1, 1, 1)]` and
extraction in any arbitrary order (e.g. `[c(3, 1, 5, 3)]`). While both
methods subset the object,
[`extractByIndex()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
only supports to subset with an `integer` index, while `[`, to be
compliant with the base R implementation, should support to subset by
indices or logical vectors. An error should be thrown if indices are out
of bounds, but the method should also support returning an empty backend
with `[integer()]`. Note that the
[`MsCoreUtils::i2index`](https://rdrr.io/pkg/MsCoreUtils/man/i2index.html)
function can be used to check for correct input (and convert the input
to an `integer` index).

The
[`extractByIndex()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method is used by the data operation and analysis methods on `Spectra`
objects, while the `[` is intended to be used by the end user (if
needed). Below we implement
[`extractByIndex()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
for our backend:

``` r

setMethod("extractByIndex", c("MsBackendTest", "ANY"), function(object, i) {
    object@spectraVars <- object@spectraVars[i, ]
    object@mz <- object@mz[i]
    object@intensity <- object@intensity[i]
    object
})
```

The `[` does not need to be defined because a default implementation for
the base `MsBackend` exists.

We can now subset our backend to the last two spectra.

``` r

a <- extractByIndex(be, 2:3)
spectraData(a)
```

    ## DataFrame with 2 rows and 19 columns
    ##     msLevel     rtime acquisitionNum scanIndex                mz
    ##   <integer> <numeric>      <integer> <integer>     <NumericList>
    ## 1         2       1.3             NA        NA         45.1,45.2
    ## 2         2       1.4             NA        NA  64.4,123.1,124.1
    ##           intensity dataStorage      dataOrigin centroided  smoothed  polarity
    ##       <NumericList> <character>     <character>  <logical> <logical> <integer>
    ## 1       100.0, 80.1    <memory> <user provided>         NA        NA        NA
    ## 2  12.3, 35.2,100.0    <memory> <user provided>         NA        NA        NA
    ##   precScanNum precursorMz precursorIntensity precursorCharge collisionEnergy
    ##     <integer>   <numeric>          <numeric>       <integer>       <numeric>
    ## 1          NA          NA                 NA              NA              NA
    ## 2          NA          NA                 NA              NA              NA
    ##   isolationWindowLowerMz isolationWindowTargetMz isolationWindowUpperMz
    ##                <numeric>               <numeric>              <numeric>
    ## 1                     NA                      NA                     NA
    ## 2                     NA                      NA                     NA

Or extracting the second spectrum multiple times.

``` r

a <- be[c(2, 2, 2)]
spectraData(a)
```

    ## DataFrame with 3 rows and 19 columns
    ##       msLevel     rtime acquisitionNum scanIndex            mz     intensity
    ##     <integer> <numeric>      <integer> <integer> <NumericList> <NumericList>
    ## 2           2       1.3             NA        NA     45.1,45.2   100.0, 80.1
    ## 2.1         2       1.3             NA        NA     45.1,45.2   100.0, 80.1
    ## 2.2         2       1.3             NA        NA     45.1,45.2   100.0, 80.1
    ##     dataStorage      dataOrigin centroided  smoothed  polarity precScanNum
    ##     <character>     <character>  <logical> <logical> <integer>   <integer>
    ## 2      <memory> <user provided>         NA        NA        NA          NA
    ## 2.1    <memory> <user provided>         NA        NA        NA          NA
    ## 2.2    <memory> <user provided>         NA        NA        NA          NA
    ##     precursorMz precursorIntensity precursorCharge collisionEnergy
    ##       <numeric>          <numeric>       <integer>       <numeric>
    ## 2            NA                 NA              NA              NA
    ## 2.1          NA                 NA              NA              NA
    ## 2.2          NA                 NA              NA              NA
    ##     isolationWindowLowerMz isolationWindowTargetMz isolationWindowUpperMz
    ##                  <numeric>               <numeric>              <numeric>
    ## 2                       NA                      NA                     NA
    ## 2.1                     NA                      NA                     NA
    ## 2.2                     NA                      NA                     NA

#### `backendMerge()`

The `backendMerge()` method merges (combines) `MsBackend` objects (of
the same type!) into a single instance. For our test backend we thus
need to combine the values in the `@spectraVars`, `@mz` and `@intensity`
slots. To support also merging of `data.frame`s with different set of
columns we use the
[`MsCoreUtils::rbindFill()`](https://rdrr.io/pkg/MsCoreUtils/man/rbindFill.html)
function instead of a simple
[`rbind()`](https://rdrr.io/pkg/BiocGenerics/man/cbind.html) (this
function joins data frames making an union of all available columns).

``` r

setMethod("backendMerge", "MsBackendTest", function(object, ...) {
    res <- object
    object <- unname(c(object, ...))
    res@mz <- do.call(c, lapply(object, function(z) z@mz))
    res@intensity <- do.call(c, lapply(object, function(z) z@intensity))
    res@spectraVars <- do.call(MsCoreUtils::rbindFill,
                               lapply(object, function(z) z@spectraVars))
    validObject(res)
    res
})
```

Again, this implementation which requires 3 loops might not be the most
efficient - but it allows to merge backends of the type `MsBackendTest`.

``` r

a <- backendMerge(be, be[2], be)
a
```

    ## An object of class "MsBackendTest"
    ## Slot "spectraVars":
    ##    msLevel rtime dataStorage      dataOrigin
    ## 1        1   1.2    <memory> <user provided>
    ## 2        2   1.3    <memory> <user provided>
    ## 3        2   1.4    <memory> <user provided>
    ## 21       2   1.3    <memory> <user provided>
    ## 11       1   1.2    <memory> <user provided>
    ## 22       2   1.3    <memory> <user provided>
    ## 31       2   1.4    <memory> <user provided>
    ## 
    ## Slot "mz":
    ## NumericList of length 7
    ## [[1]] 12.3 13.5 16.5 17.5
    ## [[2]] 45.1 45.2
    ## [[3]] 64.4 123.1 124.1
    ## [[4]] 45.1 45.2
    ## [[5]] 12.3 13.5 16.5 17.5
    ## [[6]] 45.1 45.2
    ## [[7]] 64.4 123.1 124.1
    ## 
    ## Slot "intensity":
    ## NumericList of length 7
    ## [[1]] 123.3 153.6 2354.3 243.4
    ## [[2]] 100 80.1
    ## [[3]] 12.3 35.2 100
    ## [[4]] 100 80.1
    ## [[5]] 123.3 153.6 2354.3 243.4
    ## [[6]] 100 80.1
    ## [[7]] 12.3 35.2 100
    ## 
    ## Slot "readonly":
    ## [1] FALSE
    ## 
    ## Slot "version":
    ## [1] "0.1"

#### `intensity()`

Extract the intensity values for each spectrum in the backend. The
result is expected to be a `NumericList` of length equal to the number
of spectra represented by the backend. For our test backend we can
simply return the `@intensity` slot since the data is already stored
within a `NumericList`.

``` r

setMethod("intensity", "MsBackendTest", function(object) {
    object@intensity
})
intensity(be)
```

    ## NumericList of length 3
    ## [[1]] 123.3 153.6 2354.3 243.4
    ## [[2]] 100 80.1
    ## [[3]] 12.3 35.2 100

#### `mz()`

Extract the m/z values for each spectrum in the backend. The result is
expected to be a `NumericList` of length equal to the number of spectra
represented by the backend. Also, the m/z values are expected to be
ordered increasingly for each element (spectrum).

``` r

setMethod("mz", "MsBackendTest", function(object) {
    object@mz
})
mz(be)
```

    ## NumericList of length 3
    ## [[1]] 12.3 13.5 16.5 17.5
    ## [[2]] 45.1 45.2
    ## [[3]] 64.4 123.1 124.1

#### `spectraNames()`

The
[`spectraNames()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
can be used to extract (optional) names (or IDs) for individual spectra
of a backend, or `NULL` if not set. For our test backend we can use the
`rownames` of the `@spectraVars` slot to store spectra names.

``` r

setMethod("spectraNames", "MsBackendTest", function(object) {
    rownames(object@spectraVars)
})
spectraNames(be)
```

    ## [1] "1" "2" "3"

These are all the methods that need to be implemented for a valid
*read-only* `MsBackend` class and running a test on such an object as
described in section *Testing the validity of the backend* should not
produce any errors. For backends that support also data replacement also
the methods listed in the next section need to be implemented.

### Data replacement methods

As stated in the general description, `MsBackend` implementations can
also be purely *read-only* resources allowing to just access data, but
not to replace the data. Thus, it is not strictly required to implement
these methods, but for a fully functional backend it is suggested (as
much as possible). A backend for a purely read-only MS data resource
might even extend the `MsBackendCached` backend defined in the `Spectra`
package that provides a mechanism to cache (spectra variable) data in a
`data.frame` within the object. The `MsBackendMassbankSql` implemented
in the
*[MsBackendMassbank](https://bioconductor.org/packages/3.23/MsBackendMassbank)*
package extends for example this backend and thus allows modifying some
spectra variables without changing the original data in the MassBank SQL
database.

#### `$<-`

The `$<-` method should allow to replace values for spectra variables or
also to add additional spectra variables to the backend. As with all
replacement methods, the
[`length()`](https://rdrr.io/r/base/length.html) of `value` has to match
the number of spectra represented by the backend. In addition to the
method we implement also a simple helper function `.match_length()` that
checks for the correct length of `value`. Each data replacement method
needs to check for that and this function thus reduces code duplication.

``` r

.match_length <- function(x, y) {
    if (length(x) != length(y))
        stop("Length of 'value' has to match the length of 'object'")
}

setReplaceMethod("$", "MsBackendTest", function(x, name, value) {
    .match_length(x, value)
    if (name == "mz") {
        mz(x) <- value
    } else if (name == "intensity") {
       intensity(x) <- value
    } else {
        x@spectraVars[[name]] <- value
    }
    .sv_valid_data_type(x@spectraVars, name)
    x
})
```

We can now replace for example existing spectra variables:

``` r

msLevel(be)
```

    ## [1] 1 2 2

``` r

be$msLevel <- c(2L, 1L, 2L)
msLevel(be)
```

    ## [1] 2 1 2

Or even add new spectra variables.

``` r

be$new_var <- c("a", "b", "c")
be$new_var
```

    ## [1] "a" "b" "c"

For all *core* spectra variables default replacement methods exist that
use the `$<-` to replace values.

#### `spectraData<-`

The `spectraData<-` method should allow to replace the data within a
backend. The method should take a `DataFrame` with the full data as
input value and is expected to replace the **full** data within the
backend, i.e. all spectra variables as well as peak data. Also,
importantly, the number of spectra before and after calling the
`spectraData<-` method on an object has to be the same. For our
implementation we can make use of the optional parameter `data` that we
added to the
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method and that allows to fill a `MsBackendTest` object with the full
data.

``` r

setReplaceMethod("spectraData", "MsBackendTest", function(object, value) {
    if (!inherits(value, "DataFrame"))
        stop("'value' is expected to be a 'DataFrame'")
    if (length(object) && length(object) != nrow(value))
        stop("'value' has to be a 'DataFrame' with ", length(object), " rows")
    object <- backendInitialize(MsBackendTest(), data = value)
    object
})
```

To test this new method we extract the full spectra data, add an
additional column (spectra variable) and replace the data again.

``` r

d <- spectraData(be)
d$new_col <- c("a", "b", "c")

spectraData(be) <- d
be$new_col
```

    ## [1] "a" "b" "c"

#### `intensity<-`

The `intensity<-` method should allow to replace the intensity values of
all spectra in a backend. This method is expected to only replace the
*values* of the intensities, but must not change the number of
intensities (and hence peaks) of a spectrum (that could be done with the
`peaksData<-` method that allows to replace intensity **and** m/z values
at the same time). The `value` for the method should ideally be a
`NumericList` to ensure that all intensity values are indeed `numeric`.

``` r

setReplaceMethod("intensity", "MsBackendTest", function(object, value) {
    .match_length(object, value)
    if (!(is.list(value) || inherits(value, "NumericList")))
        stop("'value' has to be a list or NumericList")
    if (!all(lengths(value) == lengths(mz(object))))
        stop("lengths of 'value' has to match the number of peaks per spectrum")
    if (!inherits(value, "NumericList"))
        value <- NumericList(value, compress = FALSE)
    object@intensity <- value
    object
})
```

We could now use this method to replace the intensities in our backend
with modified intensities.

``` r

intensity(be)
```

    ## NumericList of length 3
    ## [[1]] 123.3 153.6 2354.3 243.4
    ## [[2]] 100 80.1
    ## [[3]] 12.3 35.2 100

``` r

intensity(be) <- intensity(be) - 10
intensity(be)
```

    ## NumericList of length 3
    ## [[1]] 113.3 143.6 2344.3 233.4
    ## [[2]] 90 70.1
    ## [[3]] 2.3 25.2 90

#### `mz<-`

The `mz<-` method should allow to replace the m/z values of all spectra
in a backend. The implementation can be the same as for the
`intensity<-` method. m/z values within each spectrum need to be
increasingly ordered. We thus also check that this is the case for the
provided m/z values. We take here the advantage that a efficient
[`is.unsorted()`](https://rdrr.io/pkg/BiocGenerics/man/is.unsorted.html)
implementation for `NumericList` is already available, which is faster
than e.g. calling `vapply(mz(be), is.unsorted, logical(1))`.

``` r

setReplaceMethod("mz", "MsBackendTest", function(object, value) {
    .match_length(object, value)
    if (!(is.list(value) || inherits(value, "NumericList")))
        stop("'value' has to be a list or NumericList")
    if (!all(lengths(value) == lengths(mz(object))))
        stop("lengths of 'value' has to match the number of peaks per spectrum")
    if (!inherits(value, "NumericList"))
        value <- NumericList(value, compress = FALSE)
    if (any(is.unsorted(value)))
        stop("m/z values need to be increasingly sorted within each spectrum")
    object@mz <- value
    object
})
```

#### `peaksData<-`

The `peaksData<-` should allow to replace the peaks data (m/z and
intensity values) of all spectra in a backend. In contrast to the `mz<-`
and `intensity<-` methods this method should also support changing the
number of peaks per spectrum (e.g. due to filtering). Parameter `value`
has to be a `list` of `matrix` objects with columns `"mz"` and
`"intensity"`. The length of this list has to match the number of
spectra in the backend. In the implementation for our backend class we
need to loop over this list to extract the m/z and intensity values and
assign them to the `@mz` and `@intensity` slots.

``` r

setReplaceMethod("peaksData", "MsBackendTest", function(object, value) {
    if (!(is.list(value) || inherits(value, "SimpleList")))
        stop("'value' has to be a list-like object")
    .match_length(object, value)
    object@mz <- NumericList(lapply(value, "[", , "mz"), compress = FALSE)
    object@intensity <- NumericList(lapply(value, "[", , "intensity"),
                                    compress = FALSE)
    validObject(object)
    object
})
```

Using the `peaksData<-` method we can now also for example remove peaks.

``` r

pd <- peaksData(be)
## Remove the first peak from the first spectrum
pd[[1L]] <- pd[[1L]][-1L, ]

lengths(be)
```

    ## [1] 4 2 3

``` r

peaksData(be) <- pd
lengths(be)
```

    ## [1] 3 2 3

#### `selectSpectraVariables()`

The
[`selectSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
function should allow to reduce the information within the backend
(parameter `object`) to the selected spectra variables (parameter
`spectraVariables`). This is equivalent to a subset by
columns/variables. For core spectra variables, if not specified by
parameter `spectraVariables`, only their values are expected to be
removed (since core spectra variables are expected to be available even
if they are not defined within a backend). The implementation for our
backend will remove any columns in the `@spectraVars` data frame not
defined in the `spectraVariables` parameter. Special care is given to
the `"mz"` and `"intensity"` spectra variables: if they are not
selected, the `@mz` and `@intensity` slots are initialized with empty
`NumericList` (of length matching the number of spectra). Note also that
some backends might throw an error if a spectra variable required for
the backend is removed (such as `"dataStorage"` for a `MsBackendMzR`
backend, which is required by the backend to allow retrieval of m/z and
intensity values).

``` r

setMethod(
    "selectSpectraVariables", "MsBackendTest",
    function(object, spectraVariables = spectraVariables(object)) {
        keep <- colnames(object@spectraVars) %in% spectraVariables
        object@spectraVars <- object@spectraVars[, keep, drop = FALSE]
        if (!any(spectraVariables == "mz"))
            object@mz <- NumericList(vector("list", length(object)),
                                     compress = FALSE)
        if (!any(spectraVariables == "intensity"))
            object@intensity <- NumericList(vector("list", length(object)),
                                            compress = FALSE)
        validObject(object)
        object
    })
```

We can now use
[`selectSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
to remove for example the spectra variable `"new_var"` added above.

``` r

be2 <- be
be2 <- selectSpectraVariables(be2, c("msLevel", "rtime", "mz",
                                     "intensity", "dataStorage"))
spectraVariables(be2)
```

    ##  [1] "msLevel"                 "rtime"                  
    ##  [3] "acquisitionNum"          "scanIndex"              
    ##  [5] "mz"                      "intensity"              
    ##  [7] "dataStorage"             "dataOrigin"             
    ##  [9] "centroided"              "smoothed"               
    ## [11] "polarity"                "precScanNum"            
    ## [13] "precursorMz"             "precursorIntensity"     
    ## [15] "precursorCharge"         "collisionEnergy"        
    ## [17] "isolationWindowLowerMz"  "isolationWindowTargetMz"
    ## [19] "isolationWindowUpperMz"

The spectra variable `"new_var"` is now no longer be available. Note
however that still **all** core spectra variables are listed, even if
they were not selected with the `spectraVariables` parameter. While
these variables (such as `"dataOrigin"`) are still listed by
`spectraVariables(be2)`, their actual values have been removed:

``` r

dataOrigin(be)
```

    ## [1] "<user provided>" "<user provided>" "<user provided>"

``` r

dataOrigin(be2)
```

    ## [1] NA NA NA

If `"mz"` and `"intensitity"` are not selected, the m/z and intensity
values get removed.

``` r

be2 <- selectSpectraVariables(be2, c("msLevel", "rtime", "dataStorage"))
mz(be2)
```

    ## NumericList of length 3
    ## [[1]] numeric(0)
    ## [[2]] numeric(0)
    ## [[3]] numeric(0)

``` r

intensity(be2)
```

    ## NumericList of length 3
    ## [[1]] numeric(0)
    ## [[2]] numeric(0)
    ## [[3]] numeric(0)

#### `dataStorage<-`

Replace the values for the data storage spectra variable. Parameter
`value` has to be of type `character`. Since our backend does not really
make any use of this spectra variable, we can accept any character
value. For other backends, that for example need to load data on-the-fly
from data files, this spectra variable could be used to store the name
of the data files and hence we would need to perform some additional
checks within this replacement function.

``` r

setReplaceMethod("dataStorage", "MsBackendTest", function(object, value) {
    object$dataStorage <- value
    object
})
```

``` r

dataStorage(be)
```

    ## [1] "<memory>" "<memory>" "<memory>"

``` r

dataStorage(be) <- c("", "", "")
dataStorage(be)
```

    ## [1] "" "" ""

#### `spectraNames<-`

Replace the names of individual spectras within the backend. Same as for
`names`, `colnames` or `rownames`, `spectraNames` are expected to be of
type `character`. In our backend implementation we store the spectra
names into the `rownames` of the `@spectraVars` data frame.

``` r

setReplaceMethod("spectraNames", "MsBackendTest", function(object, value) {
    rownames(object@spectraVars) <- value
    object
})
```

``` r

spectraNames(be) <- c("a", "b", "c")
spectraNames(be)
```

    ## [1] "a" "b" "c"

### Optional methods

Default implementations for these methods are available for `MsBackend`
classes, thus these methods don’t have to be implemented for each new
backend. For some backends, depending on how the data is represented or
accessed within it, different implementations might however be more
efficient.

#### `$`

The `$` method is expected to extract a single spectra variable from a
backend. Parameter `name` should allow to name the spectra variable to
return. Each `MsBackend` **must** support extracting the core spectra
variables with this method (even if no data might be available for that
variable). In our example implementation below we make use of the
[`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method, but more efficient implementations might be available as well
(that would not require to first subset/create a `DataFrame` with the
full data and to then subset that again). Also, the `$` method should
check if the requested spectra variable is available and should throw an
error otherwise. The default implementation of the method is shown
below.

``` r

setMethod("$", "MsBackend", function(x, name) {
    spectraData(x, columns = name)[, 1L]
})
```

#### `acquisitionNum()`

Extract the `acquisitionNum` core spectra variable. The method is
expected to return an `integer` vector with the same length as there are
spectra represented by the backend. The default implementation is shown
below.

``` r

setMethod("acquisitionNum", "MsBackend", function(object) {
    spectraData(object, "acquisitionNum")[, 1L]
})
```

#### `backendBpparam()`

The
[`backendBpparam()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method is supposed to evaluate whether a provided (or the default)
parallel processing setup is supported by the backend. Backends that do
not support parallel processing should return
[`SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
instead.

The default implementation is shown below.

``` r

setMethod("backendBpparam", signature = "MsBackend",
          function(object, BPPARAM = bpparam()) {
              ## Return SerialParam() instead to disable parallel processing
              BPPARAM
          })
```

#### `backendParallelFactor()`

The
[`backendParallelFactor()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
allows a backend to suggest a preferred way how the backend could be
split for parallel processing. See also the notes on parallel processing
above for more information. The default implementation returns
[`factor()`](https://rdrr.io/r/base/factor.html) (i.e. a `factor` of
length 0) hence not suggesting any splitting:

``` r

setMethod("backendParallelFactor", "MsBackend", function(object, ...) {
    factor()
})
```

#### `backendRequiredSpectraVariables()`

The
[`backendRequiredSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method can be implemented if a backend needs specific spectra variables
to work. The default implementation is:

``` r

setMethod("backendRequiredSpectraVariables", "MsBackend",
          function(object, ...) {
              character()
          })
```

The implementation for `MsBackendMzR` returns
`c("dataStorage", "scanIndex")` as the backend needs these two spectra
variables to load the MS data on-the-fly from the original data files.

#### `centroided()`

Extract for each spectrum the information whether it contains
*centroided* data. The method is expected to return a `logical` vector
with the same length as there are spectra represented by the backend.
The default implementation is shown below.

``` r

setMethod("centroided", "MsBackend", function(object) {
    spectraData(object, "centroided")[, 1L]
})
```

#### `centroided<-`

Replace the value for the *centroided* core spectra variable. The
provided data type **must** be a logical. The default implementation of
this method is:

``` r

setReplaceMethod("centroided", "MsBackend", function(object, value) {
    object$centroided <- value
    object
})
```

``` r

centroided(be) <- c(TRUE, FALSE, TRUE)
centroided(be)
```

    ## [1]  TRUE FALSE  TRUE

#### `collisionEnergy()`

Extract for each spectrum the collision energy applied to generate the
fragment spectrum. The method is expected to return a `numeric` vector
with the same length as there are spectra represented by the backend
(with `NA_real_` for spectra for which this information is not
available, such as MS1 spectra). The default implementation is shown
below.

``` r

setMethod("collisionEnergy", "MsBackend", function(object) {
    spectraData(object, "collisionEnergy")[, 1L]
})
```

#### `collisionEnergy<-`

Replace the values for the collision energy. Parameter `value` has to be
of type `numeric`. The default implementation of this method is:

``` r

setReplaceMethod("collisionEnergy", "MsBackend", function(object, value) {
    object$collisionEnergy <- value
    object
})
```

#### `dataOrigin()`

Extract the *data origin* spectra variable for each spectrum. This
spectra variable can be used to store the origin of each spectra. The
method is expected to return a `character` vector of length equal to the
number of spectra represented by the backend. The default implementation
is shown below.

``` r

setMethod("dataOrigin", "MsBackend", function(object) {
    spectraData(object, "dataOrigin")[, 1L]
})
```

#### `dataOrigin<-`

Replace the values for the data origin spectra variable. Parameter
`value` has to be of type `character`. The default implementation of
this method is:

``` r

setReplaceMethod("dataOrigin", "MsBackend", function(object, value) {
    object$dataOrigin <- value
    object
})
```

``` r

dataOrigin(be)
```

    ## [1] "<user provided>" "<user provided>" "<user provided>"

``` r

dataOrigin(be) <- c("unknown", "file a", "file b")
dataOrigin(be)
```

    ## [1] "unknown" "file a"  "file b"

#### `dataStorage()`

The `dataStorage` spectra variable of a spectrum provides some
information how or where the data is stored. The
[`dataStorage()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method should therefor return a `character` vector with length equal to
the number of spectra of a backend object with that information. For
most backends the data storage information can be a simple string such
as `"memory"` or `"database"` to specify that the data of a spectrum is
stored within the object itself or in a database, respectively.

Backend classes that keep only a subset of the MS data in memory and
need to load data from data files upon request will use this spectra
variable to store and keep track of the original data file for each
spectrum. An example is the `MsBackendMzR` backend that retrieves the MS
data on-the-fly from the original data file(s) whenever m/z or intensity
values are requested from the backend. Calling
[`dataStorage()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
on an `MsBackendMzR` returns thus the names from the originating files.

The default implementation of the method is shown below.

``` r

setMethod("dataOrigin", "MsBackend", function(object) {
    spectraData(object, "dataStorage")[, 1L]
})
```

#### `dropNaSpectraVariables()`

The
[`dropNaSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
is supposed to allow removing all spectra variables from a data set
(storage) that contain only missing values (i.e. where the value of a
spectra variable for each spectrum is `NA`). This function is intended
to reduce memory requirements of backends such as the `MsBackendMzR`
that load values from all core spectra variables from the original data
files, even if their values are only `NA`. Removing these missing values
from the backend can hence reduce the size in memory of a backend
without data loss (because methods extracting core spectra variables are
supposed to always return `NA` values even if no data is available for
them - in such cases the `NA` values are supposed to be created
on-the-fly.

The default implementation is shown below.

``` r

setMethod("dropNaSpectraVariables", "MsBackend", function(object) {
    svs <- spectraVariables(object)
    svs <- svs[!(svs %in% c("mz", "intensity"))]
    spd <- spectraData(object, columns = svs)
    keep <- !vapply1l(spd, function(z) all(is.na(z)))
    selectSpectraVariables(object, c(svs[keep], "mz", "intensity"))
})
```

#### `isEmpty()`

The
[`isEmpty()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method is expected to return for each spectrum the information whether
it is *empty*, i.e. does not contain any MS peaks (and hence m/z or
intensity values). The result of the method has to be a `logical` of
length equal to the number of spectra represented by the backend with
`TRUE` indicating whether a spectrum is empty and `FALSE` otherwise. The
default implementation is shown below. It requires the
[`lengths()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method to be implemented.

``` r

setMethod("isEmpty", "MsBackend", function(x) {
    lengths(x) == 0L
})
```

#### `isolationWindowLowerMz()`

Extract the core spectra variable `isolationWindowLowerMz` from the
backend. This information is usually provided for each spectrum in the
raw mzML files. The method is expected to return a `numeric` vector of
length equal to the number of spectra represented by the backend. The
default implementation is shown below.

``` r

setMethod("isolationWindowLowerMz", "MsBackend", function(object) {
    spectraData(object, "isolationWindowLowerMz")[, 1L]
})
```

#### `isolationWindowLowerMz<-`

Replace the values for the *isolation window lower m/z* spectra
variable. Parameter `value` has to be of type `numeric` (`NA_real_`
missing values are supported, e.g. for MS1 spectra). The default
implementation of this method is:

``` r

setReplaceMethod(
    "isolationWindowLowerMz", "MsBackend", function(object, value) {
        object$isolationWindowLowerMz <- value
        object
    })
```

``` r

isolationWindowLowerMz(be) <- c(NA_real_, 245.3, NA_real_)
isolationWindowLowerMz(be)
```

    ## [1]    NA 245.3    NA

#### `isolationWindowTargetMz()`

Extract the core spectra variable `isolationWindowTargetMz` from the
backend. This information is usually provided for each spectrum in the
raw mzML files. The method is expected to return a `numeric` vector of
length equal to the number of spectra represented by the backend. The
default implementation is shown below.

``` r

setMethod("isolationWindowTargetMz", "MsBackend", function(object) {
    spectraData(object, "isolationWindowTargetMz")[, 1L]
})
```

#### `isolationWindowTargetMz<-`

Replace the values for the *isolation window target m/z* spectra
variable. Parameter `value` has to be of type `numeric` (`NA_real_`
missing values are supported, e.g. for MS1 spectra). The default
implementation of this method is:

``` r

setReplaceMethod(
    "isolationWindowTargetMz", "MsBackend", function(object, value) {
        object$isolationWindowTargetMz <- value
        object
    })
```

``` r

isolationWindowTargetMz(be) <- c(NA_real_, 245.4, NA_real_)
isolationWindowTargetMz(be)
```

    ## [1]    NA 245.4    NA

#### `isolationWindowUpperMz()`

Extract the core spectra variable `isolationWindowUpperMz` from the
backend. This information is usually provided for each spectrum in the
raw mzML files. The method is expected to return a `numeric` vector of
length equal to the number of spectra represented by the backend. The
default implementation is shown below.

``` r

setMethod("isolationWindowUpperMz", "MsBackend", function(object) {
    spectraData(object, "isolationWindowUpperMz")[, 1L]
})
```

#### `isolationWindowUpperMz<-`

Replace the values for the *isolation window upper m/z* spectra
variable. Parameter `value` has to be of type `numeric` (`NA_real_`
missing values are supported, e.g. for MS1 spectra). The default
implementation of this method is:

``` r

setReplaceMethod(
    "isolationWindowUpperMz", "MsBackend", function(object, value) {
        object$isolationWindowUpperMz <- value
        object
    })
```

``` r

isolationWindowUpperMz(be) <- c(NA_real_, 245.5, NA_real_)
isolationWindowUpperMz(be)
```

    ## [1]    NA 245.5    NA

#### `isReadOnly()`

`isReadOnly()` is expected to return a `logical(1)` with either `TRUE`
or `FALSE` indicating whether the backend supports replacing data or
not. The default implementation is shown below.

``` r

setMethod("isReadOnly", "MsBackend", function(object) {
    object@readonly
})
```

#### `length()`

[`length()`](https://rdrr.io/r/base/length.html) is expected to return a
single `integer` with the total number of spectra that are available
through the backend class. The default implementation is shown below.
Note that the default implementation first accesses the **full** data,
thus it is suggested, for new backends, to implement a more efficient
function.

``` r

setMethod("length", "MsBackend", function(x) {
    nrow(spectraData(x))
})
```

#### `lengths()`

The
[`lengths()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method is expected to return an `integer` vector (same length as the
number of spectra in the backend) with the total number of peaks per
spectrum. The default implementation is shown below.

``` r

setMethod("lengths", "MsBackend", function(x, use.names = FALSE) {
    vapply(peaksData(x), nrow, NA_integer_)
})
```

For our `MsBackendTest` we can simply use the
[`lengths()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method on the m/z or intensity values for that.

``` r

setMethod("lengths", "MsBackendTest", function(x, use.names = FALSE) {
    lengths(x@mz, use.names = use.names)
})
```

And we can now get the peaks count per spectrum:

``` r

lengths(be)
```

    ## [1] 3 2 3

#### `msLevel()`

Extract the MS level for each spectrum in the backend. This method is
expected to return an `integer` of length equal to the number of spectra
represented by the backend. The default implementation is shown below.

``` r

setMethod("msLevel", "MsBackend", function(object) {
    spectraData(object, "msLevel")[, 1L]
})
```

#### `msLevel<-`

Replace the MS level of spectra in a backend. Parameter `value` has to
be of type `integer`. Missing values (`NA_integer_`) are supported.

``` r

setReplaceMethod("msLevel", "MsBackend", function(object, value) {
    object$msLevel <- value
    object
})
```

``` r

msLevel(be)
```

    ## [1] 2 1 2

``` r

msLevel(be) <- c(1L, 1L, 2L)
msLevel(be)
```

    ## [1] 1 1 2

#### `peaksVariables()`

The
[`peaksVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
is expected to return a `character` vector with the names of the *peaks
variables* (i.e. information and properties of individual mass peaks)
available in the backend. The default implementation for `MsBackend`
returns by default `c("mz", "intensity")`. This method should only be
implemented for backends that (eventually) also provide additional peaks
variables. The default implementation is shown below.

``` r

setMethod("peaksVariables", "MsBackend", function(object) {
    c("mz", "intensity")
})
```

#### `polarity()`

Extract the `polarity` core spectra variable for each spectrum in the
backend. This method is expected to return an `integer` of length equal
to the number of spectra represented by the backend. Negative and
positive polarity are expected to be encoded by `0L` and `1L`,
respectively.

``` r

setMethod("polarity", "MsBackend", function(object) {
    spectraData(object, "polarity")[, 1L]
})
```

#### `polarity<-`

Replace the values for the polarity spectra variables. Parameter `value`
has to be of type `integer` and should ideally also use the *standard*
encoding `0L` and `1L` for negative and positive polarity (and
`NA_integer` for missing).

``` r

setReplaceMethod("polarity", "MsBackend", function(object, value) {
    object$polarity <- value
    object
})
```

``` r

polarity(be) <- c(0L, 0L, 0L)
polarity(be)
```

    ## [1] 0 0 0

#### `precScanNum()`

Extract the acquisition number of the precursor for each spectrum. This
method is expected to return an `integer` of length equal to the number
of spectra represented by the backend. For MS1 spectra (or if the
acquisition number of the precursor is not provided) `NA_integer_` has
to be returned. The default implementation is shown below.

``` r

setMethod("precScanNum", "MsBackend", function(object) {
    spectraData(object, "precScanNum")[, 1L]
})
```

#### `precursorCharge()`

Extract the charge of the precursor for each spectrum. This method is
expected to return an `integer` of length equal to the number of spectra
represented by the backend. For MS1 spectra (or if the charge of the
precursor is not provided) `NA_integer_` has to be returned. The default
implementation is shown below.

``` r

setMethod("precursorCharge", "MsBackend", function(object) {
    spectraData(object, "precursorCharge")[, 1L]
})
```

#### `precursorIntensity()`

Extract the intensity of the precursor for each spectrum. This method is
expected to return an `numeric` of length equal to the number of spectra
represented by the backend. For MS1 spectra (or if the precursor
intensity for a fragment spectrum is not provided) `NA_real_` has to be
returned. The default implementation is shown below.

``` r

setMethod("precursorIntensity", "MsBackend", function(object) {
    spectraData(object, "precursorIntensity")[, 1L]
})
```

#### `precursorMz()`

Extract the precursor m/z for each spectrum. This method is expected to
return an `numeric` of length equal to the number of spectra represented
by the backend. For MS1 spectra (or if the precursor m/z for a fragment
spectrum is not provided) `NA_real_` has to be returned. The default
implementation is shown below.

``` r

setMethod("precursorMz", "MsBackend", function(object) {
    spectraData(object, "precursorMz")[, 1L]
})
```

#### `precursorMz<-`

Replace the values for the *precursor m/z* spectra variable. Parameter
`value` has to be of type `numeric` (`NA_real_` missing values are
supported, e.g. for MS1 spectra). The default implementation uses the
`$<-` method:

``` r

setReplaceMethod("precursorMz", "MsBackend", function(object, ..., value) {
    object$precursorMz <- value
    object
})
```

#### `ionCount()`

The
[`ionCount()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method should return a `numeric` (length equal to the number of spectra
represented by the backend) with the sum of all intensities within each
spectrum. For empty spectra `NA_real_` should be returned. The method
below is the default implementation of the method.

``` r

setMethod("ionCount", "MsBackend", function(object) {
    vapply1d(intensity(object), sum, na.rm = TRUE)
})
```

#### `isCentroided()`

This method should return the information for each spectrum whether it
is *centroided*. In contrast to the
[`centroided()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method a heuristic approach is used. The default implementation is shown
below.

``` r

setMethod("isCentroided", "MsBackend", function(object, ...) {
    vapply1l(peaksData(object), .peaks_is_centroided)
})
```

#### `longForm()`

This method returns the MS data of a `MsBackend` as a `data.frame` with
columns being individual spectra or peaks variables and each row
representing one mass peak. A default implementation is available that
converts the `DataFrame` returned by a
[`spectraData()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
call into a *long form* `data.frame` by repeating spectra variable
values to match the number of peaks per spectrum.

#### `reset()`

This is a special method that backends may implement or support, but
don’t necessary have to. This method will be called by the
`reset,Spectra` method and is supposed to restore the data to its
original state. The default implementation for `MsBackend` shown below
simply returns the backend as-is. The `MsBackendSql` backend from the
*[MsBackendSql](https://bioconductor.org/packages/3.23/MsBackendSql)*
package in contrast re-initializes the data using the data from the
database.

``` r

setMethod("reset", "MsBackend", function(object) {
    object
})
```

#### `export()`

This method should *export* the data from a `MsBackend`. The method is
called by the `export,Spectra` method that passes itself as a second
argument to the function. The `export,MsBackend` implementation is thus
expected to take a `Spectra` object as second argument from which all
data should be taken and exported. Implementation of this method is
optional. The implementation of the method for the `MsBackendMzR`
backend is shown below.

``` r

setMethod("export", "MsBackendMzR", function(object, x, file = tempfile(),
                                             format = c("mzML", "mzXML"),
                                             copy = FALSE,
                                             BPPARAM = bpparam()) {
    l <- length(x)
    file <- sanitize_file_name(file)
    if (length(file) == 1)
        file <- rep_len(file, l)
    if (length(file) != l)
        stop("Parameter 'file' has to be either of length 1 or ",
             length(x), ", i.e. 'length(x)'.", call. = FALSE)
    f <- factor(file, levels = unique(file))
    tmp <- bpmapply(.write_ms_data_mzR, split(x, f), levels(f),
                    MoreArgs = list(format = format, copy = copy),
                    BPPARAM = BPPARAM)
})
```

See alternatively also the
*[MsBackendMgf](https://bioconductor.org/packages/3.23/MsBackendMgf)*
package for an implementation for the `MsBackendMgf` backend.

#### `rtime()`

Extract the retention time of each spectrum. This method is expected to
return a `numeric` of length equal to the number of spectra represented
by the backend. the default implementation is shown below.

``` r

setMethod("rtime", "MsBackend", function(object) {
    spectraData(object, "rtime")[, 1L]
})
```

#### `rtime<-`

Replace the retention times for the spectra represented by the backend.
Parameter `value` must be of type `numeric`. Also, although it is not a
strict requirement, retention times should ideally be ordered
increasingly *per sample* and their unit should be seconds. The default
implementation is:

``` r

setReplaceMethod("rtime", "MsBackend", function(object, value) {
    object$rtime <- value
    object
})
```

``` r

rtime(be)
```

    ## [1] 1.2 1.3 1.4

``` r

rtime(be) <- rtime(be) + 2
rtime(be)
```

    ## [1] 3.2 3.3 3.4

#### `scanIndex()`

Extract the *scan index* core spectra variable. The scan index
represents the relative index of the spectrum within the respective raw
data file and can be different than the `acquisitionNum` (which is the
index of a spectrum as recorded by the MS instrument). This method is
expected to return a `integer` of length equal to the number of spectra
represented by the backend. The default implementation is shown below.

``` r

setMethod("scanIndex", "MsBackend", function(object) {
    spectraData(object, "scanIndex")[, 1L]
})
```

#### `smoothed()`

Extract the `smoothed` core spectra variable that indicates whether a
spectrum was *smoothed*. This variable is supported for backward
compatibility but seldomly used. The method is expected to return a
`logical` with length equal to the number of spectra represented by the
backend. The default implementation is shown below.

``` r

setMethod("smoothed", "MsBackend", function(object) {
    spectraData(object, "smoothed")[, 1L]
})
```

#### `smoothed<-`

Replace the spectra variable *smoothed* that indicates whether some data
smoothing operation was performed on the spectra. Parameter `value` must
be of type `logical`. The default implementation is:

``` r

setReplaceMethod("smoothed", "MsBackend", function(object, value) {
    object$smoothed <- value
    object
})
```

``` r

smoothed(be) <- rep(TRUE, 3)
smoothed(be)
```

    ## [1] TRUE TRUE TRUE

#### `split()`

The
[`split()`](https://rformassspectrometry.github.io/Spectra/reference/combineSpectra.md)
method should allow to split a `MsBackend` into a `list` of `MsBackend`
objects. The default implementation is shown below.

``` r

setMethod("split", "MsBackend", function(x, f, drop = FALSE, ...) {
    split.default(x, f, drop = drop, ...)
})
```

#### `supportsSetBackend()`

Whether a `MsBackend` supports the
[`setBackend()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
method that allows to change the backend of a `Spectra` object from one
to another backend. To support
[`setBackend()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
the backend needs to have a parameter `data` in its
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
method that allows it to be initialized with a `DataFrame` containing
the full spectra data. The default implementation is shown below.

``` r

setMethod("supportsSetBackend", "MsBackend", function(object, ...) {
    !isReadOnly(object)
})
```

#### `tic()`

The
[`tic()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
method should return the total ion count (i.e. the sum of intensities)
for each spectrum. This information is usually also provided by the raw
MS data files, but can also be calculated on the fly from the data. The
parameter `initial` (which is by default `TRUE`) allows to define
whether the provided *original* tic should be returned (for
`initial = TRUE`) or whether the tic should be calculated on the actual
data (`initial = FALSE`). The original tic values are usually provided
by a spectra variable `"totIonCurrent"`. The default implementation is
shown below.

``` r

setMethod("tic", "MsBackend", function(object, initial = TRUE) {
    if (initial) {
        if (any(spectraVariables(object) == "totIonCurrent"))
            spectraData(object, "totIonCurrent")[, 1L]
        else rep(NA_real_, length(object))
    } else vapply(intensity(object), sum, numeric(1), na.rm = TRUE)
})
```

#### `uniqueMsLevels()`

This method should return the unique MS level(s) of all spectra within
the backend. The default implementation is shown below.

``` r

setMethod("uniqueMsLevels", "MsBackend", function(object, ...) {
    unique(msLevel(object))
})
```

This method thus retrieves first the MS levels of all spectra and then
calls [`unique()`](https://rdrr.io/pkg/BiocGenerics/man/unique.html) on
them. Database-based backends might avoid such an eventually heavy
operation by selecting the unique MS levels directly using an SQL call.

#### `filterDataOrigin()`

The `filter*` methods are expected to take a `MsBackend` and to subset
it based on some criteria. While also the `[` method could be used to
perform such subset operation, these methods might allow more efficient
ways to subset the data e.g. by performing the operation within a
database with a dedicated SQL call. A default implementation is
available for every filter function and thus a method needs only to be
implemented if the data storage/representation within a backend would
allow a more efficient operation.

All filter methods are expected to return the subset backend (i.e. an
instance of the same backend class with the same, or less spectra).

The
[`filterDataOrigin()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
should subset the backend to spectra with their `dataOrigin` spectra
variable matching the values provided with the `dataOrigin` parameter.
The default implementation for `MsBackend` is shown below.

``` r

setMethod("filterDataOrigin", "MsBackend",
          function(object, dataOrigin = character()) {
              if (length(dataOrigin)) {
                  object <- object[dataOrigin(object) %in% dataOrigin]
                  if (is.unsorted(dataOrigin))
                      object[order(match(dataOrigin(object), dataOrigin))]
                  else object
              } else object
          })
```

#### `filterDataStorage()`

Similar to the
[`filterDataOrigin()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md),
the
[`filterDataStorage()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
should subset a backend to spectra with their `dataStorage` spectra
variable matching the values provided with the `dataStorage` parameter.

``` r

setMethod("filterDataStorage", "MsBackend",
          function(object, dataStorage = character()) {
              if (length(dataStorage)) {
                  object <- object[dataStorage(object) %in% dataStorage]
                  if (is.unsorted(dataStorage))
                      object[order(match(dataStorage(object), dataStorage))]
                  else object
              } else object
          })
```

#### `filterEmptySpectra()`

The
[`filterEmptySpectra()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
should remove all *empty* spectra (i.e. spectra without any mass peaks)
from the backend. The method is expected to return the subset backend.
The default implementation for `MsBackend` is shown below.

``` r

setMethod("filterEmptySpectra", "MsBackend", function(object, ...) {
    if (!length(object)) return(object)
    object[as.logical(lengths(object))]
})
```

#### `filterIsolationWindow()`

The
[`filterIsolationWindow()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
filters the backend to spectra with the provided `mz` value being within
their
[`isolationWindowLowerMz()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
and
[`isolationWindowUpperMz()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md).
The parameter `mz` defining this target m/z is expected to be a
`numeric` of length 1. The default implementation for `MsBackend` is
shown below.

``` r

setMethod("filterIsolationWindow", "MsBackend",
          function(object, mz = numeric(), ...) {
              if (length(mz)) {
                  if (length(mz) > 1)
                      stop("'mz' is expected to be a single m/z value")
                  keep <- which(isolationWindowLowerMz(object) <= mz &
                                isolationWindowUpperMz(object) >= mz)
                  object[keep]
              } else object
          })
```

#### `filterMsLevel()`

The
[`filterMsLevel()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
method is expected to reduce the backend to spectra with the provided MS
level(s). Parameter `msLevel` has to be an `integer` (any length). The
default implementation for `MsBackend` is shown below.

``` r

setMethod("filterMsLevel", "MsBackend",
          function(object, msLevel = integer()) {
              if (length(msLevel)) {
                  object[msLevel(object) %in% msLevel]
              } else object
          })
```

#### `filterPolarity()`

The
[`filterPolarity()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
method is expected to subset the backend to spectra matching the
provided polarity (or polarities). Parameter `polarity` has to be an
`integer` (of any length). The default implementation for `MsBackend` is
shown below.

``` r

setMethod("filterPolarity", "MsBackend",
          function(object, polarity = integer()) {
              if (length(polarity))
                  object[polarity(object) %in% polarity]
              else object
          })
```

#### `filterPrecursorMzRange()`

The
[`filterPrecursorMzRange()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
method filters the backend to spectra with their `precursorMz` being
between the provided m/z range (parameter `mz`). This method was
previously named `filterPrecursorMz()`. Parameter `mz` is expected to be
a `numeric` of length 2 defining the lower and upper limit of this
precursor m/z range. The default implementation for `MsBackend` is shown
below.

``` r

library(MsCoreUtils)
```

    ## 
    ## Attaching package: 'MsCoreUtils'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     reduce

    ## The following objects are masked from 'package:Spectra':
    ## 
    ##     bin, entropy, smooth

    ## The following object is masked from 'package:stats':
    ## 
    ##     smooth

``` r

setMethod("filterPrecursorMzRange", "MsBackend",
          function(object, mz = numeric()) {
              if (length(mz)) {
                  mz <- range(mz)
                  keep <- which(between(precursorMz(object), mz))
                  object[keep]
              } else object
          })
```

#### `filterPrecursorMzValues()`

The
[`filterPrecursorMzValues()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
method filters the backend to spectra with their m/z values matching to
the provided m/z value(s). Parameters `ppm` and `tolerance` (both
expected to be `numeric` of length 1) allow to define the conditions for
the relaxed matching. Parameter `mz` has to be a `numeric` (of any
length). The default implementation for `MsBackend` is shown below.

``` r

setMethod("filterPrecursorMzValues", "MsBackend",
          function(object, mz = numeric(), ppm = 20, tolerance = 0) {
              if (length(mz)) {
                  object[.values_match_mz(precursorMz(object), mz = mz,
                                          ppm = ppm, tolerance = tolerance)]
              } else object
          })
```

The `.values_match_mz` function used by this function is defined as:

``` r

.values_match_mz <- function(x, mz, ppm = 20, tolerance = 0) {
    o <- order(x, na.last = NA)
    cmn <- common(x[o], sort(mz), tolerance = tolerance, ppm = ppm,
                  duplicates = "keep", .check = FALSE)
    sort(o[cmn])
}
```

#### `filterPrecursorCharge()`

The
[`filterPrecursorCharge()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
method filters the backend to spectra with matching precursor charge.
Parameter `z` defining the requested precursor charge has to be an
`integer` (of any length). The default implementation for `MsBackend` is
shown below.

``` r

setMethod("filterPrecursorCharge", "MsBackend",
          function(object, z = integer()) {
              if (length(z)) {
                  keep <- which(precursorCharge(object) %in% z)
                  object[keep]
              } else object
          })
```

#### `filterPrecursorScan()`

The
[`filterPrecursorScan()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
method filters the backend to parent (e.g. MS1) and children scans
(e.g. MS2) of acquisition number `acquisitionNum`. Parameter `f` defines
how the backend should be split (by default by original data file) to
avoid selecting spectra from different samples/files. The default
implementation for `MsBackend` is shown below.

``` r

setMethod("filterPrecursorScan", "MsBackend",
          function(object, acquisitionNum = integer(), f = dataOrigin(object)) {
              if (length(acquisitionNum) && length(f)) {
                  if (!is.factor(f))
                      f <- factor(f, exclude = character())
                  keep <- unsplit(lapply(split(object, f = f), function(z, an) {
                      .filterSpectraHierarchy(acquisitionNum(z),
                                              precScanNum(z),
                                              an)
                  }, an = acquisitionNum), f = f)
                  object[keep]
              } else object
          })
```

#### `filterRt()`

The
[`filterRt()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
method filters the backend to spectra with their retention time being
between the provided rt range. Parameter `rt` is expected to be a
`numeric` of length 2 defining the lower and upper bound of this range.
Parameter `msLevel.` (note the `.` in the name of the parameter!) can be
optionally used to restrict the filtering to the selected MS levels
(i.e. the RT filter is only applied to spectra of the selected MS levels
and all spectra with a different MS level are returned as well). The
default implementation for `MsBackend` is shown below.

``` r

setMethod("filterRt", "MsBackend",
          function(object, rt = numeric(), msLevel. = integer()) {
              if (length(rt)) {
                  rt <- range(rt)
                  sel_ms <- msLevel(object) %in% msLevel.
                  sel_rt <- between(rtime(object), rt) & sel_ms
                  object[sel_rt | !sel_ms]
              } else object
          })
```

### Implementation notes

In this tutorial we implemented a simple *in-memory* `MsBackend` from
scratch. For many real-life situation it might however be better to
extend some of the pre-defined backend classes from the `Spectra`
package to avoid duplicating functionality. A good starting point might
be the `MsBackendMemory` backend for any *in-memory* data
representation, or the `MsBackendCached` for backends that retrieve data
from inherently read-only resources (such as database connection or raw
data files) but still would need to support adding spectra variables or
changing values of spectra variables. Similarly, if the only purpose of
a backend is to import or export data in a specific format, the
`MsBackendMemory` might be extended and a single method
([`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md))
would need to be implemented for the new class: this new
[`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
would then call the code to import the data from the new file format and
store it within the available slots of the `MsBackendMemory` object.
Examples would be the backends provided by the
*[MsBackendMgf](https://bioconductor.org/packages/3.23/MsBackendMgf)*
and
*[MsBackendMsp](https://bioconductor.org/packages/3.23/MsBackendMsp)*
classes.

## Testing the validity of the backend

The `Spectra` package provides a set of unit tests that allow to check a
backend for compliance with `MsBackend`. Below we load this test suite
and call the tests. The tests will be performed on a variable `be` in
the current workspace (which in our case is an instance of our
`MsBackendTest` class).

``` r

library(testthat)
test_suite <- system.file("test_backends", "test_MsBackend",
                          package = "Spectra")
test_dir(test_suite, stop_on_failure = TRUE)
```

## Session information

``` r

sessionInfo()
```

    ## R Under development (unstable) (2025-12-21 r89216)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] MsCoreUtils_1.23.2  IRanges_2.45.0      Spectra_1.21.1     
    ## [4] BiocParallel_1.45.0 S4Vectors_0.49.0    BiocGenerics_0.57.0
    ## [7] generics_0.1.4      BiocStyle_2.39.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] jsonlite_2.0.0         compiler_4.6.0         BiocManager_1.30.27   
    ##  [4] parallel_4.6.0         cluster_2.1.8.1        jquerylib_0.1.4       
    ##  [7] systemfonts_1.3.1      textshaping_1.0.4      yaml_2.3.12           
    ## [10] fastmap_1.2.0          R6_2.6.1               ProtGenerics_1.43.0   
    ## [13] knitr_1.51             htmlwidgets_1.6.4      MASS_7.3-65           
    ## [16] bookdown_0.46          desc_1.4.3             bslib_0.9.0           
    ## [19] rlang_1.1.6            cachem_1.1.0           xfun_0.55             
    ## [22] fs_1.6.6               sass_0.4.10            otel_0.2.0            
    ## [25] cli_3.6.5              pkgdown_2.2.0.9000     digest_0.6.39         
    ## [28] MetaboCoreUtils_1.19.1 lifecycle_1.0.4        clue_0.3-66           
    ## [31] evaluate_1.0.5         codetools_0.2-20       ragg_1.5.0            
    ## [34] rmarkdown_2.30         tools_4.6.0            htmltools_0.5.9

## References

Gatto, Laurent, Sebastian Gibb, and Johannes Rainer. 2020. “MSnbase,
Efficient and Elegant R-Based Processing and Visualization of Raw Mass
Spectrometry Data.” *Journal of Proteome Research*, ahead of print,
September. <https://doi.org/10.1021/acs.jproteome.0c00313>.
