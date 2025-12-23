# Mapping between spectra variables and data file fields

The `spectraVariableMapping` function provides the mapping between
*spectra variables* of a
[`Spectra()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
object with data fields from a data file. Such name mapping is expected
to enable an easier import of data files with specific *dialects*, e.g.
files in MGF format that use a different naming convention for core
spectra variables.

[`MsBackend()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
implementations are expected to implement this function (if needed) to
enable import of data from file formats with non-standardized data
fields.

## Usage

``` r
spectraVariableMapping(object, ...)

spectraVariableMapping(object, ...) <- value

# S4 method for class 'MsBackend'
spectraVariableMapping(object)

# S4 method for class 'MsBackend'
spectraVariableMapping(object) <- value

# S4 method for class 'Spectra'
spectraVariableMapping(object) <- value

# S4 method for class 'Spectra'
spectraVariableMapping(object)
```

## Arguments

- object:

  An instance of an object extending
  [`MsBackend()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md).

- ...:

  Optional parameters.

- value:

  For `spectraVariableMapping<-`: a named `character` vector.

## Value

A named `character` with names being spectra variable names (use
[`spectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
for a list of supported names) and values being the data field names.

## Author

Johannes Rainer
