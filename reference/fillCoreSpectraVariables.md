# Fill spectra data with columns for missing core variables

`fillCoreSpectraVariables()` fills a provided `data.frame` with columns
for eventually missing *core* spectra variables. The missing core
variables are added as new columns with missing values (`NA`) of the
correct data type. Use
[`coreSpectraVariables()`](https://rformassspectrometry.github.io/Spectra/reference/spectraData.md)
to list the set of core variables and their data types.

## Usage

``` r
fillCoreSpectraVariables(
  x = data.frame(),
  columns = names(coreSpectraVariables())
)
```

## Arguments

- x:

  `data.frame` or `DataFrame` with potentially present core variable
  columns.

- columns:

  `character` with the names of the (core) spectra variables that should
  be added if not already present in `x`. Defaults to
  `columns = names(coreSpectraVariables())`.

## Value

input data frame `x` with missing core variables added (with the correct
data type).

## Examples

``` r

## Define a data frame
a <- data.frame(msLevel = c(1L, 1L, 2L), other_column = "b")

## Add missing core chromatogram variables to this data frame
fillCoreSpectraVariables(a)
#> Warning: corrupt data frame: columns will be truncated or padded with NAs
#>   msLevel other_column rtime acquisitionNum scanIndex dataStorage dataOrigin
#> 1       1            b    NA             NA        NA        <NA>       <NA>
#> 2       1            b    NA             NA        NA        <NA>       <NA>
#> 3       2            b    NA             NA        NA        <NA>       <NA>
#>   centroided smoothed polarity precScanNum precursorMz precursorIntensity
#> 1         NA       NA       NA          NA          NA                 NA
#> 2         NA       NA       NA          NA          NA                 NA
#> 3         NA       NA       NA          NA          NA                 NA
#>   precursorCharge collisionEnergy isolationWindowLowerMz
#> 1              NA              NA                     NA
#> 2              NA              NA                     NA
#> 3              NA              NA                     NA
#>   isolationWindowTargetMz isolationWindowUpperMz
#> 1                      NA                     NA
#> 2                      NA                     NA
#> 3                      NA                     NA
#>                                                                mz
#> 1 <S4 class ‘SimpleNumericList’ [package “IRanges”] with 4 slots>
#> 2                                                            <NA>
#> 3                                                            <NA>
#>                                                         intensity
#> 1 <S4 class ‘SimpleNumericList’ [package “IRanges”] with 4 slots>
#> 2                                                            <NA>
#> 3                                                            <NA>

## The data.frame thus contains columns for all core spectra
## variables in the respective expected data type (but filled with
## missing values).
```
