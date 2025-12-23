# Count the number of identifications per scan

The function takes a `Spectra` object containing identification results
as input. It then counts the number of identifications each scan (or
their descendants) has lead to - this is either 0 or 1 for MS2 scans,
or, for MS1 scans, the number of MS2 scans originating from any MS1 peak
that lead to an identification.

This function can be used to generate id-annotated total ion
chromatograms, as can illustrated
[here](https://rformassspectrometry.github.io/docs/sec-id.html#an-identification-annotated-chromatogram).

## Usage

``` r
countIdentifications(
  object,
  identification = "sequence",
  f = dataStorage(object),
  BPPARAM = bpparam()
)
```

## Arguments

- object:

  An instance of class `Spectra` that contains identification data, as
  defined by the `sequence` argument.

- identification:

  `character(1)` with the name of the spectra variable that defines
  whether a scan lead to an identification (typically containing the
  idenfified peptides sequence in proteomics). The absence of
  identification is encode by an `NA`. Default is `"sequence"`.

- f:

  A `factor` defining how to split `object` for parallelized processing.
  Default is `dataOrigin(x)`, i.e. each raw data files is processed in
  parallel.

- BPPARAM:

  Parallel setup configuration. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for details.

## Value

An updated
[`Spectra()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
object that now contains an integer spectra variable
`countIdentifications` with the number of identification for each scan.

## Details

The computed number of identifications is stored in a new spectra
variables named `"countIdentifications"`. If it already exists, the
function throws a message and returns the object unchanged. To force the
recomputation of the `"countIdentifications"` variable, users should
either delete or rename it.

## See also

[`addProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
for other data analysis functions.

## Author

Laurent Gatto

## Examples

``` r
spdf <- new("DFrame", rownames = NULL, nrows = 86L,
   listData = list(
       msLevel = c(1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                   2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L,
                   2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                   2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                   2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                   2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L,
                   2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                   2L, 2L),
       acquisitionNum = 8975:9060,
       precScanNum = c(NA, 8956L, 8956L, 8956L, 8956L, 8956L, 8956L,
                       8956L, 8956L, 8956L, 8956L, 8956L, 8956L,
                       8956L, 8956L, 8956L, 8956L, 8956L, 8956L, NA,
                       8975L, 8975L, 8975L, 8975L, 8975L, 8975L,
                       8975L, 8975L, 8975L, 8975L, 8975L, 8975L,
                       8975L, 8975L, 8975L, 8975L, 8975L, NA, 8994L,
                       8994L, 8994L, 8994L, 8994L, 8994L, 8994L,
                       8994L, 8994L, 8994L, 8994L, 8994L, 8994L, NA,
                       9012L, 9012L, 9012L, 9012L, 9012L, 9012L,
                       9012L, 9012L, 9012L, 9012L, 9012L, 9012L,
                       9012L, 9012L, 9012L, 9012L, 9012L, 9012L, NA,
                       9026L, 9026L, 9026L, 9026L, 9026L, 9026L,
                       9026L, 9026L, 9026L, 9026L, 9026L, 9026L,
                       9026L, 9026L, 9026L),
       sequence = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                    "LSEHATAPTR", NA, NA, NA, NA, NA, NA, NA,
                    "EGSDATGDGTK", NA, NA, "NEDEDSPNK", NA, NA, NA,
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                    NA, NA, NA, NA, NA, NA, NA, NA, NA, "GLTLAQGGVK",
                    NA, NA, NA, NA, "STLPDADRER", NA, NA, NA, NA, NA,
                    NA, NA, NA)),
   elementType = "ANY", elementMetadata = NULL, metadata = list())

sp <- Spectra(spdf)

## We have in this data 5 MS1 and 81 MS2 scans
table(msLevel(sp))
#> 
#>  1  2 
#>  5 81 

## The acquisition number of the MS1 scans
acquisitionNum(filterMsLevel(sp, 1))
#> [1] 8975 8994 9012 9026 9045

## And the number of MS2 scans with precursor ions selected
## from MS1 scans (those in the data and others)
table(precScanNum(sp))
#> 
#> 8956 8975 8994 9012 9026 
#>   18   17   13   18   15 

## Count number of sequences/identifications per scan
sp <- countIdentifications(sp)

## MS2 scans either lead to an identification (5 instances) or none
## (76). Among the five MS1 scans in the experiment, 3 lead to MS2
## scans being matched to no peptides and two MS1 scans produced two
## and three PSMs respectively.
table(sp$countIdentifications, sp$msLevel)
#>    
#>      1  2
#>   0  3 76
#>   1  0  5
#>   2  1  0
#>   3  1  0
```
