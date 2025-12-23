# Merging, aggregating and splitting Spectra

Various functions are availabe to combine, aggregate or split data from
one of more `Spectra` objects. These are:

- [`c()`](https://rdrr.io/r/base/c.html) and `concatenateSpectra()`:
  combines several `Spectra` objects into a single object. The resulting
  `Spectra` contains all data from all individual `Spectra`, i.e. the
  union of all their spectra variables. Concatenation will fail if the
  processing queue of any of the `Spectra` objects is not empty or if
  different backends are used for the `Spectra` objects. In such cases
  it is suggested to first change the backends of all `Spectra` to the
  same type of backend (using the
  [`setBackend()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  function and to eventually (if needed) apply the processing queue
  using the
  [`applyProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
  function.

- `cbind2()`: Appends multiple spectra variables from a `data.frame`,
  `DataFrame` or `matrix` to the `Spectra` object at once. The order of
  the values (rows) in `y` has to match the order of spectra in `x`. The
  function does not allow to replace existing spectra variables.
  `cbind2()` returns a `Spectra` object with the appended spectra
  variables. For a more controlled way of adding spectra variables, see
  the `joinSpectraData()` function.

- `combineSpectra()`: combines sets of spectra (defined with parameter
  `f`) into a single spectrum per set aggregating their MS data (i.e.
  their *peaks data* matrices with the *m/z* and intensity values of
  their mass peaks). The spectra variable values of the first spectrum
  per set are reported for the combined spectrum. The peak matrices of
  the spectra per set are combined using the function specified with
  parameter `FUN` which uses by default the
  [`combinePeaksData()`](https://rformassspectrometry.github.io/Spectra/reference/combinePeaksData.md)
  function. See the documentation of
  [`combinePeaksData()`](https://rformassspectrometry.github.io/Spectra/reference/combinePeaksData.md)
  for details on the aggregation of the peak data and the package
  vignette for examples. The sets of spectra can be specified with
  parameter `f` which is expected to be a `factor` or `vector` of length
  equal to the length of the `Spectra` specifying to which set a
  spectrum belongs to. The function returns a `Spectra` of length equal
  to the unique levels of `f`. The optional parameter `p` allows to
  define how the `Spectra` should be split for potential parallel
  processing. The default is `p = x$dataStorage` and hence a per storage
  file parallel processing is applied for `Spectra` with on disk data
  representations (such as the
  [`MsBackendMzR()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)).
  This also prevents that spectra from different data files/samples are
  combined (eventually use e.g. `p = x$dataOrigin` or any other spectra
  variables defining the originating samples for a spectrum). Before
  combining the peaks data, all eventual present processing steps are
  applied (by calling
  [`applyProcessing()`](https://rformassspectrometry.github.io/Spectra/reference/addProcessing.md)
  on the `Spectra`). This function will replace the original *m/z* and
  intensity values of a `Spectra` hence it can not be called on a
  `Spectra` with a *read-only* backend. In such cases, the backend
  should be changed to a *writeable* backend before using the
  [`setBackend()`](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  function (to e.g. a
  [`MsBackendMemory()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  backend).

- `joinSpectraData()`: Individual spectra variables can be directly
  added with the `$<-` or `[[<-` syntax. The `joinSpectraData()`
  function allows to merge a `DataFrame` to the existing spectra data of
  a `Spectra`. This function diverges from the
  [`merge()`](https://rdrr.io/r/base/merge.html) method in two main
  ways:

  - The `by.x` and `by.y` column names must be of length 1.

  - If variable names are shared in `x` and `y`, the spectra variables
    of `x` are not modified. It's only the `y` variables that are
    appended with the suffix defined in `suffix.y`. This is to avoid
    modifying any core spectra variables that would lead to an invalid
    object.

  - Duplicated Spectra keys (i.e. `x[[by.x]]`) are not allowed.
    Duplicated keys in the `DataFrame` (i.e `y[[by.y]]`) throw a warning
    and only the last occurrence is kept. These should be explored and
    ideally be removed using for `QFeatures::reduceDataFrame()`,
    `PMS::reducePSMs()` or similar functions. For a more general
    function that allows to append `data.frame`, `DataFrame` and
    `matrix` see `cbind2()`.

- `split()`: splits the `Spectra` object based on parameter `f` into a
  `list` of `Spectra` objects.

## Usage

``` r
concatenateSpectra(x, ...)

combineSpectra(
  x,
  f = x$dataStorage,
  p = x$dataStorage,
  FUN = combinePeaksData,
  ...,
  BPPARAM = bpparam()
)

joinSpectraData(x, y, by.x = "spectrumId", by.y, suffix.y = ".y")

# S4 method for class 'Spectra'
c(x, ...)

# S4 method for class 'Spectra,dataframeOrDataFrameOrmatrix'
cbind2(x, y, ...)

# S4 method for class 'Spectra,ANY'
split(x, f, drop = FALSE, ...)
```

## Arguments

- x:

  A `Spectra` object.

- ...:

  Additional arguments.

- f:

  For `split()`: factor defining how to split `x`. See
  [`base::split()`](https://rdrr.io/r/base/split.html) for details. For
  `combineSpectra()`: `factor` defining the grouping of the spectra that
  should be combined. Defaults to `x$dataStorage`.

- p:

  For `combineSpectra()`: `factor` defining how to split the input
  `Spectra` for parallel processing. Defaults to `x$dataStorage`, i.e.,
  depending on the used backend, per-file parallel processing will be
  performed.

- FUN:

  For `combineSpectra()`: function to combine the (peak matrices) of the
  spectra. Defaults to
  [`combinePeaksData()`](https://rformassspectrometry.github.io/Spectra/reference/combinePeaksData.md).

- BPPARAM:

  Parallel setup configuration. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information. This is passed directly to the
  [`backendInitialize()`](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md)
  method of the
  [MsBackend](https://rformassspectrometry.github.io/Spectra/reference/MsBackend.md).

- y:

  For `joinSpectraData()`: `DataFrame` with the spectra variables to
  join/add. For `cbind2()`: a `data.frame`, `DataFrame` or `matrix`. The
  number of rows and their order has to match the number of spectra in
  `x`, respectively their order.

- by.x:

  A `character(1)` specifying the spectra variable used for merging.
  Default is `"spectrumId"`.

- by.y:

  A `character(1)` specifying the column used for merging. Set to `by.x`
  if missing.

- suffix.y:

  A `character(1)` specifying the suffix to be used for making the names
  of columns in the merged spectra variables unique. This suffix will be
  used to amend `names(y)`, while `spectraVariables(x)` will remain
  unchanged.

- drop:

  For `split()`: not considered.

## See also

- [`combinePeaks()`](https://rformassspectrometry.github.io/Spectra/reference/combinePeaks.md)
  for functions to aggregate mass peaks data.

- [Spectra](https://rformassspectrometry.github.io/Spectra/reference/Spectra.md)
  for a general description of the `Spectra` object.

## Author

Sebastian Gibb, Johannes Rainer, Laurent Gatto

## Examples

``` r

## Create a Spectra providing a `DataFrame` containing a MS data.

spd <- DataFrame(msLevel = c(1L, 2L), rtime = c(1.1, 1.2))
spd$mz <- list(c(100, 103.2, 104.3, 106.5), c(45.6, 120.4, 190.2))
spd$intensity <- list(c(200, 400, 34.2, 17), c(12.3, 15.2, 6.8))

s <- Spectra(spd)
s
#> MSn data (Spectra) with 2 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1       1.1        NA
#> 2         2       1.2        NA
#>  ... 16 more variables/columns.

## Create a second Spectra from mzML files and use the `MsBackendMzR`
## on-disk backend.
sciex_file <- dir(system.file("sciex", package = "msdata"),
    full.names = TRUE)
sciex <- Spectra(sciex_file, backend = MsBackendMzR())
sciex
#> MSn data (Spectra) with 1862 spectra in a MsBackendMzR backend:
#>        msLevel     rtime scanIndex
#>      <integer> <numeric> <integer>
#> 1            1     0.280         1
#> 2            1     0.559         2
#> 3            1     0.838         3
#> 4            1     1.117         4
#> 5            1     1.396         5
#> ...        ...       ...       ...
#> 1858         1   258.636       927
#> 1859         1   258.915       928
#> 1860         1   259.194       929
#> 1861         1   259.473       930
#> 1862         1   259.752       931
#>  ... 34 more variables/columns.
#> 
#> file(s):
#> 20171016_POOL_POS_1_105-134.mzML
#> 20171016_POOL_POS_3_105-134.mzML

## Subset to the first 100 spectra to reduce running time of the examples
sciex <- sciex[1:100]


##  --------  COMBINE SPECTRA  --------

## Combining the `Spectra` object `s` with the MS data from `sciex`.
## Calling directly `c(s, sciex)` would result in an error because
## both backends use a different backend. We thus have to first change
## the backends to the same backend. We change the backend of the `sciex`
## `Spectra` to a `MsBackendMemory`, the backend used by `s`.

sciex <- setBackend(sciex, MsBackendMemory())

## Combine the two `Spectra`
all <- c(s, sciex)
all
#> MSn data (Spectra) with 102 spectra in a MsBackendMemory backend:
#>       msLevel     rtime scanIndex
#>     <integer> <numeric> <integer>
#> 1           1     1.100        NA
#> 2           2     1.200        NA
#> 3           1     0.280         1
#> 4           1     0.559         2
#> 5           1     0.838         3
#> ...       ...       ...       ...
#> 98          1    26.786        96
#> 99          1    27.065        97
#> 100         1    27.344        98
#> 101         1    27.623        99
#> 102         1    27.902       100
#>  ... 34 more variables/columns.
#> Processing:
#>  Switch backend from MsBackendMzR to MsBackendMemory [Tue Dec 23 13:20:26 2025]
#>  Merge 2 Spectra into one [Tue Dec 23 13:20:26 2025] 

## The new `Spectra` objects contains the union of spectra variables from
## both:
spectraVariables(all)
#>  [1] "msLevel"                  "rtime"                   
#>  [3] "acquisitionNum"           "scanIndex"               
#>  [5] "dataStorage"              "dataOrigin"              
#>  [7] "centroided"               "smoothed"                
#>  [9] "polarity"                 "precScanNum"             
#> [11] "precursorMz"              "precursorIntensity"      
#> [13] "precursorCharge"          "collisionEnergy"         
#> [15] "isolationWindowLowerMz"   "isolationWindowTargetMz" 
#> [17] "isolationWindowUpperMz"   "peaksCount"              
#> [19] "totIonCurrent"            "basePeakMZ"              
#> [21] "basePeakIntensity"        "electronBeamEnergy"      
#> [23] "ionisationEnergy"         "lowMZ"                   
#> [25] "highMZ"                   "mergedScan"              
#> [27] "mergedResultScanNum"      "mergedResultStartScanNum"
#> [29] "mergedResultEndScanNum"   "injectionTime"           
#> [31] "filterString"             "spectrumId"              
#> [33] "ionMobilityDriftTime"     "scanWindowLowerLimit"    
#> [35] "scanWindowUpperLimit"    

## The spectra variables that were not present in `s`:
setdiff(spectraVariables(all), spectraVariables(s))
#>  [1] "peaksCount"               "totIonCurrent"           
#>  [3] "basePeakMZ"               "basePeakIntensity"       
#>  [5] "electronBeamEnergy"       "ionisationEnergy"        
#>  [7] "lowMZ"                    "highMZ"                  
#>  [9] "mergedScan"               "mergedResultScanNum"     
#> [11] "mergedResultStartScanNum" "mergedResultEndScanNum"  
#> [13] "injectionTime"            "filterString"            
#> [15] "spectrumId"               "ionMobilityDriftTime"    
#> [17] "scanWindowLowerLimit"     "scanWindowUpperLimit"    

## The values for these were filled with missing values for spectra from
## `s`:
all$peaksCount |> head()
#> [1]   NA   NA  578 1529 1600 1664


##  --------  AGGREGATE SPECTRA  --------

## Sets of spectra can be combined into a single, representative spectrum
## per set using `combineSpectra()`. This aggregates the peaks data (i.e.
## the spectra's m/z and intensity values) while using the values for all
## spectra variables from the first spectrum per set. Below we define the
## sets as all spectra measured in the *same second*, i.e. rounding their
## retention time to the next closer integer value.
f <- round(rtime(sciex))
head(f)
#> [1] 0 1 1 1 1 2

cmp <- combineSpectra(sciex, f = f)

## The length of `cmp` is now equal to the length of unique levels in `f`:
length(cmp)
#> [1] 29

## The spectra variable value from the first spectrum per set is used in
## the representative/combined spectrum:
cmp$rtime
#>  [1]  0.280  0.559  1.675  2.512  3.628  4.744  5.581  6.697  7.534  8.650
#> [11]  9.766 10.603 11.719 12.556 13.672 14.509 15.625 16.741 17.578 18.694
#> [21] 19.531 20.647 21.763 22.601 23.717 24.554 25.670 26.507 27.623

## The peaks data was aggregated: the number of mass peaks of the first six
## spectra from the original `Spectra`:
lengths(sciex) |> head()
#> [1]  578 1529 1600 1664 1417 1602

## and for the first aggreagated spectra:
lengths(cmp) |> head()
#>    0    1    2    3    4    5 
#>  578 3928 3177 3597 3928 3190 

## The default peaks data aggregation method joins all mass peaks. See
## documentation of the `combinePeaksData()` function for more options.


##  --------  SPLITTING DATA  --------

## A `Spectra` can be split into a `list` of `Spectra` objects using the
## `split()` function defining the sets into which the `Spectra` should
## be splitted into with parameter `f`.
sciex_split <- split(sciex, f)

length(sciex_split)
#> [1] 29
sciex_split |> head()
#> $`0`
#> MSn data (Spectra) with 1 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1      0.28         1
#>  ... 34 more variables/columns.
#> Processing:
#>  Switch backend from MsBackendMzR to MsBackendMemory [Tue Dec 23 13:20:26 2025] 
#> 
#> $`1`
#> MSn data (Spectra) with 4 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1     0.559         2
#> 2         1     0.838         3
#> 3         1     1.117         4
#> 4         1     1.396         5
#>  ... 34 more variables/columns.
#> Processing:
#>  Switch backend from MsBackendMzR to MsBackendMemory [Tue Dec 23 13:20:26 2025] 
#> 
#> $`2`
#> MSn data (Spectra) with 3 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1     1.675         6
#> 2         1     1.954         7
#> 3         1     2.233         8
#>  ... 34 more variables/columns.
#> Processing:
#>  Switch backend from MsBackendMzR to MsBackendMemory [Tue Dec 23 13:20:26 2025] 
#> 
#> $`3`
#> MSn data (Spectra) with 4 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1     2.512         9
#> 2         1     2.791        10
#> 3         1     3.070        11
#> 4         1     3.349        12
#>  ... 34 more variables/columns.
#> Processing:
#>  Switch backend from MsBackendMzR to MsBackendMemory [Tue Dec 23 13:20:26 2025] 
#> 
#> $`4`
#> MSn data (Spectra) with 4 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1     3.628        13
#> 2         1     3.907        14
#> 3         1     4.186        15
#> 4         1     4.465        16
#>  ... 34 more variables/columns.
#> Processing:
#>  Switch backend from MsBackendMzR to MsBackendMemory [Tue Dec 23 13:20:26 2025] 
#> 
#> $`5`
#> MSn data (Spectra) with 3 spectra in a MsBackendMemory backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1     4.744        17
#> 2         1     5.023        18
#> 3         1     5.302        19
#>  ... 34 more variables/columns.
#> Processing:
#>  Switch backend from MsBackendMzR to MsBackendMemory [Tue Dec 23 13:20:26 2025] 
#> 


##  --------  ADDING SPECTRA DATA  --------

## Adding new spectra variables
sciex1 <- filterDataOrigin(sciex, dataOrigin(sciex)[1])
spv <- DataFrame(spectrumId = sciex1$spectrumId[3:12], ## used for merging
                 var1 = rnorm(10),
                 var2 = sample(letters, 10))
spv
#> DataFrame with 10 rows and 3 columns
#>                spectrumId       var1        var2
#>               <character>  <numeric> <character>
#> 1  sample=1 period=1 cy..  1.1970777           p
#> 2  sample=1 period=1 cy.. -0.5486274           o
#> 3  sample=1 period=1 cy..  0.3030457           r
#> 4  sample=1 period=1 cy.. -0.0569705           q
#> 5  sample=1 period=1 cy.. -0.9578494           z
#> 6  sample=1 period=1 cy..  0.5910619           e
#> 7  sample=1 period=1 cy..  0.1731049           c
#> 8  sample=1 period=1 cy..  1.3997834           f
#> 9  sample=1 period=1 cy..  0.1174596           k
#> 10 sample=1 period=1 cy.. -0.3315458           x

sciex2 <- joinSpectraData(sciex1, spv, by.y = "spectrumId")

spectraVariables(sciex2)
#>  [1] "msLevel"                  "rtime"                   
#>  [3] "acquisitionNum"           "scanIndex"               
#>  [5] "dataStorage"              "dataOrigin"              
#>  [7] "centroided"               "smoothed"                
#>  [9] "polarity"                 "precScanNum"             
#> [11] "precursorMz"              "precursorIntensity"      
#> [13] "precursorCharge"          "collisionEnergy"         
#> [15] "isolationWindowLowerMz"   "isolationWindowTargetMz" 
#> [17] "isolationWindowUpperMz"   "peaksCount"              
#> [19] "totIonCurrent"            "basePeakMZ"              
#> [21] "basePeakIntensity"        "electronBeamEnergy"      
#> [23] "ionisationEnergy"         "lowMZ"                   
#> [25] "highMZ"                   "mergedScan"              
#> [27] "mergedResultScanNum"      "mergedResultStartScanNum"
#> [29] "mergedResultEndScanNum"   "injectionTime"           
#> [31] "filterString"             "spectrumId"              
#> [33] "ionMobilityDriftTime"     "scanWindowLowerLimit"    
#> [35] "scanWindowUpperLimit"     "var1"                    
#> [37] "var2"                    
spectraData(sciex2)[1:13, c("spectrumId", "var1", "var2")]
#> DataFrame with 13 rows and 3 columns
#>                 spectrumId      var1        var2
#>                <character> <numeric> <character>
#> 1   sample=1 period=1 cy..        NA          NA
#> 2   sample=1 period=1 cy..        NA          NA
#> 3   sample=1 period=1 cy..  1.197078           p
#> 4   sample=1 period=1 cy.. -0.548627           o
#> 5   sample=1 period=1 cy..  0.303046           r
#> ...                    ...       ...         ...
#> 9   sample=1 period=1 cy..  0.173105           c
#> 10  sample=1 period=1 cy..  1.399783           f
#> 11  sample=1 period=1 cy..  0.117460           k
#> 12  sample=1 period=1 cy.. -0.331546           x
#> 13  sample=1 period=1 cy..        NA          NA

## Append new spectra variables with cbind2()
df <- data.frame(cola = seq_len(length(sciex1)), colb = "b")
data_append <- cbind2(sciex1, df)
```
