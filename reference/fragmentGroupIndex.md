# Mass fragmentation collections of each full scan

This function generates an `integer` index grouping MS^n spectra (MS
level \> 1) with their corresponding MS1 spectra based on acquisition
order. Each group contains exactly one MS1 spectrum and all subsequent
higher-level spectra (MS2, MS3, ...) acquired until the next MS1 scan.
MS1-only spectra are also assigned sequential group IDs.

Note that this function:

- does not consider the direct relationship between a precursor scan and
  the associated product scans,

- and does not distinguish between different fragmentation trees.

For example, all MS3 scans measured after a given MS1 are grouped
together with all MS2 scans from that MS1, regardless of which MS2
spectrum they originated from. See
[`filterPrecursorScan()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
for a function that considers relationships between fragment and
precursor scans.

## Usage

``` r
fragmentGroupIndex(object, BPPARAM = SerialParam())
```

## Arguments

- object:

  A `Spectra` object (from the **Spectra** package) containing MS data.
  Must include at least two MS levels (`msLevel`) and be ordered by
  `acquisitionNum` within each `dataOrigin`.

- BPPARAM:

  A `BiocParallelParam` object for parallel execution. Defaults to
  [`SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html).

## Value

An `integer` vector of the same length as `object`. Each element gives
the group index associated with the corresponding spectrum (MS1 or
MS^n). Group indices are unique across all files (`dataOrigin` values).

## Note

- Each file (`dataOrigin`) must contain at least one MS1 spectrum.

- If a group contains only MS1 spectra, each MS1 is assigned a unique
  group ID.

- The user is responsible for ensuring that spectra are correctly
  ordered. Improper ordering may lead to incorrect groupings.

## See also

[`filterPrecursorScan()`](https://rformassspectrometry.github.io/Spectra/reference/filterMsLevel.md)
for a function that instead returns a `Spectra` object containing each
parent (e.g., MS1) and its direct child scans (e.g., MS2) according to
their acquisition numbers.

## Author

Philippine Louail

## Examples

``` r

fl_ms3 <- system.file("proteomics", "MS3TMT11.mzML", package = "msdata")
sps_dda <- Spectra(fl_ms3)
idx <- fragmentGroupIndex(sps_dda)
head(idx)
#>             
#> 1 1 1 1 1 1 
```
