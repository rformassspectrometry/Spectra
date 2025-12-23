# Apply a function stepwise to chunks of data

`chunkapply()` splits `x` into chunks and applies the function `FUN`
stepwise to each of these chunks. Depending on the object it is called,
this function might reduce memory demand considerably, if for example
only the full data for a single chunk needs to be loaded into memory at
a time (e.g., for `Spectra` objects with on-disk or similar backends).

## Usage

``` r
chunkapply(x, FUN, ..., chunkSize = 1000L, chunks = factor())
```

## Arguments

- x:

  object to which `FUN` should be applied. Can be any object that
  supports `split`.

- FUN:

  the function to apply to `x`.

- ...:

  additional parameters to `FUN`.

- chunkSize:

  `integer(1)` defining the size of each chunk into which `x` should be
  splitted.

- chunks:

  optional `factor` or length equal to `length(x)` defining the chunks
  into which `x` should be splitted.

## Value

Depending on `FUN`, but in most cases a vector/result object of length
equal to `length(x)`.

## Author

Johannes Rainer

## Examples

``` r

## Apply a function (`sqrt`) to each element in `x`, processed in chunks of
## size 200.
x <- rnorm(n = 1000, mean = 500)
res <- chunkapply(x, sqrt, chunkSize = 200)
length(res)
#> [1] 1000
head(res)
#> [1] 22.37457 22.38634 22.31991 22.35515 22.35522 22.35436

## For such a calculation the vectorized `sqrt` would however be recommended
system.time(sqrt(x))
#>    user  system elapsed 
#>   0.001   0.000   0.000 
system.time(chunkapply(x, sqrt, chunkSize = 200))
#>    user  system elapsed 
#>   0.001   0.000   0.001 

## Simple example splitting a numeric vector into chunks of 200 and
## aggregating the values within the chunk using the `mean`. Due to the
## `unsplit` the result has the same length than the input with the mean
## value repeated.
x <- 1:1000
res <- chunkapply(x, mean, chunkSize = 200)
length(res)
#> [1] 1000
head(res)
#> [1] 100.5 100.5 100.5 100.5 100.5 100.5
```
