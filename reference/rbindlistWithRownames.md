# Fast *rbind-ing* `data.frame`s preserving row names

The `rbindlistWithRownames()` function uses the
[`data.table::rbindlist()`](https://rdrr.io/pkg/data.table/man/rbindlist.html)
function for a fast concatenation (row-binding) of a `list` of provided
`data.frame`s and in addition also preserves their row names (if set).

The parameters are the same as for `rbindlist()`.

## Usage

``` r
rbindlistWithRownames(
  l,
  use.names = TRUE,
  fill = FALSE,
  idcol = NULL,
  ignore.attr = FALSE
)
```

## Arguments

- l:

  A `list` containing `data.frame`s that should be combined.

- use.names:

  If `TRUE` binds by matching column names, if `FALSE` by position.
  `"check"` (default) warns if not all items have the same names in the
  same order. See
  [`data.table::rbindlist()`](https://rdrr.io/pkg/data.table/man/rbindlist.html)
  for details.

- fill:

  `logical(1)`: if `TRUE` fills missing columns with NAs or `NULL` for
  missing list columns. By defalt `FALSE`.

- idcol:

  Creates a column in the result showing which list item those rows cam
  from. See
  [`data.table::rbindlist()`](https://rdrr.io/pkg/data.table/man/rbindlist.html)
  for details.

- ignore.attr:

  `logical(1)`, default `FALSE`. When `TRUE`, allows binding columns
  with different attributes (e.g. class).

## Value

A combined `data.frame` with row names (if present in all `data.frame`s
provided).

## Note

Row names are dropped if duplicated row names are present, or if row
names are not defined for all `data.frame`s in `l`.

The function uses
[`.row_names_info()`](https://rdrr.io/r/base/base-internal.html) to
guess wheter row names are *really* set or just the default
`seq_len(nrow(x))` are used.

## See also

[`data.table::rbindlist()`](https://rdrr.io/pkg/data.table/man/rbindlist.html)

## Author

Johannes Rainer
