# Retrieve Metabolic Pathways

Retrieves the pathways representing metabolic interactions between
species.

By default, returns a concise summary with columns `consumed`,
`produced`, and `n_species`. For `ConsortiumMetabolismSet` objects,
`n_cons` is included as well. Set `verbose = TRUE` to return the full
pathway data including flux statistics, indices, and per-species detail.

The argument `type` can be used to return only specific types of
pathways from a `ConsortiumMetabolismSet` object:

- `"all"` returns all pathways

- `"pan-cons"` returns pathways present in most consortia

- `"niche"` returns niche pathways specific to few consortia

- `"core"` returns core metabolic pathways shared across most species

- `"aux"` returns auxiliary pathways found in few species

## Usage

``` r
pathways(object, ...)

# S4 method for class 'ConsortiumMetabolism'
pathways(object, verbose = FALSE)

# S4 method for class 'ConsortiumMetabolismSet'
pathways(
  object,
  type = c("all", "pan-cons", "niche", "core", "aux"),
  quantileCutoff = 0.1,
  verbose = FALSE
)

# S4 method for class 'ConsortiumMetabolismAlignment'
pathways(object)
```

## Arguments

- object:

  A `ConsortiumMetabolism`, `ConsortiumMetabolismSet`, or
  `ConsortiumMetabolismAlignment` object.

- ...:

  Object specific arguments. See methods for details.

- verbose:

  Logical scalar. If `FALSE` (default), returns a concise summary with
  columns `consumed`, `produced`, and `n_species`. If `TRUE`, returns
  the full pathway data including flux statistics, effective values, and
  per-species detail.

- type:

  Character scalar giving the type of pathways to output.

- quantileCutoff:

  Numeric scalar between 0 and 1 giving the quantile threshold to use
  for filtering pathways. For `"pan-cons"` and `"core"` types, pathways
  above `1 - quantileCutoff` are returned. For `"niche"` and `"aux"`
  types, pathways below `quantileCutoff` are returned. Defaults to 0.1
  (i.e., top/bottom 10 percent).

## Value

A data.frame of pathway information. With `verbose = FALSE` (default):
`consumed`, `produced`, `n_species` (and `n_cons` for CMS objects). With
`verbose = TRUE`: all available columns including flux statistics and
indices.

## Methods (by class)

- `pathways(ConsortiumMetabolism)`: Get Pathways From a
  `ConsortiumMetabolism` Object

- `pathways(ConsortiumMetabolismSet)`: Get Pathways From a
  `ConsortiumMetabolismSet` Object

- `pathways(ConsortiumMetabolismAlignment)`: Get Pathways From a
  `ConsortiumMetabolismAlignment` Object

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
pathways(cm)
#> # A tibble: 5 × 3
#>   consumed produced n_species
#>   <chr>    <chr>        <dbl>
#> 1 met5     met2             2
#> 2 met4     met3             1
#> 3 met5     met3             1
#> 4 met2     met3             1
#> 5 met4     met2             1
pathways(cm, verbose = TRUE)
#> # A tibble: 5 × 12
#>   consumed produced data     n_species  c_sum p_sum c_prob    p_prob c_eff p_eff
#>   <chr>    <chr>    <list>       <dbl>  <dbl> <dbl> <list>    <list> <dbl> <dbl>
#> 1 met5     met2     <tibble>         2 11.8    5.01 <dbl [2]> <dbl>   1.92  1.86
#> 2 met4     met3     <tibble>         1  0.762  1.76 <dbl [1]> <dbl>   1     1   
#> 3 met5     met3     <tibble>         1  1.09   1.76 <dbl [1]> <dbl>   1     1   
#> 4 met2     met3     <tibble>         1  5.25   1.76 <dbl [1]> <dbl>   1     1   
#> 5 met4     met2     <tibble>         1  4.65   3.45 <dbl [1]> <dbl>   1     1   
#> # ℹ 2 more variables: c_ind <int>, p_ind <int>
```
