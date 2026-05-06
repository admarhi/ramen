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

For `ConsortiumMetabolismAlignment` objects, `type` selects the pathway
subset:

- `"all"` returns the union of all pathways

- `"shared"` (pairwise only) returns pathways shared between query and
  reference

- `"unique"` (pairwise only) returns pathways unique to query and
  reference as a list

- `"consensus"` (multiple only) returns consensus network pathways with
  prevalence

## Usage

``` r
pathways(object, ...)

# S4 method for class 'ConsortiumMetabolism'
pathways(object, verbose = FALSE)

# S4 method for class 'ConsortiumMetabolismAlignment'
pathways(
  object,
  type = c("all", "shared", "unique", "consensus"),
  verbose = FALSE
)

# S4 method for class 'ConsortiumMetabolismSet'
pathways(
  object,
  type = c("all", "pan-cons", "niche", "core", "aux"),
  quantileCutoff = 0.1,
  verbose = FALSE
)
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

  Character scalar giving the type of pathways to output. Note that the
  four filtered types use two different rulers: `"pan-cons"` and
  `"niche"` rank by the integer per-pathway consortium count and
  quantile over a uniform integer ruler (`seq_len(total_cons)`);
  `"core"` and `"aux"` rank by per-pathway species count and quantile
  over the empirical `n_species` distribution across pathways.

- quantileCutoff:

  Numeric scalar between 0 and 1 giving the quantile threshold to use
  for filtering pathways. For `"pan-cons"` and `"core"` types, pathways
  strictly above `1 - quantileCutoff` are returned. For `"niche"` and
  `"aux"` types, pathways at or below `quantileCutoff` are returned (the
  boundary is included on the lower tail so that ties at the quantile
  floor are not silently dropped). Defaults to 0.1 (i.e., top/bottom 10
  percent).

## Value

A data.frame of pathway information. With `verbose = FALSE` (default):
`consumed`, `produced`, `n_species` (and `n_cons` for CMS objects). With
`verbose = TRUE`: all available columns including flux statistics and
indices. For CMA objects, the return depends on `type`: a data.frame for
`"all"`, `"shared"`, and `"consensus"`; a list of data.frames for
`"unique"`.

## Methods (by class)

- `pathways(ConsortiumMetabolism)`: Get Pathways From a
  `ConsortiumMetabolism` Object

- `pathways(ConsortiumMetabolismAlignment)`: Get Pathways From a
  `ConsortiumMetabolismAlignment` Object

- `pathways(ConsortiumMetabolismSet)`: Get Pathways From a
  `ConsortiumMetabolismSet` Object

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
pathways(cm)
#> # A tibble: 10 × 3
#>    consumed produced n_species
#>    <chr>    <chr>        <dbl>
#>  1 met4     media            1
#>  2 met4     met5             1
#>  3 met4     met1             1
#>  4 met5     met4             1
#>  5 met5     met1             1
#>  6 met3     met5             1
#>  7 met3     met1             2
#>  8 met3     met4             1
#>  9 met2     met4             1
#> 10 met2     met1             1
pathways(cm, verbose = TRUE)
#> # A tibble: 10 × 14
#>    consumed produced data     n_species  c_sum p_sum c_prob p_prob c_neff p_neff
#>    <chr>    <chr>    <list>       <dbl>  <dbl> <dbl> <list> <list>  <dbl>  <dbl>
#>  1 met4     media    <tibble>         1  0.384 1     <dbl>  <dbl>    1      1   
#>  2 met4     met5     <tibble>         1  1.91  1.34  <dbl>  <dbl>    1      1   
#>  3 met4     met1     <tibble>         1  1.91  3.94  <dbl>  <dbl>    1      1   
#>  4 met5     met4     <tibble>         1  1.34  2.29  <dbl>  <dbl>    1      1   
#>  5 met5     met1     <tibble>         1  1.34  0.790 <dbl>  <dbl>    1      1   
#>  6 met3     met5     <tibble>         1  7.00  1.34  <dbl>  <dbl>    1      1   
#>  7 met3     met1     <tibble>         2 10.2   4.73  <dbl>  <dbl>    1.86   1.57
#>  8 met3     met4     <tibble>         1  3.18  2.29  <dbl>  <dbl>    1      1   
#>  9 met2     met4     <tibble>         1  0.826 2.29  <dbl>  <dbl>    1      1   
#> 10 met2     met1     <tibble>         1  0.826 0.790 <dbl>  <dbl>    1      1   
#> # ℹ 4 more variables: c_eff <dbl>, p_eff <dbl>, c_ind <int>, p_ind <int>

cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cma <- align(cm1, cm2)
pathways(cma)
#> # A tibble: 0 × 2
#> # ℹ 2 variables: consumed <chr>, produced <chr>
pathways(cma, type = "shared")
#> # A tibble: 0 × 2
#> # ℹ 2 variables: consumed <chr>, produced <chr>
pathways(cma, type = "unique")
#> $query
#> # A tibble: 5 × 2
#>   consumed produced
#>   <chr>    <chr>   
#> 1 met1     media   
#> 2 met2     met1    
#> 3 met1     met2    
#> 4 met1     met3    
#> 5 met2     met5    
#> 
#> $reference
#> # A tibble: 3 × 2
#>   consumed produced
#>   <chr>    <chr>   
#> 1 met4     met1    
#> 2 met1     met4    
#> 3 met5     met4    
#> 
```
