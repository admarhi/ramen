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
#> # A tibble: 7 × 3
#>   consumed produced n_species
#>   <chr>    <chr>        <dbl>
#> 1 met5     met3             1
#> 2 met5     met1             1
#> 3 met3     met5             1
#> 4 met2     met5             1
#> 5 met2     met3             1
#> 6 met4     met3             1
#> 7 met4     met1             1
pathways(cm, verbose = TRUE)
#> # A tibble: 7 × 12
#>   consumed produced data     n_species c_sum p_sum c_prob    p_prob c_eff p_eff
#>   <chr>    <chr>    <list>       <dbl> <dbl> <dbl> <list>    <list> <dbl> <dbl>
#> 1 met5     met3     <tibble>         1  2.76  2.89 <dbl [1]> <dbl>      1     1
#> 2 met5     met1     <tibble>         1  2.76  2.45 <dbl [1]> <dbl>      1     1
#> 3 met3     met5     <tibble>         1  6.36  1.14 <dbl [1]> <dbl>      1     1
#> 4 met2     met5     <tibble>         1  1.20  1.62 <dbl [1]> <dbl>      1     1
#> 5 met2     met3     <tibble>         1  1.20  3.45 <dbl [1]> <dbl>      1     1
#> 6 met4     met3     <tibble>         1  3.08  2.89 <dbl [1]> <dbl>      1     1
#> 7 met4     met1     <tibble>         1  3.08  2.45 <dbl [1]> <dbl>      1     1
#> # ℹ 2 more variables: c_ind <int>, p_ind <int>

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
#> # A tibble: 11 × 2
#>    consumed produced
#>    <chr>    <chr>   
#>  1 met2     met1    
#>  2 met3     met1    
#>  3 met4     met1    
#>  4 met5     met1    
#>  5 met2     met3    
#>  6 met4     met3    
#>  7 met2     met4    
#>  8 met3     met4    
#>  9 met5     met4    
#> 10 met2     met5    
#> 11 met4     met5    
#> 
#> $reference
#> # A tibble: 5 × 2
#>   consumed produced
#>   <chr>    <chr>   
#> 1 met1     met4    
#> 2 met6     met4    
#> 3 met1     met5    
#> 4 met6     met5    
#> 5 met4     met6    
#> 
```
