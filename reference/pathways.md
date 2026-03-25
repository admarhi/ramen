# Retrieve Metabolic Pathways

Retrieves the pathways representing metabolic interactions between
species. The argument `type` can be used to return only specific types
of pathways from a `ConsortiumMetabolismSet` type object.

- `all` will return all pathways

- `pan-cons` will return only pathways that exist in all consortia

- `niche` will return niche pathways specific to individual consortia

- `core` will return core metabolic pathways

- `aux` will return auxiliary pathways

## Usage

``` r
pathways(object, ...)

# S4 method for class 'ConsortiumMetabolism'
pathways(object)

# S4 method for class 'ConsortiumMetabolismSet'
pathways(
  object,
  type = c("all", "pan-cons", "niche", "core", "aux"),
  quantileCutoff = 0.1
)

# S4 method for class 'ConsortiumMetabolismAlignment'
pathways(object)
```

## Arguments

- object:

  A `ConsortiumMetabolism` object.

- ...:

  Object specific arguments.

- type:

  Character scalar giving the type of pathways to output.

- quantileCutoff:

  Numeric scalar between 0 and 1 giving the quantile threshold to use
  for filtering pathways. For "pan-cons" and "core" types, pathways
  above `1 - quantileCutoff` are returned. For "niche" and "aux" types,
  pathways below `quantileCutoff` are returned. Defaults to 0.1 (i.e.,
  top/bottom 10 percent).

## Value

A tibble containing pathway information including:

- consumed/produced metabolites

- number of species involved

- consumption/production sums

- effective consumption/production metrics

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
