# Species Summary

Returns an enriched per-species summary as a tibble. Provides more
detail than
[`species()`](https://admarhi.github.io/ramen/reference/species.md),
which returns only a character vector.

## Usage

``` r
speciesSummary(object, ...)

# S4 method for class 'ConsortiumMetabolism'
speciesSummary(object, ...)

# S4 method for class 'ConsortiumMetabolismAlignment'
speciesSummary(object, ...)

# S4 method for class 'ConsortiumMetabolismSet'
speciesSummary(object, ...)
```

## Arguments

- object:

  A `ConsortiumMetabolismSet` object.

- ...:

  Additional arguments passed to methods.

## Value

A tibble with per-species metrics. Columns depend on the class:

- `ConsortiumMetabolism`: `species`, `n_pathways`, `n_consumed`,
  `n_produced`.

- `ConsortiumMetabolismSet`: `species`, `n_consortia`, `n_pathways`.

- `ConsortiumMetabolismAlignment` (pairwise): `species`, `role`
  (`"shared"`, `"unique_query"`, or `"unique_reference"`).

## Methods (by class)

- `speciesSummary(ConsortiumMetabolism)`: Species summary for a
  `ConsortiumMetabolism`

- `speciesSummary(ConsortiumMetabolismAlignment)`: Species summary for a
  `ConsortiumMetabolismAlignment`. Only available for pairwise
  alignments.

- `speciesSummary(ConsortiumMetabolismSet)`: Species summary for a
  `ConsortiumMetabolismSet`

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
speciesSummary(cm)
#> # A tibble: 3 × 4
#>   species  n_pathways n_consumed n_produced
#>   <chr>         <int>      <int>      <int>
#> 1 FXY9545K          6          2          3
#> 2 RIJ505L           6          3          2
#> 3 MRP6146F          1          1          1
```
