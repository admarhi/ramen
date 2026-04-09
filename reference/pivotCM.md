# Pivot `ConsortiumMetabolism` Input Data

Wrapper function around tidyr's `pivot_longer()` function to transform
wide-format community data (one column per direction) into the long
format required by
[`ConsortiumMetabolism`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md).

## Usage

``` r
pivotCM(tb, species, from, to, flux)
```

## Arguments

- tb:

  A data.frame or tibble with one row per species-metabolite pair in
  wide format (separate columns for consumed and produced metabolites).

- species:

  Column name of the species column.

- from:

  Column name specifying the metabolite consumed.

- to:

  Column name specifying the metabolite produced.

- flux:

  Column name of the flux column.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
three columns:

- species:

  Character, the species identifier.

- met:

  Character, the metabolite name.

- flux:

  Numeric, the metabolic flux (negative for consumption, positive for
  production).

## Examples

``` r
tb <- tibble::tribble(
  ~uptake, ~secretion, ~flux, ~species,
  "m1", "m2", 1, "s1",
  "m2", "m3", 1, "s2",
  "m3", "m4", 1, "s3",
  "m4", "m1", 1, "s4"
)

pivotCM(
  tb = tb,
  species = "species",
  from = "uptake",
  to = "secretion",
  flux = "flux"
)
#> # A tibble: 8 × 3
#>   species met    flux
#>   <chr>   <chr> <dbl>
#> 1 s1      m1       -1
#> 2 s1      m2        1
#> 3 s2      m2       -1
#> 4 s2      m3        1
#> 5 s3      m3       -1
#> 6 s3      m4        1
#> 7 s4      m4       -1
#> 8 s4      m1        1
```
