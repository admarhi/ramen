# Pivot `ConsortiaMetabolism` Input Data

Wrapper function around tidyr's `pivot_longer()` function to facilitate
the easy transformation into the correct data format for MiCo
construction.

## Usage

``` r
pivotCM(tb, species, from, to, flux)
```

## Arguments

- tb:

  Tibble with data on a microbial community in long or short format to
  be

- species:

  Column name of the species column

- from:

  the name of the column specifying the met consumed

- to:

  the name of the column specifying the met excreted

- flux:

  Column name of the flux column

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
