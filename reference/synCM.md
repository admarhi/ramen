# Generate Synthetic Consortium Metabolism

Creates a synthetic metabolic network with random species-metabolite
interactions and flux values. Useful for testing and demonstrating
`ramen` package functionality.

## Usage

``` r
synCM(
  name,
  n_species,
  max_met,
  scale_fac = 2,
  seed = FALSE,
  dead_ends = FALSE,
  cm = TRUE
)
```

## Arguments

- name:

  Character string. Name for the consortium.

- n_species:

  Integer. Number of species in the consortium.

- max_met:

  Integer. Size of the metabolite pool to sample from.

- scale_fac:

  Integer. Multiplier for the initial species name pool from which
  `n_species` are sampled. Defaults to 2.

- seed:

  Integer or `FALSE`. Random seed for reproducibility. If `FALSE`
  (default), a random seed is chosen.

- dead_ends:

  Logical. If `FALSE` (default), ensures each species has both consumed
  and produced metabolites by flipping one flux value when all fluxes
  share the same sign.

- cm:

  Logical. If `TRUE` (default), returns a
  [`ConsortiumMetabolism`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
  object. If `FALSE`, returns the raw edge list as a tibble.

## Value

A
[`ConsortiumMetabolism`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
object when `cm = TRUE`, or a
[`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
columns `species`, `metabolites`, and `fluxes` when `cm = FALSE`.

## Details

Each species is assigned a random subset of metabolites (between 2 and
`max_met`) with normally distributed flux values. By default, species
are guaranteed to have both positive and negative fluxes (i.e., no dead
ends).

## Examples

``` r
synCM("Ex. Community", n_species = 5, max_met = 10)
#> 
#> ── ConsortiumMetabolism 
#> Name: "Ex. Community"
#> Weighted metabolic network with 10 metabolites.

# Return raw edge list instead
synCM("Test", n_species = 3, max_met = 8, cm = FALSE)
#> # A tibble: 15 × 3
#>    species  metabolites fluxes
#>    <chr>    <chr>        <dbl>
#>  1 SCE6853Y met3         1.03 
#>  2 SCE6853Y met1        -1.91 
#>  3 PVO1450E met6         5.15 
#>  4 PVO1450E met3         2.32 
#>  5 PVO1450E met7         3.67 
#>  6 PVO1450E met8         0.309
#>  7 PVO1450E met4         2.33 
#>  8 PVO1450E met1         0.338
#>  9 PVO1450E met5        -4.77 
#> 10 ZRK8603K met4        -1.38 
#> 11 ZRK8603K met7         5.40 
#> 12 ZRK8603K met8         4.35 
#> 13 ZRK8603K met6        -0.899
#> 14 ZRK8603K met5        -0.685
#> 15 ZRK8603K met2        -0.965
```
