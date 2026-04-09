# Generate Synthetic Consortium Metabolism

Creates a synthetic metabolic network with biologically realistic
structure: hub metabolites (Zipf-distributed degree), log-normal species
degrees, cyclic cross-feeding backbone for connectivity, approximate
mass balance, and no dead-end species by default.

## Usage

``` r
synCM(name, n_species, max_met, seed = FALSE, dead_ends = FALSE, cm = TRUE)
```

## Arguments

- name:

  Character string. Name for the consortium.

- n_species:

  Integer. Number of species in the consortium.

- max_met:

  Integer. Size of the metabolite pool to sample from.

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

## Examples

``` r
synCM("Ex. Community", n_species = 5, max_met = 10)
#> 
#> ── ConsortiumMetabolism 
#> Name: "Ex. Community"
#> Weighted metabolic network: 5 species, 10 metabolites, 28 pathways.

# Return raw edge list instead
synCM("Test", n_species = 3, max_met = 8, cm = FALSE)
#> # A tibble: 11 × 3
#>    species  metabolites fluxes
#>    <chr>    <chr>        <dbl>
#>  1 WQN4027C met7         1.20 
#>  2 ADI9989E met7        -1.20 
#>  3 ADI9989E met2         2.29 
#>  4 LMR532J  met2        -2.29 
#>  5 LMR532J  met7         2.91 
#>  6 WQN4027C met7        -2.91 
#>  7 ADI9989E met5        -7.50 
#>  8 ADI9989E met4        -1.88 
#>  9 ADI9989E met3         0.679
#> 10 LMR532J  met5        -3.28 
#> 11 LMR532J  met6         1.66 
```
