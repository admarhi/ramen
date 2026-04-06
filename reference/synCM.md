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
#> Weighted metabolic network with 7 metabolites.

# Return raw edge list instead
synCM("Test", n_species = 3, max_met = 8, cm = FALSE)
#> # A tibble: 10 × 3
#>    species  metabolites fluxes
#>    <chr>    <chr>        <dbl>
#>  1 KYS3742B met5         3.77 
#>  2 PCI233C  met5        -3.77 
#>  3 PCI233C  met4         0.875
#>  4 MKI7263L met4        -0.875
#>  5 MKI7263L met5         1.48 
#>  6 KYS3742B met5        -1.48 
#>  7 KYS3742B met6        -0.894
#>  8 MKI7263L met7         0.904
#>  9 MKI7263L met8         2.87 
#> 10 MKI7263L met6        -4.32 
```
