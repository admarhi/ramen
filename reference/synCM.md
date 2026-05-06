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
#> Weighted metabolic network: 5 species, 6 metabolites, 10 pathways.
#> Pathways per species: min 1, mean 2.2, max 6.

# Return raw edge list instead
synCM("Test", n_species = 3, max_met = 8, cm = FALSE)
#> # A tibble: 8 × 3
#>   species  metabolites fluxes
#>   <chr>    <chr>        <dbl>
#> 1 OQZ8515G met4        -0.713
#> 2 APF6467H met4        -0.888
#> 3 APF6467H met2         0.886
#> 4 UOH9021W met2        -0.886
#> 5 UOH9021W met4         1.60 
#> 6 UOH9021W met8        -2.57 
#> 7 UOH9021W met3         4.64 
#> 8 UOH9021W met6        -3.19 
```
