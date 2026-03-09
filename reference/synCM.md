# Generate Random Synthetic Consortium Metabolism

A function that creates synthetic data suitable for demonstration
purposes of the `ramen` package.

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

  Character string giving the desired name of the community.

- n_species:

  Number of species in the community

- max_met:

  Maximum number of metabolites in the communities

- scale_fac:

  Scaling factor

- seed:

  Seed for reproducibility

- dead_ends:

  Logical value to toggle dead ends in data

- cm:

  Logical value to toggle return of `ConsortiumMetabolism` object or
  tibble.

## Value

List with `n_co` number of communities.

## Examples

``` r
synCM("Ex. Community", n_species = 5, max_met = 10)
#> Ex. Community: ConsortiumMetabolism Object
#> Weighted metabolic network with 10 metabolites.
```
