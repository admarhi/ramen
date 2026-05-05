# Get Metabolites

Retrieves the metabolites involved in the metabolic network. For
`ConsortiumMetabolism` objects, the result can optionally be restricted
to a specific species and/or direction (`"consumed"` or `"produced"`).

## Usage

``` r
metabolites(object, ...)

# S4 method for class 'ConsortiumMetabolism'
metabolites(
  object,
  species = NULL,
  direction = c("all", "consumed", "produced")
)

# S4 method for class 'ConsortiumMetabolismAlignment'
metabolites(object)

# S4 method for class 'ConsortiumMetabolismSet'
metabolites(object)
```

## Arguments

- object:

  A `ConsortiumMetabolism`, `ConsortiumMetabolismSet`, or
  `ConsortiumMetabolismAlignment` object.

- ...:

  Additional arguments. For `ConsortiumMetabolism`: `species` (character
  scalar; restrict to metabolites involved with this species) and
  `direction` (one of `"all"`, `"consumed"`, or `"produced"`; defaults
  to `"all"`).

- species:

  Optional length-1 character scalar. If supplied, restrict the result
  to metabolites involved in pathways that include this species.
  Defaults to `NULL` (all species).

- direction:

  One of `"all"` (default), `"consumed"`, or `"produced"`. Restricts the
  result to metabolites in the given role across the (possibly
  species-filtered) pathways.

## Value

A character vector containing the names of metabolites in the network.

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
metabolites(cm)
#> [1] "met5" "met4" "met3" "met2" "met1"
## Metabolites consumed by a specific species:
sp <- species(cm)[1]
metabolites(cm, species = sp, direction = "consumed")
#> [1] "met3" "met5"
```
