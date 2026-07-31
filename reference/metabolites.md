# Get Metabolites

Retrieves the metabolites involved in the metabolic network. For
`ConsortiumMetabolism` objects, the result can optionally be restricted
to a specific species and/or direction (`"consumed"` or `"produced"`).

For a `ConsortiumMetabolism`, the default (`species = NULL`,
`direction = "all"`) returns the full metabolite axis of the assays,
i.e. it is identical to `rownames(assay(object))`. This includes the
`"media"` boundary node whenever the consortium retains it (a consortium
in which every species both consumes and produces has no `"media"`
node). Indexing an assay with the result, e.g.
`assay(object)[metabolites(object), metabolites(object)]`, is therefore
always well-formed.

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
For a `ConsortiumMetabolism` with the default arguments this equals
`rownames(assay(object))` and may include the `"media"` boundary node.

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
metabolites(cm)
#> [1] "met1"  "met2"  "met5"  "media"
## Metabolites consumed by a specific species:
sp <- species(cm)[1]
metabolites(cm, species = sp, direction = "consumed")
#> [1] "media"
```
