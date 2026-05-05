# Return Species in a Consortium

Returns the species present in a consortium or set of consortia.

## Usage

``` r
species(object, ...)

# S4 method for class 'ConsortiumMetabolism'
species(object)

# S4 method for class 'ConsortiumMetabolismAlignment'
species(object)

# S4 method for class 'ConsortiumMetabolismSet'
species(
  object,
  type = c("all", "generalists", "specialists"),
  quantileCutoff = 0.15,
  ...
)
```

## Arguments

- object:

  A `ConsortiumMetabolismSet` object.

- ...:

  Additional arguments passed to methods.

- type:

  Character scalar. One of `"all"` (default), `"generalists"` (top
  fraction by pathway count), or `"specialists"` (bottom fraction).

- quantileCutoff:

  Numeric scalar in (0, 1) giving the fraction of species to label as
  generalists or specialists. Defaults to `0.15`.

## Value

A character vector of species names.

## Methods (by class)

- `species(ConsortiumMetabolism)`: Return Species in a
  `ConsortiumMetabolism`

- `species(ConsortiumMetabolismAlignment)`: Return species from a
  `ConsortiumMetabolismAlignment`

- `species(ConsortiumMetabolismSet)`: Return Species in a
  `ConsortiumMetabolismSet`, optionally filtered to metabolic
  generalists or specialists by pathway count.

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
species(cm)
#> [1] "ZEQ5398Z" "UOG489Y"  "MXO4181U"
```
