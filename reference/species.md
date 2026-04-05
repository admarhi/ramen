# Return Species in a Consortium

Returns the species present in a consortium or set of consortia.

## Usage

``` r
species(object, ...)

# S4 method for class 'ConsortiumMetabolism'
species(object)

# S4 method for class 'ConsortiumMetabolismSet'
species(
  object,
  type = c("all", "generalists", "specialists"),
  quantileCutoff = 0.15
)

# S4 method for class 'ConsortiumMetabolismAlignment'
species(object)
```

## Arguments

- object:

  A `ConsortiumMetabolismAlignment` object.

- ...:

  Additional arguments passed to methods.

- type:

  Character scalar giving the type of species to output.

- quantileCutoff:

  Numeric scalar between 0 and 1 specifying the fraction of species to
  return when `type` is "generalists" or "specialists". For
  "generalists", the top `quantileCutoff` fraction of species with the
  most pathways is returned. For "specialists", the bottom
  `quantileCutoff` fraction with the fewest pathways is returned.
  Defaults to 0.15 (i.e., 15 percent). Ignored when `type = "all"`.

## Value

A character vector of species names, or a tibble with species and
pathway counts (for
[ConsortiumMetabolismSet](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md)).

A character vector representing the microorganisms.

A character vector representing the microorganisms.

A character vector of species names.

## Methods (by class)

- `species(ConsortiumMetabolism)`: Return Species in a Microbiome

- `species(ConsortiumMetabolismSet)`: Return Species in a Microbiome

- `species(ConsortiumMetabolismAlignment)`: Return species from a
  `ConsortiumMetabolismAlignment`

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
species(cm)
#> [1] "RZQ1660H" "TCB696S"  "IRL7860P"
```
