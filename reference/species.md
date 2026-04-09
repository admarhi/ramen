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
  quantileCutoff = 0.15
)
```

## Arguments

- object:

  a `ConsortiumMetabolismSet` Object

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

For
[ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
and
[ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md),
a character vector of species names. For
[ConsortiumMetabolismSet](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md),
a tibble with columns `species` and `n_pathways`.

## Methods (by class)

- `species(ConsortiumMetabolism)`: Return Species in a Microbiome

- `species(ConsortiumMetabolismAlignment)`: Return species from a
  `ConsortiumMetabolismAlignment`

- `species(ConsortiumMetabolismSet)`: Return Species in a Microbiome

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
species(cm)
#> [1] "XEW6572A" "SHJ4890C" "MQL2872Q"
```
