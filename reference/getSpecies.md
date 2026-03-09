# Return Species in a Consortium

## Usage

``` r
getSpecies(object, ...)

# S4 method for class 'ConsortiumMetabolism'
getSpecies(object)

# S4 method for class 'ConsortiumMetabolismSet'
getSpecies(
  object,
  type = c("all", "generalists", "specialists"),
  quantileCutoff = 0.15
)

# S4 method for class 'ConsortiumMetabolismAlignment'
getSpecies(object)
```

## Arguments

- object:

  a `ConsortiumMetabolismAlignment` Object

- ...:

  Object specific arguments.

- type:

  Character scalar giving the type of species to output.

- quantileCutoff:

  Numeric scalar between 0 and 1 specifying the fraction of species to
  return when `type` is "generalists" or "specialists". For
  "generalists", the top `quantileCutoff` fraction of species with the
  most edges is returned. For "specialists", the bottom `quantileCutoff`
  fraction with the fewest edges is returned. Defaults to 0.15 (i.e.,
  15\\

A character vector containing the names of species in the consortiumA
character vector representing the microorganisms.A character vector
representing the microorganisms.A character vector representing the
microorganisms. Return Species in a Consortium Methods (by class)

- `getSpecies(ConsortiumMetabolism)`: Return Species in a Microbiome

- `getSpecies(ConsortiumMetabolismSet)`: Return Species in a Microbiome

- `getSpecies(ConsortiumMetabolismAlignment)`: Return Species in a
  Microbiome
