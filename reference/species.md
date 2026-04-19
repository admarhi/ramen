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
species(object, ...)
```

## Arguments

- object:

  A `ConsortiumMetabolismSet` object.

- ...:

  Additional arguments passed to methods.

## Value

A character vector of species names.

## Methods (by class)

- `species(ConsortiumMetabolism)`: Return Species in a
  `ConsortiumMetabolism`

- `species(ConsortiumMetabolismAlignment)`: Return species from a
  `ConsortiumMetabolismAlignment`

- `species(ConsortiumMetabolismSet)`: Return Species in a
  `ConsortiumMetabolismSet`

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
species(cm)
#> [1] "AJU3814B" "COK7774C" "YPS2470L"
```
