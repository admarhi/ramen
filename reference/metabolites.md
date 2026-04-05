# Get Metabolites

Retrieves the metabolites involved in the metabolic network.

## Usage

``` r
metabolites(object)

# S4 method for class 'ConsortiumMetabolism'
metabolites(object)

# S4 method for class 'ConsortiumMetabolismSet'
metabolites(object)

# S4 method for class 'ConsortiumMetabolismAlignment'
metabolites(object)
```

## Arguments

- object:

  A `ConsortiumMetabolism` or `ConsortiumMetabolismAlignment` object.

## Value

A character vector containing the names of metabolites in the network.

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
metabolites(cm)
#> [1] "met5" "met3" "met4" "met1" "met2"
```
