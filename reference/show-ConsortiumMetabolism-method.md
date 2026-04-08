# Show Method for `ConsortiumMetabolism` Object

Show Method for `ConsortiumMetabolism` Object

## Usage

``` r
# S4 method for class 'ConsortiumMetabolism'
show(object)
```

## Arguments

- object:

  An object of class `ConsortiumMetabolism`

## Value

The object, invisibly.

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
show(cm)
#> 
#> ── ConsortiumMetabolism 
#> Name: "test"
#> Weighted metabolic network with 5 metabolites.
```
