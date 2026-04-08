# Show method for `ConsortiumMetabolismAlignment` Objects

Show method for `ConsortiumMetabolismAlignment` Objects

## Usage

``` r
# S4 method for class 'ConsortiumMetabolismAlignment'
show(object)
```

## Arguments

- object:

  A `ConsortiumMetabolismAlignment` object.

## Value

The object, invisibly.

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cma <- align(cm1, cm2)
show(cma)
#> 
#> ── ConsortiumMetabolismAlignment 
#> Name: NA
#> Type: "pairwise"
#> Metric: "FOS"
#> Score: 0.8
#> Query: "comm_1", Reference: "comm_2"
#> Coverage: query 0.8, reference 0.25
```
