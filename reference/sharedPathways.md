# Get Shared Pathways

Returns pathways shared between query and reference in a pairwise
alignment.

## Usage

``` r
sharedPathways(object)

# S4 method for class 'ConsortiumMetabolismAlignment'
sharedPathways(object)
```

## Arguments

- object:

  A
  [ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)
  object of type `"pairwise"`.

## Value

A data.frame of shared pathways.

## Methods (by class)

- `sharedPathways(ConsortiumMetabolismAlignment)`: Shared pathways from
  a pairwise
  [ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cma <- align(cm1, cm2)
sharedPathways(cma)
#> # A tibble: 2 × 4
#>   consumed produced querySpecies     referenceSpecies
#>   <chr>    <chr>    <list>           <list>          
#> 1 met1     met3     <tibble [1 × 3]> <tibble [1 × 3]>
#> 2 met2     met4     <tibble [1 × 3]> <tibble [1 × 3]>
```
