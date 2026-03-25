# Get Unique Pathways

Returns pathways unique to query and reference in a pairwise alignment.

## Usage

``` r
uniquePathways(object)

# S4 method for class 'ConsortiumMetabolismAlignment'
uniquePathways(object)
```

## Arguments

- object:

  A
  [ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)
  object of type `"pairwise"`.

## Value

A list with `query` and `reference` data.frames.

## Methods (by class)

- `uniquePathways(ConsortiumMetabolismAlignment)`: Unique pathways from
  a pairwise
  [ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cma <- align(cm1, cm2)
uniquePathways(cma)
#> $query
#> # A tibble: 9 × 2
#>   consumed produced
#>   <chr>    <chr>   
#> 1 met1     met2    
#> 2 met4     met2    
#> 3 met5     met2    
#> 4 met1     met3    
#> 5 met3     met4    
#> 6 met1     met5    
#> 7 met2     met5    
#> 8 met3     met5    
#> 9 met4     met5    
#> 
#> $reference
#> # A tibble: 11 × 2
#>    consumed produced
#>    <chr>    <chr>   
#>  1 met2     met1    
#>  2 met5     met1    
#>  3 met6     met1    
#>  4 met6     met2    
#>  5 met2     met3    
#>  6 met6     met3    
#>  7 met6     met4    
#>  8 met6     met5    
#>  9 met2     met6    
#> 10 met4     met6    
#> 11 met5     met6    
#> 
```
