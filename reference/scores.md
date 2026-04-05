# Get Alignment Scores

Returns the scores from a
[ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)
object.

## Usage

``` r
scores(object)

# S4 method for class 'ConsortiumMetabolismAlignment'
scores(object)
```

## Arguments

- object:

  A
  [ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)
  object.

## Value

A named list of scores.

## Methods (by class)

- `scores(ConsortiumMetabolismAlignment)`: Scores from a
  [ConsortiumMetabolismAlignment](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismAlignment.md)

## Examples

``` r
cm1 <- synCM("comm_1", n_species = 3, max_met = 5)
cm2 <- synCM("comm_2", n_species = 4, max_met = 6)
cma <- align(cm1, cm2)
scores(cma)
#> $FOS
#> [1] 0.1666667
#> 
#> $jaccard
#> [1] 0.07692308
#> 
#> $brayCurtis
#> [1] 0.05323122
#> 
#> $redundancyOverlap
#> [1] 0.05882353
#> 
```
