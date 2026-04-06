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
#> [1] 0.5714286
#> 
#> $jaccard
#> [1] 0.25
#> 
#> $brayCurtis
#> [1] 0.3327621
#> 
#> $redundancyOverlap
#> [1] 0.2272727
#> 
#> $coverageQuery
#> [1] 0.3076923
#> 
#> $coverageReference
#> [1] 0.5714286
#> 
```
