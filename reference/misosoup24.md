# Example Data from MiSoSoup 2024

List of 56 solutions from MiSoSoup for different focal strains in
different media.

## Format

A list of 56 tibbles, each with columns:

- species:

  Character, species identifier.

- metabolite:

  Character, metabolite name (BiGG identifiers).

- flux:

  Numeric, metabolic flux (negative for consumption, positive for
  production).

## Value

A list of 56 tibbles.

## Examples

``` r
data("misosoup24")
head(misosoup24[[1]])
#> # A tibble: 6 × 3
#>   metabolite species    flux
#>   <chr>      <chr>     <dbl>
#> 1 ac         A1R12     0.773
#> 2 ac         I2R16   -10.8  
#> 3 acald      A1R12    -1.12 
#> 4 acald      I2R16     1.12 
#> 5 ala__D     A1R12     0.760
#> 6 ala__D     I2R16    -0.760
```
