# Get the Consortia

Returns the consortia data in a tabular format.

## Usage

``` r
consortia(object)

# S4 method for class 'ConsortiumMetabolism'
consortia(object)

# S4 method for class 'ConsortiumMetabolismAlignment'
consortia(object)
```

## Arguments

- object:

  A `ConsortiumMetabolismAlignment` object.

## Value

For `ConsortiumMetabolism` objects, returns a tibble with species,
metabolite and flux information. For `ConsortiumMetabolismSet` and
`ConsortiumMetabolismAlignment` objects, returns a list of such tibbles.

A list with the community data.

## Methods (by class)

- `consortia(ConsortiumMetabolism)`: Get the Community

- `consortia(ConsortiumMetabolismAlignment)`: Not applicable for
  alignments

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
consortia(cm)
#>     met  species      flux
#> 1  met1 DUT1312V -1.178753
#> 2  met2 DUT1312V  1.231698
#> 3  met5 DUT1312V -1.686556
#> 4  met4 DUT1312V  3.573680
#> 5  met3 BSX9691E  3.006446
#> 6  met2 BSX9691E -0.970617
#> 7  met3 MSZ3715I -5.336016
#> 8  met5 MSZ3715I  2.592752
#> 9  met2 MSZ3715I -1.037214
#> 10 met1 MSZ3715I  5.514784
```
