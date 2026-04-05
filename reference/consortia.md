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
#>     met  species       flux
#> 1  met2 EMC1269Y  2.0882296
#> 2  met2 JCA8468A -2.0882296
#> 3  met3 JCA8468A  0.9595214
#> 4  met3 HRC1916Z -0.9595214
#> 5  met2 HRC1916Z  0.6600125
#> 6  met2 EMC1269Y -0.6600125
#> 7  met4 JCA8468A  0.4940870
#> 8  met5 HRC1916Z  1.6160980
#> 9  met1 HRC1916Z  7.9723679
#> 10 met4 HRC1916Z  3.2459186
```
