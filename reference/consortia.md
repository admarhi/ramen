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
#>     met  species        flux
#> 1  met1 DRQ2970J -2.01564102
#> 2  met4 DRQ2970J -0.09895164
#> 3  met2 DRQ2970J -3.54985640
#> 4  met5 DRQ2970J  3.27273818
#> 5  met3 DRQ2970J -0.08681708
#> 6  met5 GNO3319G -3.46695470
#> 7  met1 GNO3319G  1.26745490
#> 8  met3 GNO3319G -5.29857257
#> 9  met2 GNO3319G -1.42685091
#> 10 met4 RVV5814P -7.76167419
#> 11 met1 RVV5814P  3.36388087
```
