# Get the Consortia

Returns the consortia data in a tabular format.

## Usage

``` r
consortia(object)

# S4 method for class 'ConsortiumMetabolism'
consortia(object)

# S4 method for class 'ConsortiumMetabolismAlignment'
consortia(object)

# S4 method for class 'ConsortiumMetabolismSet'
consortia(object)
```

## Arguments

- object:

  A `ConsortiumMetabolismSet` object.

## Value

For `ConsortiumMetabolism` objects, returns a tibble with species,
metabolite and flux information. For `ConsortiumMetabolismSet` and
`ConsortiumMetabolismAlignment` objects, returns a list of such tibbles.

A list with the community data.

A named list of `ConsortiumMetabolism` objects.

## Methods (by class)

- `consortia(ConsortiumMetabolism)`: Get the Community

- `consortia(ConsortiumMetabolismAlignment)`: Not applicable for
  alignments

- `consortia(ConsortiumMetabolismSet)`: Get the list of
  `ConsortiumMetabolism` objects from a `ConsortiumMetabolismSet`

## Examples

``` r
cm <- synCM("test", n_species = 3, max_met = 5)
consortia(cm)
#>    met  species       flux
#> 1 met2 QWG2904R -2.3108858
#> 2 met2  GRW821Q -1.8850686
#> 3 met5  GRW821Q  8.7532211
#> 4 met5 GIH4406S -2.7387147
#> 5 met2 GIH4406S  4.1959544
#> 6 met3 QWG2904R  0.4356225
#> 7 met1 QWG2904R  0.9238720
#> 8 met5 QWG2904R -5.8595975
```
