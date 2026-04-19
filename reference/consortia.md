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
#>     met  species       flux
#> 1  met2 DUO5910U  2.5107857
#> 2  met2 CRP7583I -3.5003419
#> 3  met5 CRP7583I  3.1752836
#> 4  met5 TGH8437W -3.1752836
#> 5  met2 TGH8437W  0.9895562
#> 6  met1 DUO5910U  1.0452324
#> 7  met3 DUO5910U -0.8894377
#> 8  met1 CRP7583I  1.2633339
#> 9  met3 CRP7583I -1.2579714
#> 10 met4 CRP7583I  5.5036742
#> 11 met1 TGH8437W -2.3015127
#> 12 met3 TGH8437W  1.7197916
```
