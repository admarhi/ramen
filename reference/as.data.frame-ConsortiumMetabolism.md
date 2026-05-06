# Coerce a ConsortiumMetabolism to a data.frame

Returns the underlying edge list of a `ConsortiumMetabolism` as a plain
`data.frame` with three columns: `met`, `species`, and `flux`.
Equivalent to
`as.data.frame(object@InputData[, c("met", "species", "flux")])`.

## Usage

``` r
# S4 method for class 'ConsortiumMetabolism'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A
  [`ConsortiumMetabolism`](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
  object.

- row.names:

  Ignored.

- optional:

  Ignored.

- ...:

  Additional arguments (currently unused).

## Value

A `data.frame` with columns `met`, `species`, and `flux`.

## Examples

``` r
cm <- synCM("demo", n_species = 3, max_met = 5)
head(as.data.frame(cm))
#>    met  species       flux
#> 1 met2 DUO5910U  2.5107857
#> 2 met2 CRP7583I -3.5003419
#> 3 met5 CRP7583I  3.1752836
#> 4 met5 TGH8437W -3.1752836
#> 5 met2 TGH8437W  0.9895562
#> 6 met1 DUO5910U  1.0452324
```
