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
#> 1 met3  AUP039Y -0.4343764
#> 2 met3 KYS6207S -0.9798091
#> 3 met5 KYS6207S  0.9716773
#> 4 met5  QZV739R -0.9716773
#> 5 met3  QZV739R  1.4141855
```
