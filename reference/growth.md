# Get Growth Rates

Returns per-species growth rates (e.g. FBA objective values) stored in a
ramen object. For
[ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
objects, returns the named numeric vector supplied at construction. For
[ConsortiumMetabolismSet](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md)
objects, returns a named list with one entry per consortium.

## Usage

``` r
growth(object)

# S4 method for class 'ConsortiumMetabolism'
growth(object)

# S4 method for class 'ConsortiumMetabolismSet'
growth(object)
```

## Arguments

- object:

  A
  [ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md)
  or
  [ConsortiumMetabolismSet](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md)
  object.

## Value

For
[ConsortiumMetabolism](https://admarhi.github.io/ramen/reference/ConsortiumMetabolism.md):
a named numeric vector of growth rates (names = species), or `NULL` if
no growth data was supplied. For
[ConsortiumMetabolismSet](https://admarhi.github.io/ramen/reference/ConsortiumMetabolismSet.md):
a named list of such vectors.

## Methods (by class)

- `growth(ConsortiumMetabolismSet)`: Get Growth Rates From a
  `ConsortiumMetabolismSet`

## Examples

``` r
data <- data.frame(
    species = c("s1", "s1", "s2", "s2"),
    metabolite = c("m1", "m2", "m1", "m3"),
    flux = c(-1, 1, -1, 1)
)
cm <- ConsortiumMetabolism(
    data, name = "test",
    growth = c(s1 = 0.5, s2 = 0.3)
)
growth(cm)
#>  s1  s2 
#> 0.5 0.3 
```
